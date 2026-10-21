//! Curve classes as fixed-size keys.
//!
//! The modular pipeline never stores a curve class as a vector. It packs the coordinates into
//! a [`Key`], a small `Copy` value with $\text{key}(a + b) = \text{key}(a) + \text{key}(b)$, so
//! that multiplying monomials is one addition, and unpacks it only where the coordinates are
//! needed (the period step, the multicover and the output). The encoding is lossless, not a
//! hash: [`Key::hash64`] picks table slots separately. A coordinate that does not fit its
//! field would silently corrupt the key, so a [`KeyLayout`] is sized up front from a bound
//! that holds for every class the run can produce ([`coordinate_bounds`]).
//!
//! How the coordinates are packed is chosen per run, and the tables are generic over it. There
//! is one layout, [`PackedLayout`], over integers of 128, 256 or 512 bits ([`Word`]), and the
//! run takes the narrowest whose fields fit. (Keys of 64 bits, which most of the benchmark
//! models fit, made no measurable difference to their time or memory, so they are not worth
//! another copy of the pipeline.)

use super::level::LevelWord;
use core::fmt::{self, Debug};
use core::ops::Add;

/// A packed curve class.
///
/// [`Default`] must give the zero class, and addition must add the classes. The order is
/// arbitrary but total, and only serves to make sorted output deterministic.
pub trait Key:
    Copy + Eq + Ord + core::hash::Hash + Default + Add<Output = Self> + Send + Sync + Debug
{
    /// A hash of the key, for picking a slot in a table.
    fn hash64(self) -> u64;

    /// The key of `n` times the class. Keys add like the classes, so this is `n` times the
    /// key; the result is only meaningful when the class fits the layout.
    fn times(self, n: i64) -> Self;

    /// The widest word that the extraction's level tables pack offsets into for these keys:
    /// `u128`, or a wider one for keys wider than that.
    type Level: LevelWord;
}

/// How curve classes are packed into keys.
pub trait KeyLayout: Send + Sync {
    /// The key type.
    type Key: Key;

    /// The number of coordinates of a curve class.
    fn dim(&self) -> usize;

    /// Whether every coordinate of `coords` fits its field.
    fn fits(&self, coords: &[i32]) -> bool;

    /// Packs a curve class. The coordinates must fit (see [`KeyLayout::fits`]), which is only
    /// checked in debug builds.
    fn pack(&self, coords: &[i32]) -> Self::Key;

    /// Packs a curve class, or fails if a coordinate does not fit.
    fn try_pack(&self, coords: &[i32]) -> Result<Self::Key, KeyError> {
        if self.fits(coords) {
            Ok(self.pack(coords))
        } else {
            Err(KeyError::CoordinateOutOfRange)
        }
    }

    /// Unpacks a curve class into `coords`, which must have length [`KeyLayout::dim`].
    fn unpack(&self, key: Self::Key, coords: &mut [i32]);

    /// A key that no curve class packs to, for marking empty table entries.
    fn sentinel(&self) -> Self::Key;
}

/// Errors in setting up or using a key layout.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KeyError {
    /// The generators, grading vector and bounds do not all have the same dimension.
    DimensionMismatch,
    /// The grading vector is not positive on a nonzero generator.
    NonPositiveDegree,
    /// A coordinate bound exceeds the range of an `i32`.
    BoundTooLarge,
    /// The fields need more bits than the key has.
    TooWide { bits: u32, max_bits: u32 },
    /// A coordinate does not fit its field.
    CoordinateOutOfRange,
}

impl fmt::Display for KeyError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::DimensionMismatch => write!(f, "the dimensions of the key data do not match"),
            Self::NonPositiveDegree => {
                write!(
                    f,
                    "the grading vector must be positive on every nonzero generator"
                )
            }
            Self::BoundTooLarge => write!(f, "a coordinate bound does not fit in an i32"),
            Self::TooWide { bits, max_bits } => write!(
                f,
                "the curve classes need {bits} bits, but keys have at most {max_bits}"
            ),
            Self::CoordinateOutOfRange => write!(f, "a coordinate does not fit its key field"),
        }
    }
}

impl std::error::Error for KeyError {}

/// Bounds on the absolute value of each coordinate of every class in the cone spanned by
/// `generators` with degree at most `max_deg`.
///
/// A point $x = \sum_i \lambda_i g_i$ with $\lambda_i \ge 0$ has
/// $\lvert x_t \rvert \le \sum_i \lambda_i \lvert g_{it} \rvert
/// \le \deg(x) \max_i \lvert g_{it} \rvert / \deg(g_i)$, and $x_t$ is an integer, so with
/// $D$ = `max_deg` the bound is $\max_i \lfloor D \lvert g_{it} \rvert / \deg(g_i) \rfloor$. The
/// generators only need to span the cone, so its rays are enough, and zero generators are
/// ignored. Every class that the pipeline forms is a sum of classes in the cone with degree
/// at most `max_deg`, and is therefore covered.
pub fn coordinate_bounds<'a>(
    generators: impl IntoIterator<Item = &'a [i32]>,
    grading_vector: &[i32],
    max_deg: u32,
) -> Result<Vec<u32>, KeyError> {
    let dim = grading_vector.len();
    let mut bounds = vec![0u64; dim];
    for g in generators {
        if g.len() != dim {
            return Err(KeyError::DimensionMismatch);
        }
        if g.iter().all(|&x| x == 0) {
            continue;
        }
        let deg: i64 = g
            .iter()
            .zip(grading_vector)
            .map(|(&x, &w)| x as i64 * w as i64)
            .sum();
        if deg <= 0 {
            return Err(KeyError::NonPositiveDegree);
        }
        for (b, &x) in bounds.iter_mut().zip(g) {
            *b = (*b).max(max_deg as u64 * x.unsigned_abs() as u64 / deg as u64);
        }
    }
    bounds
        .into_iter()
        .map(|b| u32::try_from(b).ok().filter(|&b| b <= i32::MAX as u32))
        .collect::<Option<_>>()
        .ok_or(KeyError::BoundTooLarge)
}

/// The integer that a [`PackedLayout`] packs curve classes into: `BITS` bits of two's
/// complement with wrapping arithmetic. Implemented by `u128` and [`Wide`].
pub trait Word:
    Copy + Eq + Ord + Default + core::hash::Hash + Send + Sync + Debug + 'static
{
    /// The number of bits.
    const BITS: u32;

    /// `x`, sign-extended.
    fn from_i64(x: i64) -> Self;

    fn wrapping_add(self, other: Self) -> Self;

    fn wrapping_sub(self, other: Self) -> Self;

    /// `n` times `self`.
    fn wrapping_mul_i64(self, n: i64) -> Self;

    /// Shifts left by `b < 64` bits.
    fn shl(self, b: u32) -> Self;

    /// Shifts right by `b < 64` bits, filling with zeros.
    fn shr(self, b: u32) -> Self;

    /// The low 64 bits.
    fn low(self) -> u64;

    /// $2^{\text{BITS} - 1}$.
    fn top() -> Self;

    /// A hash, for picking a slot in a table.
    fn hash64(self) -> u64;

    /// See [`Key::Level`].
    type Level: LevelWord;
}

/// The finalizer of the hashes, after splitmix64.
#[inline]
fn mix(mut x: u64) -> u64 {
    x ^= x >> 31;
    x = x.wrapping_mul(0xD6E8_FEB8_6659_FD93);
    x ^ (x >> 32)
}

impl Word for u128 {
    const BITS: u32 = u128::BITS;
    type Level = u128;

    #[inline]
    fn from_i64(x: i64) -> Self {
        x as i128 as u128
    }

    #[inline]
    fn wrapping_add(self, other: Self) -> Self {
        u128::wrapping_add(self, other)
    }

    #[inline]
    fn wrapping_sub(self, other: Self) -> Self {
        u128::wrapping_sub(self, other)
    }

    #[inline]
    fn wrapping_mul_i64(self, n: i64) -> Self {
        u128::wrapping_mul(self, n as i128 as u128)
    }

    #[inline]
    fn shl(self, b: u32) -> Self {
        self << b
    }

    #[inline]
    fn shr(self, b: u32) -> Self {
        self >> b
    }

    #[inline]
    fn low(self) -> u64 {
        self as u64
    }

    #[inline]
    fn top() -> Self {
        1 << 127
    }

    #[inline]
    fn hash64(self) -> u64 {
        // The same mix as hash128 in the C program.
        mix((self as u64) ^ ((self >> 64) as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15))
    }
}

/// An integer of `64 * N` bits, as `N` words from the least significant one, for the keys
/// that `u128` cannot hold (256 and 512 bits, as in the C program of PR #102).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Wide<const N: usize>(pub [u64; N]);

impl<const N: usize> Default for Wide<N> {
    #[inline]
    fn default() -> Self {
        Self([0; N])
    }
}

impl<const N: usize> Word for Wide<N> {
    const BITS: u32 = 64 * N as u32;
    type Level = Self;

    #[inline]
    fn from_i64(x: i64) -> Self {
        let mut w = [if x < 0 { u64::MAX } else { 0 }; N];
        w[0] = x as u64;
        Self(w)
    }

    #[inline]
    fn wrapping_add(self, other: Self) -> Self {
        let mut w = [0; N];
        let mut carry = false;
        for ((x, &a), &b) in w.iter_mut().zip(&self.0).zip(&other.0) {
            let (s, c1) = a.overflowing_add(b);
            let (s, c2) = s.overflowing_add(carry as u64);
            *x = s;
            carry = c1 || c2;
        }
        Self(w)
    }

    #[inline]
    fn wrapping_sub(self, other: Self) -> Self {
        let mut w = [0; N];
        let mut borrow = false;
        for ((x, &a), &b) in w.iter_mut().zip(&self.0).zip(&other.0) {
            let (d, b1) = a.overflowing_sub(b);
            let (d, b2) = d.overflowing_sub(borrow as u64);
            *x = d;
            borrow = b1 || b2;
        }
        Self(w)
    }

    #[inline]
    fn wrapping_mul_i64(self, n: i64) -> Self {
        let m = n.unsigned_abs() as u128;
        let mut w = [0; N];
        let mut carry = 0u128;
        for (x, &a) in w.iter_mut().zip(&self.0) {
            let p = a as u128 * m + carry;
            *x = p as u64;
            carry = p >> 64;
        }
        let w = Self(w);
        if n < 0 {
            Self::default().wrapping_sub(w)
        } else {
            w
        }
    }

    #[inline]
    fn shl(self, b: u32) -> Self {
        debug_assert!(b < 64);
        if b == 0 {
            return self;
        }
        let mut w = [0; N];
        for i in (0..N).rev() {
            w[i] = self.0[i] << b;
            if i > 0 {
                w[i] |= self.0[i - 1] >> (64 - b);
            }
        }
        Self(w)
    }

    #[inline]
    fn shr(self, b: u32) -> Self {
        debug_assert!(b < 64);
        if b == 0 {
            return self;
        }
        let mut w = [0; N];
        for (i, x) in w.iter_mut().enumerate() {
            *x = self.0[i] >> b;
            if i + 1 < N {
                *x |= self.0[i + 1] << (64 - b);
            }
        }
        Self(w)
    }

    #[inline]
    fn low(self) -> u64 {
        self.0[0]
    }

    #[inline]
    fn top() -> Self {
        let mut w = [0; N];
        w[N - 1] = 1 << 63;
        Self(w)
    }

    #[inline]
    fn hash64(self) -> u64 {
        let mut h = 0u64;
        for &x in &self.0 {
            h = (h ^ x).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(29);
        }
        mix(h)
    }
}

/// A curve class packed into a [`Word`] by a [`PackedLayout`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
#[repr(transparent)]
pub struct PackedKey<W = u128>(pub W);

impl<W: Word> Add for PackedKey<W> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self(self.0.wrapping_add(other.0))
    }
}

impl<W: Word> Key for PackedKey<W> {
    type Level = W::Level;

    #[inline]
    fn hash64(self) -> u64 {
        self.0.hash64()
    }

    #[inline]
    fn times(self, n: i64) -> Self {
        Self(self.0.wrapping_mul_i64(n))
    }
}

/// The widths of the fields that hold the coordinates in `-bound..=bound`, each with its own
/// bound: a signed field of the bit length of the bound plus one, and none for a coordinate
/// that is always zero.
pub fn field_widths(bounds: &[u32]) -> Vec<u32> {
    bounds
        .iter()
        .map(|&b| {
            if b == 0 {
                0
            } else {
                u32::BITS - b.leading_zeros() + 1
            }
        })
        .collect()
}

/// Packs the coordinates into one [`Word`] `W`, each in a signed field of its own width.
///
/// Coordinate $t$ occupies bits $o_t$ to $o_t + b_t - 1$ with $o_t = \sum_{s < t} b_s$, and the
/// key is $\sum_t x_t 2^{o_t} \bmod 2^{n}$, $n$ the bits of `W`, i.e. two's complement
/// throughout. That is linear, so keys add like the classes, with carries and borrows between
/// the fields that unpacking undoes. A field of $b$ bits holds $-2^{b-1} \le x_t < 2^{b-1}$,
/// and a coordinate that is always zero needs no bits at all.
///
/// The widths are sized from bounds on the coordinates ([`coordinate_bounds`],
/// [`field_widths`]) rather than split evenly as in the C program, which lets more
/// coordinates fit. They may total at most [`PackedLayout::MAX_BITS`], $n - 2$: then every
/// key lies strictly between $-2^{n-2}$ and $2^{n-2}$ as a signed integer, and $-2^{n-1}$ is
/// free to serve as the sentinel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedLayout<W = u128> {
    widths: Vec<u32>,
    word: core::marker::PhantomData<W>,
}

impl<W: Word> PackedLayout<W> {
    /// The most bits the fields may use together.
    pub const MAX_BITS: u32 = W::BITS - 2;

    /// A layout whose fields hold the coordinates in `-bound..=bound`, each with its own bound.
    pub fn new(bounds: &[u32]) -> Result<Self, KeyError> {
        let widths = field_widths(bounds);
        let bits = widths.iter().sum();
        if bits > Self::MAX_BITS {
            return Err(KeyError::TooWide {
                bits,
                max_bits: Self::MAX_BITS,
            });
        }
        Ok(Self {
            widths,
            word: core::marker::PhantomData,
        })
    }

    /// The width of each field in bits.
    pub fn widths(&self) -> &[u32] {
        &self.widths
    }

    /// The width of all fields together in bits.
    pub fn bits(&self) -> u32 {
        self.widths.iter().sum()
    }
}

impl<W: Word> KeyLayout for PackedLayout<W> {
    type Key = PackedKey<W>;

    fn dim(&self) -> usize {
        self.widths.len()
    }

    fn fits(&self, coords: &[i32]) -> bool {
        coords.len() == self.widths.len()
            && coords.iter().zip(&self.widths).all(|(&x, &b)| match b {
                0 => x == 0,
                _ => {
                    let half = 1i64 << (b - 1);
                    (-half..half).contains(&(x as i64))
                }
            })
    }

    #[inline]
    fn pack(&self, coords: &[i32]) -> PackedKey<W> {
        debug_assert!(self.fits(coords), "{coords:?} does not fit {self:?}");
        let mut k = W::default();
        for (&x, &b) in coords.iter().zip(&self.widths).rev() {
            k = k.shl(b).wrapping_add(W::from_i64(x as i64));
        }
        PackedKey(k)
    }

    #[inline]
    fn unpack(&self, key: PackedKey<W>, coords: &mut [i32]) {
        debug_assert_eq!(coords.len(), self.widths.len());
        let mut k = key.0;
        for (x, &b) in coords.iter_mut().zip(&self.widths) {
            if b == 0 {
                *x = 0;
                continue;
            }
            // Sign-extend the low b bits (b <= 32), then remove them, borrow included.
            let f = ((k.low() as i64) << (64 - b)) >> (64 - b);
            *x = f as i32;
            k = k.wrapping_sub(W::from_i64(f)).shr(b);
        }
    }

    fn sentinel(&self) -> PackedKey<W> {
        PackedKey(W::top())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::TestRng;

    fn random_point<W: Word>(rng: &mut TestRng, layout: &PackedLayout<W>) -> Vec<i32> {
        layout
            .widths()
            .iter()
            .map(|&b| match b {
                0 => 0,
                _ => rng.range(-(1 << (b - 1)), (1 << (b - 1)) - 1) as i32,
            })
            .collect()
    }

    /// Layouts of every word they fit.
    fn layouts<W: Word>() -> Vec<PackedLayout<W>> {
        [
            vec![1, 1],
            vec![5, 0, 3],
            vec![60; 8],
            vec![60; 11],
            vec![0, 0, 0],
            vec![i32::MAX as u32, i32::MAX as u32, 1000, 7],
            vec![2; 42],
            vec![1; 63],
            vec![1000; 20],
            vec![1; 127],
            vec![i32::MAX as u32; 7],
            vec![7; 127],
            vec![1; 255],
        ]
        .iter()
        .filter_map(|b| PackedLayout::new(b).ok())
        .collect()
    }

    #[test]
    fn test_widths() {
        let layout =
            PackedLayout::<u128>::new(&[0, 1, 2, 3, 4, 7, 8, 60, i32::MAX as u32]).unwrap();
        assert_eq!(layout.widths(), [0, 2, 3, 3, 4, 4, 5, 7, 32]);
        assert_eq!(layout.bits(), 60);
        // The bounds of 11 coordinates of 60, as at h11 = 11 and high degree, fit easily
        // where the C's even split of 11 bits each would not.
        assert_eq!(PackedLayout::<u128>::new(&[60; 11]).unwrap().bits(), 77);
        assert_eq!(PackedLayout::<u128>::new(&[1; 63]).unwrap().bits(), 126);
        assert_eq!(
            PackedLayout::<u128>::new(&[1; 64]),
            Err(KeyError::TooWide {
                bits: 128,
                max_bits: 126
            })
        );
        assert_eq!(PackedLayout::<Wide<4>>::new(&[1; 127]).unwrap().bits(), 254);
        assert!(PackedLayout::<Wide<4>>::new(&[1; 128]).is_err());
        assert_eq!(PackedLayout::<Wide<8>>::new(&[1; 255]).unwrap().bits(), 510);
        assert!(PackedLayout::<Wide<8>>::new(&[1; 256]).is_err());
    }

    /// Random words of `Wide<2>` and `u128` and their sums, differences, products and shifts
    /// agree.
    #[test]
    fn test_wide_arithmetic() {
        let mut rng = TestRng::new(9);
        let wide = |x: u128| Wide::<2>([x as u64, (x >> 64) as u64]);
        let random = |rng: &mut TestRng| match rng.range(0, 3) {
            0 => rng.next_u64() as u128,
            1 => (rng.range(-1000, 1000) as i128) as u128,
            _ => (rng.next_u64() as u128) << 64 | rng.next_u64() as u128,
        };
        for _ in 0..10000 {
            let (a, b) = (random(&mut rng), random(&mut rng));
            let n = match rng.range(0, 2) {
                0 => rng.range(-5, 5),
                1 => rng.next_u64() as i64,
                _ => i64::MIN,
            };
            let s = rng.range(0, 63) as u32;
            assert_eq!(wide(a).wrapping_add(wide(b)), wide(a.wrapping_add(b)));
            assert_eq!(wide(a).wrapping_sub(wide(b)), wide(a.wrapping_sub(b)));
            assert_eq!(
                wide(a).wrapping_mul_i64(n),
                wide(a.wrapping_mul_i64(n)),
                "{a} {n}"
            );
            assert_eq!(wide(a).shl(s), wide(a << s));
            assert_eq!(wide(a).shr(s), wide(a >> s));
            assert_eq!(
                <Wide<2> as Word>::from_i64(n),
                wide(<u128 as Word>::from_i64(n))
            );
            assert_eq!(wide(a).low(), a as u64);
        }
        assert_eq!(<Wide<2> as Word>::top(), wide(<u128 as Word>::top()));
    }

    fn round_trip_and_additivity<W: Word>(seed: u64) {
        let mut rng = TestRng::new(seed);
        for layout in layouts::<W>() {
            let n = layout.dim();
            let mut out = vec![0; n];
            for _ in 0..2000 {
                let a = random_point(&mut rng, &layout);
                let ka = layout.pack(&a);
                layout.unpack(ka, &mut out);
                assert_eq!(out, a, "{layout:?}");
                assert_ne!(ka, layout.sentinel());

                // Splitting a point into two summands that fit, packing them separately and
                // adding the keys gives back the key of the point.
                let b: Vec<i32> = a
                    .iter()
                    .zip(layout.widths())
                    .map(|(&x, &w)| {
                        if w == 0 {
                            return 0;
                        }
                        let half = 1i64 << (w - 1);
                        let lo = (-half).max(x as i64 - half + 1);
                        let hi = (half - 1).min(x as i64 + half);
                        rng.range(lo, hi) as i32
                    })
                    .collect();
                let c: Vec<i32> = a.iter().zip(&b).map(|(x, y)| x - y).collect();
                assert!(layout.fits(&b) && layout.fits(&c), "{a:?} {b:?} {c:?}");
                assert_eq!(
                    layout.pack(&b) + layout.pack(&c),
                    ka,
                    "{a:?} = {b:?} + {c:?}"
                );

                // A multiple that fits.
                let m = rng.range(-3, 3);
                let ma: Option<Vec<i32>> = a
                    .iter()
                    .map(|&x| i32::try_from(x as i64 * m).ok())
                    .collect();
                if let Some(ma) = ma.filter(|ma| layout.fits(ma)) {
                    assert_eq!(ka.times(m), layout.pack(&ma));
                }
            }

            // The extremes of every field at once.
            for sign in [-1, 1] {
                let a: Vec<i32> = layout
                    .widths()
                    .iter()
                    .map(|&b| match (b, sign) {
                        (0, _) => 0,
                        (_, -1) => (-(1i64 << (b - 1))) as i32,
                        _ => ((1i64 << (b - 1)) - 1) as i32,
                    })
                    .collect();
                let k = layout.pack(&a);
                layout.unpack(k, &mut out);
                assert_eq!(out, a);
                assert_ne!(k, layout.sentinel());
            }
            assert_eq!(layout.pack(&vec![0; n]), PackedKey::default());
        }
    }

    #[test]
    fn test_round_trip_and_additivity() {
        assert_eq!(layouts::<u128>().len(), 8);
        assert_eq!(layouts::<Wide<4>>().len(), 11);
        assert_eq!(layouts::<Wide<8>>().len(), 13);
        round_trip_and_additivity::<u128>(5);
        round_trip_and_additivity::<Wide<4>>(6);
        round_trip_and_additivity::<Wide<8>>(7);
    }

    #[test]
    fn test_overflow_detection() {
        let layout = PackedLayout::<u128>::new(&[3, 0, 100]).unwrap();
        assert!(layout.fits(&[3, 0, 100]));
        assert!(layout.fits(&[-4, 0, -128]));
        assert!(!layout.fits(&[4, 0, 0]));
        assert!(!layout.fits(&[-5, 0, 0]));
        assert!(!layout.fits(&[0, 1, 0]));
        assert!(!layout.fits(&[0, 0, 128]));
        assert!(!layout.fits(&[0, 0]));
        assert_eq!(layout.try_pack(&[1, 0, -1]), Ok(layout.pack(&[1, 0, -1])));
        assert_eq!(
            layout.try_pack(&[0, 0, 200]),
            Err(KeyError::CoordinateOutOfRange)
        );
    }

    #[test]
    fn test_coordinate_bounds() {
        // The generators of the threefold in tests/full_hkty.rs, both of degree 1.
        let gens = [[0, -1], [1, 2], [0, 0]];
        let w = [3, -1];
        let bounds = coordinate_bounds(gens.iter().map(|g| &g[..]), &w, 10).unwrap();
        assert_eq!(bounds, [10, 20]);

        // Every point of the generated semigroup up to the degree is within the bounds, and
        // they are attained.
        let mut max = [0i32; 2];
        for a in 0..=10 {
            for b in 0..=10 - a {
                let x: [i32; 2] = [b, -a + 2 * b];
                for t in 0..2 {
                    assert!(x[t].unsigned_abs() <= bounds[t]);
                    max[t] = max[t].max(x[t].abs());
                }
            }
        }
        assert_eq!(max, [10, 20]);

        // A degree that does not divide evenly rounds down.
        let gens = [[3, 1], [0, 1]];
        let bounds = coordinate_bounds(gens.iter().map(|g| &g[..]), &[1, 2], 7).unwrap();
        assert_eq!(bounds, [4, 3]);

        assert_eq!(
            coordinate_bounds([&[1, -1][..]], &[1, 1], 5),
            Err(KeyError::NonPositiveDegree)
        );
        assert_eq!(
            coordinate_bounds([&[1, 0, 0][..]], &[1, 1], 5),
            Err(KeyError::DimensionMismatch)
        );
        assert_eq!(
            coordinate_bounds([&[1, i32::MAX][..]], &[1, 0], 2),
            Err(KeyError::BoundTooLarge)
        );
    }
}
