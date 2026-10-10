//! The modular pipeline, a faster reorganization of the HKTY procedure (work in progress).
//!
//! It computes the invariants modulo a few word-sized primes and recovers them by the Chinese
//! remainder theorem, stores curve classes as packed integer keys, and never enumerates the
//! whole semigroup. The design follows the C program in PR #94 and is described in `PLAN.md`.
//! So far this module holds only the primitives the stages are built from: the arithmetic
//! ([`field`]), the curve-class keys ([`key`]) and the tables keyed by them ([`table`]).

pub mod field;
pub mod key;
pub mod table;

/// A small deterministic pseudorandom generator for the unit tests (SplitMix64).
#[cfg(test)]
pub(crate) struct TestRng(u64);

#[cfg(test)]
impl TestRng {
    pub(crate) fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A uniformly distributed integer in `lo..=hi`.
    pub(crate) fn range(&mut self, lo: i64, hi: i64) -> i64 {
        let span = (hi - lo) as u64 + 1;
        lo + (self.next_u64() % span) as i64
    }
}
