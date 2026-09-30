"""cgv in any phase of a toric hypersurface: FRST or vex ("BGV").

A vex phase is a fine regular fan on the boundary points (not interior to facets) of a reflexive 4-polytope that does
not refine its face fan (MacFadden-Sheridan, arXiv:2512.14817). Curves of negative anticanonical degree m = -K.C < 0
exist there and HKTY's coefficient has a pole; it is removed exactly by the class of the stratum V_S (S = the divisors
the curve meets negatively), which lies inside the CY (Liam McAllister's note, via J. Wang's quasimap I-function,
arXiv:1910.14440). For fine fans on reflexive polytopes S is always a 3-cone (MacFadden-Sheridan Prop. 5), so each
term needs only the class of a toric curve. cgv computes these terms when the input carries the strata section that
this module writes; without it cgv is unchanged.

    from cgv_phase import compute_gv
    gv = compute_gv(cones, q, kappa, max_deg)           # {curve tuple: int GV}

Inputs (cygv's, plus the fan):
    cones      maximal cones of the fan, as tuples of column indices of q (0-based)
    q          GLSM charges, h11 x N (column i: D_i = sum_a q[a][i] J_a)
    kappa      CY triple intersections in the same basis: {(a, b, c): v} or [[a, b, c, v], ...]
optional:
    generators Mori cone generators, or those of a subcone (lightcone GVs); default: the fan's wall curves
    saturate   True: every lattice point of the cone they span (normaliz Hilbert basis); False: exactly their
               semigroup (as cygv), full enumeration
    grading    grading vector, positive on the generators; default: an interior point of the dual cone

Derived here (numpy + normaliz; no CYTools): ray vectors (integer Gale dual of q), the fan's toric intersection ring
(kappa is checked against it), wall curves, vex cones (rays sharing no facet of conv(rays)) and their exact classes.
phase_from_polytope() builds (cones, q, kappa) from a CYTools polytope and fan.
"""
import itertools, os, subprocess, sys, tempfile
from fractions import Fraction
from functools import lru_cache
from math import gcd

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if os.path.exists(os.path.join(HERE, "cgv_run.py")):   # standalone: cgv/tools
    sys.path.insert(0, HERE)
    import cgv_run  # noqa: E402  (normaliz helpers and cgv's input writer)
else:                                                   # bundled in cygv, next to cygv._cgv_run
    from cygv import _cgv_run as cgv_run  # noqa: E402

BIN = os.path.join(HERE, "..", "cgv")


# ---------------------------------------------------------------- exact linear algebra
def _solve(rows, rhs):
    """exact solution of a consistent (possibly overdetermined) rational system with a unique solution"""
    n = len(rows[0]); A = [[Fraction(x) for x in r] + [Fraction(b)] for r, b in zip(rows, rhs)]
    piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]; A[r] = [x / A[r][c] for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c); r += 1
    if any(A[i][n] != 0 for i in range(r, len(A))): raise ValueError("inconsistent linear system")
    if len(piv) != n: raise ValueError("underdetermined linear system")
    x = [Fraction(0)] * n
    for i, c in enumerate(piv): x[c] = A[i][n]
    return x


def _particular(rows, rhs, n):
    """some rational solution of an underdetermined consistent system (free variables = 0)"""
    A = [[Fraction(x) for x in r] + [Fraction(b)] for r, b in zip(rows, rhs)]
    piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]; A[r] = [x / A[r][c] for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c); r += 1
    x = [Fraction(0)] * n
    for i, c in enumerate(piv): x[c] = A[i][n]
    return x


def _det(M):
    M = [[Fraction(int(x)) for x in r] for r in M]; n = len(M); d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return 0
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]; M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return int(d)


def gale_dual(q):
    """Z-basis of the integer kernel of q (h11 x N): N x (N - rank) matrix whose rows are the rays, up to GL(4, Z)"""
    q = [[int(x) for x in r] for r in q]; h, N = len(q), len(q[0])
    A = [row[:] for row in q]
    U = [[int(i == j) for j in range(N)] for i in range(N)]     # column operations: A = q U
    r = 0
    for i in range(h):
        # reduce row i on columns r.. to a single nonzero entry at column r (extended gcd by column ops)
        while True:
            nz = [c for c in range(r, N) if A[i][c] != 0]
            if len(nz) <= 1: break
            c0 = min(nz, key=lambda c: abs(A[i][c]))
            for c in nz:
                if c == c0: continue
                f = A[i][c] // A[i][c0]
                for rr in range(h): A[rr][c] -= f * A[rr][c0]
                for rr in range(N): U[rr][c] -= f * U[rr][c0]
        nz = [c for c in range(r, N) if A[i][c] != 0]
        if not nz: continue
        c0 = nz[0]
        for rr in range(h): A[rr][r], A[rr][c0] = A[rr][c0], A[rr][r]
        for rr in range(N): U[rr][r], U[rr][c0] = U[rr][c0], U[rr][r]
        r += 1
    K = [row[r:] for row in U]
    assert all(sum(q[a][i] * K[i][k] for i in range(N)) == 0 for a in range(h) for k in range(N - r))
    return K


# ---------------------------------------------------------------- the fan
class Fan:
    def __init__(self, cones, q):
        self.q = [[int(x) for x in r] for r in q]
        self.h, self.N = len(self.q), len(self.q[0])
        self.rays = gale_dual(self.q)                          # N x d
        self.d = len(self.rays[0])
        if self.d != 4: raise ValueError(f"expected a fourfold (rank of relations gives dimension {self.d})")
        self.cones = [tuple(sorted(int(x) for x in c)) for c in cones]
        self.coneset = {frozenset(c) for c in self.cones}
        self.faces = {k: {frozenset(s) for c in self.cones for s in itertools.combinations(c, k)} for k in (1, 2, 3)}
        for c in self.cones:
            if _det([self.rays[i] for i in c]) == 0: raise ValueError(f"cone {c} is not full-dimensional")

    def in_cone(self, s):
        return any(set(s) <= c for c in self.coneset)

    @lru_cache(maxsize=None)
    def I4(self, mono):
        """toric intersection number D_{i1} D_{i2} D_{i3} D_{i4} (mono: sorted tuple of 4 indices)"""
        dist = sorted(set(mono))
        if not self.in_cone(dist): return Fraction(0)
        if len(dist) == 4: return Fraction(1, abs(_det([self.rays[i] for i in dist])))
        i = next(x for x in dist if mono.count(x) >= 2)
        # m with <m, v_i> = 1, <m, v_j> = 0 for the other rays of the monomial; D_i = -sum_{j} <m, v_j> D_j (j outside)
        others = [j for j in dist if j != i]
        m = _particular([self.rays[i]] + [self.rays[j] for j in others], [1] + [0] * len(others), self.d)
        rest = list(mono); rest.remove(i)
        tot = Fraction(0)
        for j in range(self.N):
            if j in dist: continue
            mj = sum(a * b for a, b in zip(m, self.rays[j]))
            if mj == 0: continue
            tot -= mj * self.I4(tuple(sorted(rest + [j])))
        return tot

    def basis_divisors(self):
        """J_a as rational combinations of toric divisors: J_a = sum_i coef[a][i] D_i (via h11 independent columns)"""
        for B in itertools.combinations(range(self.N), self.h):
            M = [[self.q[a][i] for a in range(self.h)] for i in B]      # D_B = M J
            if _det(M) != 0:
                inv = [_solve(M, [int(a == b) for b in range(self.h)]) for a in range(self.h)]   # rows of M^{-1}... solve M x = e_a
                # J_a = sum_{k} (M^{-1})_{a,k} D_{B_k};  inv[a] solves M x = e_a -> x = column a of M^{-1}
                Minv = [[inv[c][r] for c in range(self.h)] for r in range(self.h)]   # Minv[r][c] = (M^{-1})_{rc}
                coef = [[Fraction(0)] * self.N for _ in range(self.h)]
                for a in range(self.h):
                    for k, i in enumerate(B): coef[a][i] = Minv[a][k]
                return coef
        raise ValueError("q has rank < h11")

    def integrate(self, *divs):
        """intersection of 4 divisors given as coefficient vectors over the toric divisors"""
        nz = [[(i, c) for i, c in enumerate(D) if c] for D in divs]
        tot = Fraction(0)
        for combo in itertools.product(*nz):
            coef = Fraction(1)
            for _, c in combo: coef *= c
            tot += coef * self.I4(tuple(sorted(i for i, _ in combo)))
        return tot

    def kappa(self):
        """CY triple intersections kappa_abc = J_a J_b J_c (sum_i D_i) in the basis of q"""
        J = self.basis_divisors(); X = [Fraction(1)] * self.N
        return {(a, b, c): self.integrate(J[a], J[b], J[c], X)
                for a in range(self.h) for b in range(a, self.h) for c in range(b, self.h)}

    def curve_coords(self, dvec):
        """curve class n (n_a = J_a . C) from its intersections with all toric divisors d_i = D_i . C"""
        return _solve([[self.q[a][i] for a in range(self.h)] for i in range(self.N)], dvec)

    def walls(self):
        return [s for s in self.faces[3] if sum(s < c for c in self.coneset) == 2]

    def wall_class(self, S):
        """exact D_i . V(S) for all i (Cox-Little-Schenck 6.4.4), as curve coordinates"""
        S = tuple(sorted(S)); Ss = set(S)
        apex = [next(iter(c - Ss)) for c in self.coneset if Ss < c]
        u1, u2 = apex
        lab = [u1, u2] + list(S)
        M = [[self.rays[l][k] for l in lab] for k in range(self.d)]          # d x 5
        lam = _particular([row[1:] for row in M], [-row[0] for row in M], 4)  # l_u1 = 1
        lam = [Fraction(1)] + lam
        mult_tau = 0
        for cols in itertools.combinations(range(self.d), 3):
            mult_tau = gcd(mult_tau, abs(_det([[self.rays[l][c] for c in cols] for l in S])))
        d_u1 = Fraction(mult_tau, abs(_det([self.rays[l] for l in [u1] + list(S)])))
        dvec = [Fraction(0)] * self.N
        for l, x in zip(lab, lam): dvec[l] = x * d_u1
        return self.curve_coords(dvec)

    def vex_faces(self):
        """faces (2 or 3 rays) of the fan whose rays share no facet of conv(rays)"""
        _, hyp = cgv_run._normaliz("cone", [list(v) + [1] for v in self.rays], self.d + 1)
        onf = [frozenset(k for k, h in enumerate(hyp) if sum(int(a) * b for a, b in zip(h, list(v) + [1])) == 0) for v in self.rays]
        return sorted((tuple(sorted(s)) for k in (2, 3) for s in self.faces[k] if not frozenset.intersection(*[onf[i] for i in s])),
                      key=lambda s: (len(s), s))

    def stratum_class(self, S, kappa):
        """|S| = 3: the exact class of the wall curve V_S (vex 2-cones are excluded, see prepare)"""
        assert len(S) == 3
        return self.wall_class(S)


# ---------------------------------------------------------------- preparing a run
def _primitive(v):
    den = 1
    for x in v: den = den * Fraction(x).denominator // gcd(den, Fraction(x).denominator)
    iv = [int(Fraction(x) * den) for x in v]; g = 0
    for x in iv: g = gcd(g, abs(x))
    return [x // g for x in iv]


def prepare(cones, q, kappa, generators=None, grading=None, saturate=True, check_kappa=True):
    fan = Fan(cones, q)
    kappa = {tuple(sorted(int(x) for x in k[:3])): int(k[3]) for k in kappa} if not isinstance(kappa, dict) else \
            {tuple(sorted(int(x) for x in k)): int(v) for k, v in kappa.items()}
    if check_kappa:
        kf = {k: v for k, v in fan.kappa().items() if v}
        bad = {k for k in set(kf) | set(kappa) if kf.get(k, 0) != kappa.get(k, 0)}
        if bad: raise ValueError(f"kappa does not match the fan's intersection ring (e.g. {sorted(bad)[:3]}); check that q and kappa use the same basis")
    walls = None
    if generators is None:
        if not saturate: raise ValueError("saturate=False needs explicit generators")
        walls = [_primitive(fan.wall_class(S)) for S in fan.walls()]
        gens = sorted({tuple(w) for w in walls if any(w)})
    else:
        gens = [tuple(int(x) for x in g) for g in generators]
    G = np.array(gens, dtype=np.int64)
    if saturate:
        G, hyp = cgv_run._normaliz("cone", G.tolist(), G.shape[1])
    else:
        _, hyp = cgv_run._normaliz("cone", G.tolist(), G.shape[1])
    if grading is None:
        grading = [int(x) for x in np.asarray(hyp).sum(0)]
    if not ((G @ np.array(grading)) > 0).all():
        raise ValueError("grading is not positive on every generator (pass grading=...)")
    vex = fan.vex_faces()
    if any(len(S) == 2 for S in vex):
        raise ValueError(f"vex 2-cone(s) {[S for S in vex if len(S) == 2]}: impossible for a fine fan on a reflexive polytope "
                         "(MacFadden-Sheridan, arXiv:2512.14817, Prop. 5); check that the fan is fine and uses every boundary "
                         "point not interior to a facet")
    strata = [(list(S), fan.stratum_class(S, kappa)) for S in vex]
    return dict(q=fan.q, generators=G.tolist(), mori_rays=[list(g) for g in gens], dual_rays=np.asarray(hyp).tolist(),
                grading_vector=list(map(int, grading)),
                intnums=[[a, b, c, v] for (a, b, c), v in kappa.items()], strata=strata,
                vex_cones=[S for S, _ in strata], saturate=saturate, fan=fan)


def cone_data_vex(d):
    base = cgv_run.cone_data(d)
    cones, tsets = list(base["cones"]), [list(T) + [-1] * (3 - len(T)) for T in base["tsets"]]
    G = np.array(d["generators"], dtype=np.int64); Q = np.array(d["q"], dtype=np.int64).T
    _, hyp = cgv_run._normaliz("cone", G.tolist(), G.shape[1])
    for T in d["vex_cones"]:
        if len(T) != 3: continue
        keep = [r for r in range(Q.shape[0]) if r not in T]
        B, _ = cgv_run._normaliz("inequalities", np.vstack([hyp, Q[keep]]).tolist(), G.shape[1])
        cones.append(B.tolist()); tsets.append(sorted(T))
    return dict(mori_hb=base["mori_hb"], cones=cones, tsets=tsets)


def write_input(d, max_deg, path, cones=True):
    cgv_run.write_input(d, max_deg, path, cones=False)
    with open(path, "a") as f:
        if cones:
            cd = cone_data_vex(d)
            f.write("2\n%d\n" % len(cd["mori_hb"]))
            f.writelines(" ".join(map(str, g)) + "\n" for g in cd["mori_hb"])
            f.write("%d\n" % -len(cd["cones"]))
            for B, T in zip(cd["cones"], cd["tsets"]):
                f.write(f"{len(B)} {T[0]} {T[1]} {T[2]}\n")
                f.writelines(" ".join(map(str, g)) + "\n" for g in B)
        else:
            f.write("0\n")
        f.write("%d\n" % len(d["strata"]))
        for S, vals in d["strata"]:
            f.write(f"{len(S)} {' '.join(map(str, S))} {' '.join(str(Fraction(v)) for v in vals)}\n")


def run(d, max_deg, threads=None, cone_data=True, binary=None, extra=()):
    """GVs {curve: value} for a prepared phase (binary: e.g. ../cgv_gpu with extra=["-g", "0"]).
    cone_data is ignored (full enumeration) when saturate=False."""
    use_cones = cone_data and d.get("saturate", True)
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        path = f.name
    write_input(d, max_deg, path, use_cones)
    p = subprocess.run([binary or BIN, "-t", str(threads or os.cpu_count()), *extra, path], capture_output=True, text=True)
    os.unlink(path)
    if p.returncode: raise RuntimeError(p.stderr[-2000:])
    out = {}
    for line in p.stdout.splitlines():
        v = line.split(); out[tuple(int(x) for x in v[:-1])] = int(v[-1])
    return out, p.stderr


def compute_gv(cones, q, kappa, max_deg, generators=None, grading=None, saturate=True, threads=None, check_kappa=True,
               binary=None, extra=()):
    d = prepare(cones, q, kappa, generators, grading, saturate, check_kappa)
    return run(d, max_deg, threads, binary=binary, extra=extra)[0]


# ---------------------------------------------------------------- CYTools helper (optional)
def ray_labels(p):
    """labels of p's boundary points not interior to facets (the rays of every phase), in label order"""
    nf = {tuple(int(x) for x in v) for v in p.points_not_interior_to_facets()}
    return [l for l in p.labels if l != p._label_origin and tuple(int(x) for x in p.points(which=l)[0]) in nf]


def phase_from_polytope(p, cones):
    """(cones as column indices, q, kappa) for the phase of CYTools polytope p (N lattice) with the given maximal
    cones (point labels, origin excluded). FRST or vex; CYTools computes q and kappa, the rest is derived above."""
    from cytools.triangulation import Triangulation
    org = p._label_origin; labels = ray_labels(p)
    cones = [tuple(sorted(int(x) for x in c)) for c in cones]
    tri = Triangulation(p, [org] + labels, simplices=[(org,) + c for c in cones], check_input_simplices=False)
    tv, cy = tri.get_toric_variety(), tri.get_cy()
    G = np.array(tv.glsm_charge_matrix()); lab2col = {l: i for i, l in enumerate(p.labels)}
    q = G[:, [lab2col[l] for l in labels]]
    col = {l: i for i, l in enumerate(labels)}
    kappa = {(int(a), int(b), int(c)): int(v) for (a, b, c), v in cy.intersection_numbers(in_basis=True).items()}
    return [tuple(col[l] for l in c) for c in cones], q.tolist(), kappa
