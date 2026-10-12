"""Run cgv on a CYTools CalabiYau or an input dict: the generators of the Mori cone (or of a face of it), q,
the intersection numbers and a grading vector (see write_input).

Library use:  compute_gvs(cy_or_dict, max_deg) -> {curve tuple: int GV}
              write_input(d, max_deg, path): cgv's input file, for running ./cgv from the shell
CLI:          cgv_run.py INPUTS.json MAX_DEG [--threads T] [--filter S] [--write FILE]
"""
import argparse, hashlib, itertools, json, os, shutil, subprocess, sys, tempfile, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "cgv")
CACHE = os.path.expanduser("~/.cache/cgv")
NORMALIZ = shutil.which("normaliz") or os.path.join(os.path.dirname(sys.executable), "normaliz")


def _read_matrices(path):
    """normaliz output file: blocks "m k rows [name]"; returns {name: int64 array (m x k)} ("rows" if unnamed)."""
    with open(path) as f:
        tok = f.read().split()
    out, i = {}, 0
    while i + 1 < len(tok):
        m, k = int(tok[i]), int(tok[i + 1])
        a = np.array(list(map(int, tok[i + 2:i + 2 + m * k])), dtype=np.int64).reshape(m, k)
        i += 2 + m * k
        name = "rows"
        if i < len(tok) and not tok[i].lstrip("-").isdigit():
            name, i = tok[i], i + 1
        out[name] = a
    return out


def _nz(n, data, hilbert_basis=False, extreme_rays=False):
    """Run normaliz in Z^n on data {input type: rows} (e.g. "cone", "inequalities", "equations"). Returns
    (Hilbert basis or None, extreme rays or None, support hyperplanes, equations); any description may be redundant."""
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, "c.in")
        body = "".join(f"{kind} {len(rows)}\n" + "".join(" ".join(map(str, r)) + "\n" for r in rows)
                       for kind, rows in data.items() if len(rows))
        with open(f, "w") as fh:
            fh.write(f"amb_space {n}\n{body}" + ("HilbertBasis\n" if hilbert_basis else "")
                     + ("ExtremeRays\n" if extreme_rays else "") + "SupportHyperplanes\n")
        subprocess.run([NORMALIZ, "-x=1", "--cst", *(["--gen"] if hilbert_basis else []),
                        *(["--ext"] if extreme_rays else []), f], check=True, capture_output=True)
        cst = _read_matrices(os.path.join(td, "c.cst"))
        hb = _read_matrices(os.path.join(td, "c.gen"))["rows"] if hilbert_basis else None
        ext = _read_matrices(os.path.join(td, "c.ext"))["rows"] if extreme_rays else None
    return hb, ext, cst["inequalities"], cst.get("equations", np.zeros((0, n), dtype=np.int64))


def _normaliz(kind, rows, n, hilbert_basis=True):
    """(Hilbert basis or None, support hyperplanes) of a full-dimensional cone given by spanning vectors ("cone")
    or inequalities ("inequalities"), either possibly redundant."""
    hb, _, hyp, _ = _nz(n, {kind: rows}, hilbert_basis=hilbert_basis)
    return hb, hyp


def _cone(d):
    """The cone the curves live on (the Mori cone or a face of it), from d["mori_rays"] (vectors spanning it,
    redundant allowed) or d["mori_hyperplanes"] (inequalities, redundant allowed): (extreme rays, inequalities,
    equations), the equations nonempty for a cone of lower dimension, e.g. a face."""
    n = len(d["q"])
    if d.get("mori_rays") is not None:
        G = np.array(d["mori_rays"], dtype=np.int64).reshape(-1, n)
        G = G[(G != 0).any(1)]
        if not len(G): raise ValueError("mori_rays: no nonzero vector")
        _, ext, hyp, eq = _nz(n, {"cone": G.tolist()}, extreme_rays=True)
    elif d.get("mori_hyperplanes") is not None:
        _, ext, hyp, eq = _nz(n, {"inequalities": np.asarray(d["mori_hyperplanes"]).astype(int).tolist()}, extreme_rays=True)
    else:
        raise ValueError("no description of the Mori cone (or face): give mori_rays or mori_hyperplanes")
    if not len(ext):
        raise ValueError("the cone must be pointed and nonzero (a Mori cone or a face of it)")
    return ext, hyp, eq


def mori_facets(d):
    """Inequalities H x >= 0 describing exactly the cone the curves live on (its facets, and both signs of its
    equations if it is a face of lower dimension)."""
    _, hyp, eq = _cone(d)
    return np.vstack([hyp, eq, -eq])


def _lattice_basis(rows):
    """Echelon basis (python ints) of the lattice the integer rows generate, and its pivot columns."""
    A = [[int(x) for x in r] for r in rows]
    A = [r for r in A if any(r)]
    basis, piv, col, n = [], [], 0, len(A[0]) if A else 0
    while A and col < n:
        while sum(1 for r in A if r[col]) > 1:   # Euclid on column col
            b = min((r for r in A if r[col]), key=lambda r: abs(r[col]))
            A = [r if r is b or not r[col] else [x - (r[col] // b[col]) * y for x, y in zip(r, b)] for r in A]
            A = [r for r in A if any(r)]
        b = next((r for r in A if r[col]), None)
        if b is not None:
            basis.append(b if b[col] > 0 else [-x for x in b]); piv.append(col)
            A = [r for r in A if r is not b]
        col += 1
    return basis, piv


def _coords(basis, piv, x):
    """Coordinates of the lattice vector x in the echelon basis."""
    x, c = [int(v) for v in x], []
    for b, j in zip(basis, piv):
        q, r = divmod(x[j], b[j])
        if r: raise ValueError("vector not in the lattice")
        c.append(q); x = [u - q * v for u, v in zip(x, b)]
    if any(x): raise ValueError("vector not in the lattice")
    return c


def cone_data(d):
    """Hilbert bases of the cones K_T = cone cap {Q_r . C >= 0 for r not in T} for |T| = 1, 2 and the vex 3-cones
    d["vex_cones"] if any (cone: the Mori cone or a face of it, see _cone), each cone once; the cone's extreme rays;
    and inequalities describing the cone (for the lightcone section). Computed by normaliz; cached per cone and Q.

    Only divisors negative somewhere on the cone (on one of its extreme rays) are put in T, and only Q_r of those
    enter the inequalities: Q_r >= 0 holds on the whole cone for the others. The cone K_T is K_M, M the divisors
    negative on its Hilbert basis (every point of K_T is a nonnegative combination of it, so its negative set lies
    in M; and K_M is inside K_T as M is inside T). So each K_T is written once, labeled by one T with K_T = K_M,
    preferring a label that is not a pair: cgv skips a pair cone whose quadratic form vanishes (see gv.c,
    build_candidates), and the union of the lattice points of the other cones is all it uses. The cones with
    M empty lie inside every K_{r}; when no divisor is ever negative the cone itself is the only one."""
    Q = np.array(d["q"], dtype=np.int64).T  # divisor rows
    ndiv, n = Q.shape
    vex = sorted(set(tuple(sorted(int(r) for r in T)) for T in d.get("vex_cones") or [] if len(T) == 3))
    desc = next(((k, d[k]) for k in ("mori_rays", "mori_hyperplanes") if d.get(k) is not None), None)
    if desc is None:
        raise ValueError("no description of the Mori cone (or face): give mori_rays or mori_hyperplanes")
    h = hashlib.sha1(json.dumps(["kt2", desc[0], np.asarray(desc[1]).astype(int).tolist(), Q.tolist(), vex]).encode()).hexdigest()[:16]
    path = os.path.join(CACHE, f"kt2_{h}.json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    t0 = time.time()
    ext, hyp, eq = _cone(d)
    facets, rays = np.vstack([hyp, eq, -eq]), ext.tolist()
    neg = [r for r in range(ndiv) if (ext @ Q[r] < 0).any()]
    nset = set(neg)
    calls = 1
    L = None
    if len(eq):   # a face (lower dimension): work in a basis L of the lattice it spans (much faster normaliz runs)
        hb_cone = _nz(n, {"cone": ext.tolist()}, hilbert_basis=True)[0]
        calls += 1
        L, piv = _lattice_basis(hb_cone.tolist())  # the Hilbert basis generates span cap Z^n
        Q = np.array([[sum(int(a) * int(b) for a, b in zip(q, l)) for l in L] for q in Q.tolist()], dtype=np.int64)
        ext = np.array([_coords(L, piv, x) for x in ext.tolist()], dtype=np.int64)
        n = len(L)
        _, _, hyp, eq = _nz(n, {"cone": ext.tolist()})
        calls += 1
    if not neg:   # nothing is ever negative: the lattice points of the cone itself
        B = _nz(n, {"inequalities": hyp.tolist(), "equations": eq.tolist()}, hilbert_basis=True)[0]
        groups = {(): ((0,), B)}
        calls += 1
    else:
        hb = {}       # E = T cap neg -> (Hilbert basis of K_E, its negative set M)
        groups = {}   # M -> (label T, Hilbert basis)
        for T in [(r,) for r in neg] + list(itertools.combinations(neg, 2)) + vex:
            E = tuple(r for r in T if r in nset)
            if not E:
                continue
            if E not in hb:
                ineq = np.vstack([hyp, Q[[r for r in neg if r not in E]]])
                B = _nz(n, {"inequalities": ineq.tolist(), "equations": eq.tolist()}, hilbert_basis=True)[0]
                calls += 1
                hb[E] = (B, tuple(r for r in E if (B @ Q[r] < 0).any()))
            B, M = hb[E]
            if not M:
                continue
            if M not in groups or (len(groups[M][0]) == 2 and len(T) != 2):
                groups[M] = (T, B)
    lift = (lambda B: B.tolist()) if L is None else (lambda B: (B.astype(object) @ np.array(L, dtype=object)).tolist())
    out = dict(cones=[lift(B) for _, B in groups.values()], tsets=[list(T) for T, _ in groups.values()],
               rays=rays, facets=facets.tolist())
    sys.stderr.write(f"cgv_run: {len(out['cones'])} cones K_T ({len(neg)} of {ndiv} divisors negative somewhere on the "
                     f"cone), {calls} normaliz runs, {time.time() - t0:.1f}s\n")
    os.makedirs(CACHE, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def default_grading(rays):
    """A small integer grading vector positive on the cone spanned by rays (CYTools' Cone.find_grading_vector)."""
    try:
        from cytools import Cone
    except ImportError:
        raise ValueError("no grading_vector given, and CYTools (which would choose one) is not installed") from None
    return [int(x) for x in Cone(rays=rays).find_grading_vector()]


def write_input(d, max_deg, path):
    """Write cgv's input file (format: the comment at the top of gv.c) for the input dict d:

      mori_rays        vectors spanning the cone the curves live on: the Mori cone, or a face of it (curves on a
                       face are closed under decomposition, so their GVs are exact); redundant vectors allowed.
                       Or mori_hyperplanes: inequalities H x >= 0 describing it.
      q                GLSM charges, h11 rows x ndiv columns (CYTools' curve_basis(as_matrix=True) layout)
      intnums          triple intersection numbers in the same basis: [[i, j, k, value], ...] (each unordered
                       triple once) or {(i, j, k): value}
      grading_vector   integer vector positive on the cone; default: default_grading(the cone's extreme rays)
      lightcone        optional: points p; only the curves C with p - C in the cone (whose GVs are then exact)
      strata, vex_cones   optional, vex phases (written by cgv_phase.py)

    max_deg: largest degree; optional with lightcone (default: the largest degree of its points)."""
    from fractions import Fraction
    q = d["q"]                       # h11 rows x ndiv columns
    h11, ndiv = len(q), len(q[0])
    intnums = d["intnums"]
    if isinstance(intnums, dict):
        intnums = [[*k, v] for k, v in intnums.items()]
    cd = cone_data(d)
    w = d.get("grading_vector")
    if w is None:
        w = default_grading(cd["rays"])
    w = [int(x) for x in w]
    if len(w) != h11:
        raise ValueError(f"grading_vector has {len(w)} entries, q has {h11} rows")
    if (np.array(cd["rays"], dtype=np.int64) @ np.array(w, dtype=np.int64) <= 0).any():
        raise ValueError("grading_vector is not positive on the cone (on some generator of it)")
    if max_deg is None:              # lightcone GVs: up to the largest degree of the chosen points
        if d.get("lightcone") is None: raise ValueError("max_deg is needed unless lightcone points are given")
        max_deg = max(int(np.dot(p, w)) for p in d["lightcone"])
    lines = [f"{h11} {ndiv} {max_deg}", " ".join(map(str, w))]
    for r in range(ndiv):            # cgv wants divisor rows
        lines.append(" ".join(str(q[a][r]) for a in range(h11)))
    lines.append(str(len(intnums)))
    lines += [" ".join(str(int(v)) for v in x) for x in intnums]
    lines.append(str(len(cd["cones"])))
    for B, T in zip(cd["cones"], cd["tsets"]):
        lines.append(f"{len(B)} {len(T)} {' '.join(map(str, T))}")
        lines += [" ".join(map(str, g)) for g in B]
    if d.get("lightcone") is not None:
        lines.append(f"lightcone {len(d['lightcone'])}")
        lines += [" ".join(str(int(x)) for x in p) for p in d["lightcone"]]
        lines.append(str(len(cd["facets"])))
        lines += [" ".join(map(str, h)) for h in cd["facets"]]
    if d.get("strata"):
        lines.append(f"vex {len(d['strata'])}")
        lines += [f"{len(S)} {' '.join(map(str, S))} {' '.join(str(Fraction(v)) for v in vals)}" for S, vals in d["strata"]]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def run_cgv(d, max_deg, threads=1, extra=(), env=None, binary=None):
    """Write d's input (see write_input), run cgv on it (binary: default BIN; extra: more cgv arguments; env: extra
    environment variables). Returns ({curve tuple: GV}, cgv's wall time, cgv's stderr)."""
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        path = f.name
    try:
        write_input(d, max_deg, path)
        t0 = time.time()
        p = subprocess.run([binary or BIN, "-t", str(threads), *extra, path], capture_output=True, text=True,
                           env=dict(os.environ, **env) if env else None)
        dt = time.time() - t0
    finally:
        os.unlink(path)
    if p.returncode:
        raise RuntimeError(p.stderr[-2000:])
    for line in p.stderr.splitlines():
        if line.startswith("cgv: note:"):
            import warnings
            warnings.warn(line[len("cgv: note: "):], stacklevel=2)
    out = {}
    for line in p.stdout.splitlines():
        v = line.split()
        out[tuple(int(x) for x in v[:-1])] = int(v[-1])
    return out, dt, p.stderr


def cy_input(cy, grading_vec=None):
    """cgv input dict for a cytools CalabiYau (threefold hypersurface): the Mori cone (its rays; cgv needs only a
    description of the cone, not its lattice points), grading vector, GLSM charges and intersection numbers."""
    mori = cy.mori_cone_cap(in_basis=True)
    if grading_vec is None:
        grading_vec = mori.find_grading_vector()
    return dict(
        mori_rays=np.asarray(mori.rays()).astype(int).tolist(),
        grading_vector=[int(x) for x in grading_vec],
        q=np.asarray(cy.curve_basis(include_origin=False, as_matrix=True)).astype(int).tolist(),
        intnums=[[int(i), int(j), int(k), int(x)] for (i, j, k), x in
                 cy.intersection_numbers(in_basis=True, format="dok").items()],
    )


def pick_device(device, gpu_bin, hip=False):
    """cgv's device arguments: ["-g", "N"] for GPU N, [] for the CPU.

    device: "cpu", "gpu" (GPU 0), "gpu:N", or "auto": the NVIDIA GPU with the most free memory (more than
    4 GB) that does not drive a display, if the GPU build gpu_bin exists, else the CPU. With hip=True (gpu_bin
    is an AMD build, which nvidia-smi cannot see) "auto" means GPU 0; cgv itself falls back to the CPU when
    the GPU is busy, too small or drives a display.
    """
    if device == "auto":
        device = "cpu"
        if gpu_bin and os.path.exists(gpu_bin) and hip:
            device = "gpu:0"
        elif gpu_bin and os.path.exists(gpu_bin) and shutil.which("nvidia-smi"):
            # the GPU with the most free memory (cgv_gpu itself falls back to the CPU if it is too busy)
            try:
                # never pick a GPU that drives a display: filling its memory can black out the desktop
                q = subprocess.run(["nvidia-smi", "--query-gpu=index,memory.free,display_active", "--format=csv,noheader,nounits"],
                                   capture_output=True, text=True, check=True).stdout
                free = [(int(r[0]), int(r[1])) for r in (l.split(", ") for l in q.strip().splitlines())
                        if r[2].strip() != "Enabled"]
                if free:
                    best = max(free, key=lambda r: r[1])
                    if best[1] > 4000:
                        device = f"gpu:{best[0]}"
            except (subprocess.CalledProcessError, ValueError):
                pass
    if device.startswith("gpu"):
        return ["-g", device.split(":")[1] if ":" in device else "0"]
    return []


def compute_gvs(cy_or_input, max_deg=None, grading_vec=None, device="auto", threads=None, lanes=None, verbose=False,
                low_memory=False, lightcone=None):
    """GV invariants with cgv. Returns {curve tuple: int GV} (nonzero only), like
    cytools' cy.compute_gvs(...).dok, i.e. the same as cygv.

    cy_or_input: a cytools CalabiYau, or an input dict (see write_input): dict(mori_rays=G, q=Q, intnums=kappa,
            grading_vector=w), G the generators of the Mori cone or of a face of it.
    max_deg: largest degree computed; optional with lightcone (default: the largest degree of its points).
    lightcone: curves p; computes only their backward lightcones, the curves C with p - C in the Mori cone
            (exactly the curves whose GVs enter those of p), so the GVs are exact.
    device: "cpu", "gpu" (GPU 0), "gpu:N", or "auto" (GPU if a CUDA build and GPU are
            present; cgv itself keeps small or sparse-degree problems on the CPU).
    threads: CPU threads (default: all cores). Fewer threads use proportionally less host memory
            (each thread keeps its own working tables) and take longer.
    lanes:  number of ~62-bit primes per pass (default: chosen by a cheap probe).
    low_memory: low-memory mode (CGV_MEM=low): no curve-class replay (CPU and GPU), and on Linux/glibc freed
            memory goes back to the system at once. Lower peak memory for somewhat more time.
    """
    d = cy_or_input if isinstance(cy_or_input, dict) else cy_input(cy_or_input, grading_vec)
    if grading_vec is not None and isinstance(cy_or_input, dict):
        d = dict(d, grading_vector=[int(x) for x in grading_vec])
    if lightcone is not None:
        d = dict(d, lightcone=[[int(x) for x in p] for p in lightcone])
    gpu_bin = os.path.join(HERE, "..", "cgv_gpu")
    extra = pick_device(device, gpu_bin)
    binary = gpu_bin if extra else os.path.join(HERE, "..", "cgv")
    if lanes:
        extra += ["-l", str(lanes)]
    out, dt, err = run_cgv(d, max_deg, threads or os.cpu_count(), extra=extra, binary=binary,
                           env={"CGV_MEM": "low"} if low_memory else None)
    if verbose:
        sys.stderr.write(err)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="cgv on input dicts (see write_input) in a JSON file: one dict or a list")
    ap.add_argument("inputs"); ap.add_argument("max_deg", type=int)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--filter", default="", help="only the inputs whose name contains this")
    ap.add_argument("--write", metavar="FILE", help="write cgv's input file (one input) instead of running cgv")
    ap.add_argument("--stderr", action="store_true", help="show cgv's log")
    a = ap.parse_args()
    with open(a.inputs) as f:
        ds = json.load(f)
    ds = [ds] if isinstance(ds, dict) else ds
    if a.write:
        if len(ds) != 1: sys.exit("--write takes a file with one input")
        write_input(ds[0], a.max_deg, a.write)
        sys.exit()
    for i, d in enumerate(ds):
        tag = str(d.get("name", i)) + (f":{d['grading']}" if "grading" in d else "")
        if a.filter not in tag:
            continue
        mine, dt, err = run_cgv(d, a.max_deg, a.threads)
        if a.stderr:
            sys.stderr.write(err)
        print(f"{tag:28s} D={a.max_deg:<5d} cgv {dt:8.3f}s n={len(mine)}", flush=True)
