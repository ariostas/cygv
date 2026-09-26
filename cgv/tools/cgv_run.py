"""Run cgv on a cygv-style input dict, and optionally compare against cygv.

Library use:  run_cgv(d, max_deg, threads) -> {curve tuple: int GV}
CLI:          cgv_run.py INPUTS.json MAX_DEG [--threads T] [--compare] [--filter S]
"""
import argparse, hashlib, itertools, json, os, shutil, subprocess, sys, tempfile, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "cgv")
CACHE = os.path.expanduser("~/.cache/cgv")
NORMALIZ = shutil.which("normaliz") or os.path.join(os.path.dirname(sys.executable), "normaliz")


def _normaliz(kind, rows, n):
    """Run normaliz on a cone given by generators ("cone") or inequalities
    ("inequalities"); return (hilbert_basis, support_hyperplanes)."""
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, "c.in")
        body = "\n".join(" ".join(map(str, r)) for r in rows)
        open(f, "w").write(f"amb_space {n}\n{kind} {len(rows)}\n{body}\nHilbertBasis\nSupportHyperplanes\n")
        subprocess.run([NORMALIZ, "-x=1", "-f", f], check=True, capture_output=True)
        g = open(os.path.join(td, "c.gen")).read().split()
        c = open(os.path.join(td, "c.cst")).read().split()
    m, k = int(g[0]), int(g[1])
    hb = np.array(list(map(int, g[2:2 + m * k])), dtype=np.int64).reshape(m, k)
    m2, k2 = int(c[0]), int(c[1])
    hyp = np.array(list(map(int, c[2:2 + m2 * k2])), dtype=np.int64).reshape(m2, k2)
    return hb, hyp


def cone_data(d):
    """Hilbert bases of the Mori cone and of the cones
    K_T = Mori cone cap {Q_r . C >= 0 for r not in T}, |T| = 2. Cached per geometry."""
    G = np.array(d["generators"], dtype=np.int64)
    G = G[(G != 0).any(1)]
    Q = np.array(d["q"], dtype=np.int64).T  # divisor rows
    h = hashlib.sha1(G.tobytes() + Q.tobytes()).hexdigest()[:16]
    path = os.path.join(CACHE, f"cones2_{h}.json")
    if os.path.exists(path):
        return json.load(open(path))
    n = G.shape[1]
    hb, hyp = _normaliz("cone", G.tolist(), n)
    cones, tsets = [], []
    for size in (1, 2):
        for T in itertools.combinations(range(Q.shape[0]), size):
            keep = [r for r in range(Q.shape[0]) if r not in T]
            B, _ = _normaliz("inequalities", np.vstack([hyp, Q[keep]]).tolist(), n)
            cones.append(B.tolist())
            tsets.append(list(T))
    out = dict(mori_hb=hb.tolist(), cones=cones, tsets=tsets)
    os.makedirs(CACHE, exist_ok=True)
    json.dump(out, open(path, "w"))
    return out


def write_input(d, max_deg, path, cones=True):
    q = d["q"]                       # cytools layout: h11 rows x ndiv columns
    h11, ndiv = len(q), len(q[0])
    lines = [f"{h11} {ndiv} {max_deg}", " ".join(map(str, d["grading_vector"]))]
    for r in range(ndiv):            # cgv wants divisor rows
        lines.append(" ".join(str(q[a][r]) for a in range(h11)))
    lines.append(str(len(d["generators"])))
    lines += [" ".join(map(str, g)) for g in d["generators"]]
    lines.append(str(len(d["intnums"])))
    lines += [" ".join(map(str, x)) for x in d["intnums"]]
    if cones:
        cd = cone_data(d)
        lines.append("1")
        lines.append(str(len(cd["mori_hb"])))
        lines += [" ".join(map(str, g)) for g in cd["mori_hb"]]
        lines.append(str(-len(cd["cones"])))  # negative count: cones carry their divisor sets
        for B, T in zip(cd["cones"], cd["tsets"]):
            T2 = list(T) + [-1] * (2 - len(T))
            lines.append(f"{len(B)} {T2[0]} {T2[1]}")
            lines += [" ".join(map(str, g)) for g in B]
    open(path, "w").write("\n".join(lines) + "\n")


def run_cgv(d, max_deg, threads=1, extra=(), cones=True, env=None):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        path = f.name
    write_input(d, max_deg, path, cones)
    t0 = time.time()
    p = subprocess.run([BIN, "-t", str(threads), *extra, path], capture_output=True, text=True,
                       env=dict(os.environ, **env) if env else None)
    if p.returncode:
        raise RuntimeError(p.stderr[-2000:])
    dt = time.time() - t0
    os.unlink(path)
    for line in p.stderr.splitlines():
        if line.startswith("cgv: note:"):
            import warnings
            warnings.warn(line[len("cgv: note: "):], stacklevel=2)
    out = {}
    for line in p.stdout.splitlines():
        v = line.split()
        out[tuple(int(x) for x in v[:-1])] = int(v[-1])
    return out, dt, p.stderr


def cy_input(cy, grading_vec=None, min_points=None):
    """cgv/cygv input dict for a cytools CalabiYau (threefold hypersurface), built the
    same way cytools' own compute_gvs builds its cygv call."""
    mori = cy.mori_cone_cap(in_basis=True)
    pts = mori.find_lattice_points(min_points=min_points or 100 * cy.h11())
    gens = np.vstack([np.asarray(mori.rays()), np.asarray(pts)])
    if grading_vec is None:
        grading_vec = mori.find_grading_vector()
    return dict(
        generators=gens.astype(int).tolist(),
        grading_vector=[int(x) for x in grading_vec],
        q=np.asarray(cy.curve_basis(include_origin=False, as_matrix=True)).astype(int).tolist(),
        intnums=[[int(i), int(j), int(k), int(x)] for (i, j, k), x in
                 cy.intersection_numbers(in_basis=True, format="dok").items()],
    )


def compute_gvs(cy_or_input, max_deg, grading_vec=None, device="auto", threads=None, lanes=None, verbose=False,
                low_memory=False):
    """GV invariants with cgv. Returns {curve tuple: int GV} (nonzero only), like
    cytools' cy.compute_gvs(...).dok, i.e. the same as cygv.

    cy_or_input: a cytools CalabiYau, or an input dict (see cy_input).
    device: "cpu", "gpu" (GPU 0), "gpu:N", or "auto" (GPU if a CUDA build and GPU are
            present; cgv itself keeps small or sparse-degree problems on the CPU).
    lanes:  number of ~62-bit primes per pass (default: chosen by a cheap probe).
    low_memory: return freed memory to the system at once (Linux/glibc; CGV_LOW_MEM=1): lower peak host
            memory for somewhat more time. No effect on macOS.
    """
    global BIN
    d = cy_or_input if isinstance(cy_or_input, dict) else cy_input(cy_or_input, grading_vec)
    if grading_vec is not None and isinstance(cy_or_input, dict):
        d = dict(d, grading_vector=[int(x) for x in grading_vec])
    extra = []
    gpu_bin = os.path.join(HERE, "..", "cgv_gpu")
    if device == "auto":
        device = "cpu"
        if os.path.exists(gpu_bin) and shutil.which("nvidia-smi"):
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
        dev = device.split(":")[1] if ":" in device else "0"
        extra = ["-g", dev]
        binary = gpu_bin
    else:
        binary = os.path.join(HERE, "..", "cgv")
    if lanes:
        extra += ["-l", str(lanes)]
    old = BIN
    BIN = binary
    try:
        out, dt, err = run_cgv(d, max_deg, threads or os.cpu_count(), extra=extra,
                               env={"CGV_LOW_MEM": "1"} if low_memory else None)
    finally:
        BIN = old
    if verbose:
        sys.stderr.write(err)
    return out


def run_cygv(d, max_deg):
    import cygv
    t0 = time.time()
    res = cygv.compute_gv(generators=d["generators"], grading_vector=d["grading_vector"], q=d["q"],
                          intnums={(i, j, k): x for i, j, k, x in d["intnums"]}, max_deg=max_deg)
    return {tuple(v): int(g) for v, g in res}, time.time() - t0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs"); ap.add_argument("max_deg", type=int)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--filter", default="")
    ap.add_argument("--stderr", action="store_true")
    a = ap.parse_args()
    for d in json.load(open(a.inputs)):
        tag = f"{d['name']}:{d['grading']}"
        if a.filter not in tag:
            continue
        mine, dt, err = run_cgv(d, a.max_deg, a.threads)
        if a.stderr:
            sys.stderr.write(err)
        msg = f"{tag:28s} D={a.max_deg:<5d} cgv {dt:8.3f}s n={len(mine)}"
        if a.compare:
            ref, rt = run_cygv(d, a.max_deg)
            ok = ref == mine
            msg += f" | cygv {rt:8.3f}s n={len(ref)} | match={ok} speedup={rt/dt:.1f}x"
            if not ok:
                bad = [k for k in set(ref) | set(mine) if ref.get(k) != mine.get(k)]
                msg += f" ndiff={len(bad)} e.g. {[(k, ref.get(k), mine.get(k)) for k in bad[:3]]}"
        print(msg, flush=True)
