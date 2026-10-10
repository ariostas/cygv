"""Run cgv on an input dict (Mori cone description, grading, Q, intersection numbers).

Library use:  run_cgv(d, max_deg, threads) -> {curve tuple: int GV}
CLI:          cgv_run.py INPUTS.json MAX_DEG [--threads T] [--filter S]
"""
import argparse, hashlib, itertools, json, os, shutil, subprocess, sys, tempfile, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "cgv")
CACHE = os.path.expanduser("~/.cache/cgv")
NORMALIZ = shutil.which("normaliz") or os.path.join(os.path.dirname(sys.executable), "normaliz")


def _normaliz(kind, rows, n, hilbert_basis=True):
    """Run normaliz on a cone given by vectors spanning it ("cone") or by inequalities ("inequalities"); return
    (hilbert_basis or None, support_hyperplanes). Either description may be redundant."""
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, "c.in")
        body = "\n".join(" ".join(map(str, r)) for r in rows)
        open(f, "w").write(f"amb_space {n}\n{kind} {len(rows)}\n{body}\n" + ("HilbertBasis\n" if hilbert_basis else "")
                           + "SupportHyperplanes\n")
        subprocess.run([NORMALIZ, "-x=1", "-f", f], check=True, capture_output=True)
        hb = None
        if hilbert_basis:
            g = open(os.path.join(td, "c.gen")).read().split()
            m, k = int(g[0]), int(g[1])
            hb = np.array(list(map(int, g[2:2 + m * k])), dtype=np.int64).reshape(m, k)
        c = open(os.path.join(td, "c.cst")).read().split()
    m2, k2 = int(c[0]), int(c[1])
    hyp = np.array(list(map(int, c[2:2 + m2 * k2])), dtype=np.int64).reshape(m2, k2)
    return hb, hyp


def mori_facets(d):
    """Facet normals of the cone the GVs live on (inward: H x >= 0), from whichever description d has:
    "mori_rays" (any vectors spanning the Mori cone, redundant allowed), "mori_hyperplanes" (any inequalities,
    redundant allowed)."""
    n = len(d["q"])
    if d.get("mori_rays") is not None:
        G = np.array(d["mori_rays"], dtype=np.int64).reshape(-1, n)
        G = G[(G != 0).any(1)]
        return _normaliz("cone", G.tolist(), n, hilbert_basis=False)[1]
    if d.get("mori_hyperplanes") is not None:
        return _normaliz("inequalities", d["mori_hyperplanes"], n, hilbert_basis=False)[1]
    raise ValueError("no description of the Mori cone: give mori_rays or mori_hyperplanes")


def cone_data(d):
    """Hilbert bases of the cones K_T = Mori cone cap {Q_r . C >= 0 for r not in T}: |T| = 1, 2, and the vex
    3-cones in d["vex_cones"] if any, and the Mori cone's facets. Computed by normaliz; cached per cone and Q."""
    Q = np.array(d["q"], dtype=np.int64).T  # divisor rows
    vex = sorted(tuple(sorted(T)) for T in d.get("vex_cones") or [] if len(T) == 3)
    desc = next((k, d[k]) for k in ("mori_rays", "mori_hyperplanes") if d.get(k) is not None)
    h = hashlib.sha1(json.dumps([desc[0], np.asarray(desc[1]).astype(int).tolist(), Q.tolist(), vex]).encode()).hexdigest()[:16]
    path = os.path.join(CACHE, f"kt_{h}.json")
    if os.path.exists(path):
        out = json.load(open(path))
        if "facets" in out:
            return out
        hyp = mori_facets(d)   # a cache file from before the facets were stored
        out["facets"] = hyp.tolist()
        json.dump(out, open(path, "w"))
        return out
    hyp = mori_facets(d)
    n = hyp.shape[1]
    cones, tsets = [], []
    for T in [T for size in (1, 2) for T in itertools.combinations(range(Q.shape[0]), size)] + vex:
        keep = [r for r in range(Q.shape[0]) if r not in T]
        B, _ = _normaliz("inequalities", np.vstack([hyp, Q[keep]]).tolist(), n)
        cones.append(B.tolist())
        tsets.append(list(T))
    out = dict(cones=cones, tsets=tsets, facets=hyp.tolist())
    os.makedirs(CACHE, exist_ok=True)
    json.dump(out, open(path, "w"))
    return out


def write_input(d, max_deg, path):
    """cgv's input ("cgv 2", see gv.c): grading, Q, intersection numbers, the K_T Hilbert bases, and the optional
    lightcone points (d["lightcone"]: keep only the curves C with p - C in the Mori cone for some p) and vex strata."""
    from fractions import Fraction
    q = d["q"]                       # cytools layout: h11 rows x ndiv columns
    h11, ndiv = len(q), len(q[0])
    if max_deg is None:              # lightcone GVs: up to the largest degree of the chosen points
        if d.get("lightcone") is None: raise ValueError("max_deg is needed unless lightcone points are given")
        max_deg = max(int(np.dot(p, d["grading_vector"])) for p in d["lightcone"])
    lines = ["cgv 2", f"{h11} {ndiv} {max_deg}", " ".join(map(str, d["grading_vector"]))]
    for r in range(ndiv):            # cgv wants divisor rows
        lines.append(" ".join(str(q[a][r]) for a in range(h11)))
    lines.append(str(len(d["intnums"])))
    lines += [" ".join(map(str, x)) for x in d["intnums"]]
    cd = cone_data(d)
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
    open(path, "w").write("\n".join(lines) + "\n")


def run_cgv(d, max_deg, threads=1, extra=(), env=None):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        path = f.name
    write_input(d, max_deg, path)
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

    cy_or_input: a cytools CalabiYau, or an input dict (see cy_input).
    max_deg: largest degree computed; optional with lightcone (default: the largest degree of its points).
    lightcone: curves p; computes only their backward lightcones, the curves C with p - C in the Mori cone
            (exactly the curves whose GVs enter those of p), so the GVs are exact.
    device: "cpu", "gpu" (GPU 0), "gpu:N", or "auto" (GPU if a CUDA build and GPU are
            present; cgv itself keeps small or sparse-degree problems on the CPU).
    threads: CPU threads (default: all cores). Fewer threads use proportionally less host memory
            (each thread keeps its own working tables) and take longer.
    lanes:  number of ~62-bit primes per pass (default: chosen by a cheap probe).
    low_memory: return freed memory to the system at once (Linux/glibc; CGV_LOW_MEM=1): lower peak host
            memory for somewhat more time. No effect on macOS.
    """
    global BIN
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


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs"); ap.add_argument("max_deg", type=int)
    ap.add_argument("--threads", type=int, default=1)
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
        print(f"{tag:28s} D={a.max_deg:<5d} cgv {dt:8.3f}s n={len(mine)}", flush=True)
