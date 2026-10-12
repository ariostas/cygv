"""Regression for phases (FRST and vex) through tools/cgv_phase.py (needs numpy and normaliz, not CYTools).

phase refs: recompute each fan's GVs from (cones, q, kappa, grading) and compare with the stored table.
            Also: the backward lightcones of three stored curves give the stored GVs on those lightcones.
pair refs:  a vex phase and an FRST of a different polytope with the same CY: compute both under one grading
            (w_F = T^{-1} w_V) and compare class by class (n_F = T^T n_V); also compare the vex side with the stored table.
Usage: vex_regress.py [--threads T] [--filter S] [--bin PATH] [--extra "-g 0"]. Exits nonzero on any mismatch."""
import argparse, glob, gzip, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import cgv_phase

ap = argparse.ArgumentParser()
ap.add_argument("--threads", type=int, default=os.cpu_count())
ap.add_argument("--filter", default="")
ap.add_argument("--bin", default=None)
ap.add_argument("--extra", default="")
a = ap.parse_args()
kw = dict(threads=a.threads, binary=a.bin, extra=a.extra.split())
fails = 0
for path in sorted(glob.glob(os.path.join(os.path.dirname(__file__), "refs_vex", "*.json.gz"))):
    r = json.load(gzip.open(path, "rt"))
    if a.filter not in r["name"]: continue
    t = time.time()
    if r["kind"] == "phase":
        got = cgv_phase.compute_gv(r["cones"], r["q"], r["kappa"], r["max_deg"], grading=r["grading"], **kw)
        want = {tuple(k): v for k, v in r["gvs"]}
        ok = got == want; msg = f"{len(want)} GVs, vex cones {r['vex_cones']}"
        # backward lightcones of three stored curves (top, middle and low degree): equal to the stored GVs there
        W = np.array(r["grading"]); byd = sorted(want, key=lambda c: (W @ c, c))
        pts = [byd[-1], byd[len(byd) // 2], byd[len(byd) // 4]]
        d = cgv_phase.prepare(r["cones"], r["q"], r["kappa"], grading=r["grading"], lightcone=pts)
        H = np.array(cgv_phase.cgv_run.cone_data(d)["facets"])
        inl = lambda c: any((H @ (np.array(p) - c) >= 0).all() for p in pts)
        lc = cgv_phase.run(d, None, kw["threads"], binary=kw["binary"], extra=kw["extra"])[0]
        lok = lc == {c: v for c, v in want.items() if inl(np.array(c))}
        ok = ok and lok; msg += f"; lightcones of 3 curves: {len(lc)} GVs {'same' if lok else 'DIFFER'}"
    else:
        T = np.array(r["T"]); wV = np.array(r["grading_vex"]); wF = np.round(np.linalg.inv(T) @ wV).astype(int)
        gV = cgv_phase.compute_gv(r["vex"]["cones"], r["vex"]["q"], r["vex"]["kappa"], r["max_deg"], grading=wV.tolist(), **kw)
        gF = cgv_phase.compute_gv(r["frst"]["cones"], r["frst"]["q"], r["frst"]["kappa"], r["max_deg"], grading=wF.tolist(), **kw)
        mapped = {tuple(int(x) for x in T.T @ np.array(n)): v for n, v in gV.items()}
        same_cy = mapped == gF
        frozen = gV == {tuple(k): v for k, v in r["gvs_vex"]}
        ok = same_cy and frozen
        msg = f"vex vs FRST of the other polytope: {'agree' if same_cy else 'DIFFER'} on {len(set(mapped) | set(gF))} classes; vex vs stored: {'same' if frozen else 'DIFFER'}"
    fails += not ok
    print(f"{'PASS' if ok else 'FAIL'} {r['name']:18s} D={r['max_deg']:<3d} {msg} ({time.time() - t:.1f}s)", flush=True)
print("ALL PASS" if not fails else f"{fails} FAILURES")
sys.exit(1 if fails else 0)
