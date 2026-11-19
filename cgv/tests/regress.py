"""Regression test: cgv vs frozen cygv references.

Usage: regress.py [--threads T] [--max-deg-cap D] [--filter S] [--bin PATH]
Exits nonzero on any mismatch.
"""
import argparse, glob, gzip, json, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import cgv_run

ap = argparse.ArgumentParser()
ap.add_argument("--threads", type=int, default=os.cpu_count())
ap.add_argument("--max-deg-cap", type=int, default=10**9, help="skip refs with max_deg above this")
ap.add_argument("--filter", default="")
ap.add_argument("--bin", default=None)
ap.add_argument("--extra", default="", help="extra cgv arguments, e.g. '-g 1'")
ap.add_argument("--no-cones", action="store_true", help="force the full-enumeration fallback")
a = ap.parse_args()
if a.bin:
    cgv_run.BIN = a.bin
fails = 0
for path in sorted(glob.glob(os.path.join(os.path.dirname(__file__), "refs", "*.json.gz"))):
    ref = json.load(gzip.open(path, "rt"))
    if a.filter not in ref["name"]:
        continue
    if ref["max_deg"] > a.max_deg_cap and ref["name"].startswith("2d7b127a"):
        continue
    want = {tuple(k): int(v) for k, v in ref["gvs"]}
    try:
        got, dt, err = cgv_run.run_cgv(ref["input"], ref["max_deg"], a.threads, extra=a.extra.split(), cones=not a.no_cones)
    except RuntimeError as e:
        fails += 1
        print(f"ERROR {ref['name']:31s} D={ref['max_deg']:<5d} {str(e).strip().splitlines()[-1]}", flush=True)
        continue
    ok = got == want
    fails += not ok
    extra = ""
    if not ok:
        bad = [k for k in set(want) | set(got) if want.get(k) != got.get(k)]
        extra = f" ndiff={len(bad)} e.g. {[(k, want.get(k), got.get(k)) for k in bad[:2]]}"
    print(f"{'PASS' if ok else 'FAIL'} {ref['name']:32s} D={ref['max_deg']:<5d} n={len(want):6d} {dt:7.2f}s{extra}", flush=True)
print("ALL PASS" if not fails else f"{fails} FAILURES")
sys.exit(1 if fails else 0)
