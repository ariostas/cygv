"""Freeze cygv 0.2.3 outputs as regression references (tests/refs/*.json.gz).

Each ref: {"name", "max_deg", "input": <cygv-style input dict>, "gvs": [[curve, "gv"], ...]}
Existing cygv outputs for the hard case (data/2d7b127a/cygv_deg*.json) are reused.
"""
import gzip, json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from cgv_run import run_cygv

OUT = os.path.join(os.path.dirname(__file__), "refs")
G = os.environ.get("GVS_DATA", os.path.expanduser("~/gvs"))  # where the reference inputs live (only for regenerating refs)


def save(name, d, D, gvs):
    ref = dict(name=name, max_deg=D, input=d, gvs=[[list(k), str(v)] for k, v in sorted(gvs.items())])
    with gzip.open(f"{OUT}/{name}.json.gz", "wt") as f:
        json.dump(ref, f)
    print(name, D, len(gvs), flush=True)


for d in json.load(open(f"{G}/bench/inputs_quintic.json")):
    save("quintic_D10", d, 10, run_cygv(d, 10)[0])
for d in json.load(open(f"{G}/bench/inputs_h10.json")):
    degs = np.array(d["generators"]) @ np.array(d["grading_vector"])
    for K in (5, 7):
        D = int(K * degs[degs > 0].min())
        save(f"{d['name']}_{d['grading']}_K{K}", d, D, run_cygv(d, D)[0])
hard = json.load(open(f"{G}/data/2d7b127a/inputs.json"))[0]
for D in (10, 12, 14, 16, 18):
    gvs = {tuple(v): int(g) for v, g in json.load(open(f"{G}/data/2d7b127a/cygv_deg{D}.json"))}
    save(f"2d7b127a_D{D}", hard, D, gvs)
