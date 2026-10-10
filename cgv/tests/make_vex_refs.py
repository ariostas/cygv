"""Regenerate tests/refs_vex (needs CYTools and regfans; the tests themselves do not).

Phase references: every fine regular fan (FRST and vex) of two small polytopes with vex phases, from Liam McAllister's
vexhkty note (h11 = 2: one FRST, one vex fan; h11 = 3: two FRSTs, eight vex fans), stored as (cones, q, kappa, grading,
max_deg) with the GVs cgv computed. Pair references: a vex phase and an FRST of a *different* polytope giving the same
CY (identified by T: c2_F = c2_V T, kappa_F = kappa_V(T,T,T)); vex_regress.py checks that the two GV tables agree
class by class (n_F = T^T n_V) under one grading, independently of the stored values.
Usage: make_vex_refs.py [pairs.jsonl]   (pairs.jsonl: vertex/cone records as written by bgv's collect_pairs.py)"""
import gzip, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import cgv_phase
from cytools import Polytope, Cone
from regfans import VectorConfiguration

OUT = os.path.join(HERE, "refs_vex"); os.makedirs(OUT, exist_ok=True)


def fans_of(p):
    labels = cgv_phase.ray_labels(p)
    vc = VectorConfiguration(np.array([p.points(which=l)[0] for l in labels]), labels=labels)
    _, tris, _ = vc.flip_graph(only_fine=True, only_regular=True)
    return [[tuple(int(x) for x in c) for c in t.cones()] for t in tris]


def save(name, obj):
    json.dump(obj, gzip.open(os.path.join(OUT, name + ".json.gz"), "wt"))
    print("wrote", name, flush=True)


for tag, V, D in (("liam_h2", [[-1, 1, -1, -1], [0, -1, -1, 0], [0, 0, -1, -1], [0, 0, 1, 0], [0, 1, 1, 1], [1, -1, -1, 0]], 20),
                  ("liam_h3", [[0, 0, 0, -1], [0, 1, 0, 1], [-1, -1, -1, -1], [-1, -1, 1, 0], [1, 0, -1, -1], [-1, -1, 0, -1], [1, 0, 0, 1]], 14)):
    p = Polytope(V)
    for i, cones in enumerate(fans_of(p)):
        ci, q, kappa = cgv_phase.phase_from_polytope(p, cones)
        d = cgv_phase.prepare(ci, q, kappa)
        gv, _ = cgv_phase.run(d, D)
        save(f"{tag}_fan{i}", dict(kind="phase", name=f"{tag}_fan{i}", polytope_vertices=V, cones=ci, q=q,
                                   kappa=[[a, b, c, v] for (a, b, c), v in kappa.items()], grading=d["grading_vector"],
                                   max_deg=D, vex_cones=d["vex_cones"], gvs=[[list(k), v] for k, v in sorted(gv.items())]))

if len(sys.argv) > 1:
    pairs = [json.loads(l) for l in open(sys.argv[1])]
    want = {4: 2, 5: 2}; D = 20
    for pr in pairs:
        if want.get(pr["h11"], 0) == 0: continue
        sides = {}
        for s in ("vex", "frst"):
            p = Polytope(np.array(pr[s]["vertices"]))
            lab = {tuple(int(x) for x in p.points(which=l)[0]): l for l in p.labels}
            ci, q, kappa = cgv_phase.phase_from_polytope(p, [tuple(lab[tuple(c)] for c in cone) for cone in pr[s]["cones"]])
            sides[s] = dict(cones=ci, q=q, kappa=kappa, d=cgv_phase.prepare(ci, q, kappa), vertices=pr[s]["vertices"])
        T = np.array(pr["T"][0])
        GF = (np.round(np.linalg.inv(T.T)).astype(np.int64) @ np.array(sides["frst"]["d"]["generators"]).T).T
        try:
            w = [int(x) for x in Cone(np.vstack([np.array(sides["vex"]["d"]["generators"]), GF])).find_grading_vector()]
        except Exception:
            continue
        gv = cgv_phase.compute_gv(sides["vex"]["cones"], sides["vex"]["q"], sides["vex"]["kappa"], D, grading=w)
        name = f"pair_h{pr['h11']}_{want[pr['h11']]}"
        save(name, dict(kind="pair", name=name, h11=pr["h11"], h21=pr["h21"], max_deg=D, T=T.tolist(), grading_vex=w,
                        vex=dict(polytope_vertices=sides["vex"]["vertices"], cones=sides["vex"]["cones"], q=sides["vex"]["q"],
                                 kappa=[[a, b, c, v] for (a, b, c), v in sides["vex"]["kappa"].items()]),
                        frst=dict(polytope_vertices=sides["frst"]["vertices"], cones=sides["frst"]["cones"], q=sides["frst"]["q"],
                                  kappa=[[a, b, c, v] for (a, b, c), v in sides["frst"]["kappa"].items()]),
                        gvs_vex=[[list(k), v] for k, v in sorted(gv.items())]))
        want[pr["h11"]] -= 1
        if not any(want.values()): break
