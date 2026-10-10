"""Repeat reference cases with varying thread counts / modes to catch races.
Usage: stress.py N_ITER   (logs mismatches; exits nonzero if any)"""
import gzip, json, os, random, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import cgv_run
HERE = os.path.dirname(os.path.abspath(__file__))
refs = {n: json.load(gzip.open(f"{HERE}/refs/{n}.json.gz", "rt")) for n in ("2d7b127a_D20", "2d7b127a_D16", "h10_44_1857_0_default_K7")}
rng = random.Random(int(time.time()))
bad = 0
for it in range(int(sys.argv[1])):
    name = rng.choice(list(refs))
    ref = refs[name]
    want = {tuple(k): int(v) for k, v in ref["gvs"]}
    mode = rng.choice(os.environ.get("STRESS_MODES", "cpu,cpu,cpu_cs,gpu").split(","))
    thr = rng.choice([1, 2, 3, 5, 8, 13, 24])
    env = dict(os.environ)
    extra = []
    if mode == "cpu_cs":
        env["CGV_CS_NU"] = "0"; env["CGV_LSYNC"] = "100000"
    if mode == "gpu":
        cgv_run.BIN = os.path.join(HERE, "..", "cgv_gpu"); extra = ["-g", "0"]; env["CGV_GPU_MIN"] = "0"
    else:
        cgv_run.BIN = os.path.join(HERE, "..", "cgv")
    os.environ.clear(); os.environ.update(env)
    try:
        got, dt, err = cgv_run.run_cgv(ref["input"], ref["max_deg"], thr, extra=extra)
        ok = got == want
    except RuntimeError as e:
        ok, dt = None, 0
        err_text = str(e).strip().splitlines()[-1] if str(e).strip() else "?"
    if ok is None:
        bad += 1
        print(f"ERROR it={it} {name} mode={mode} threads={thr}: {err_text}", flush=True)
    elif not ok:
        bad += 1
        print(f"MISMATCH it={it} {name} mode={mode} threads={thr}", flush=True)
    elif it % 20 == 0:
        print(f"ok it={it} {name} mode={mode} threads={thr} {dt:.2f}s", flush=True)
print("STRESS DONE", "bad =", bad)
sys.exit(1 if bad else 0)
