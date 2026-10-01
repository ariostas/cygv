"""Run many cgv jobs across devices and machines.

Each job is a cgv input file (from cgv_run.write_input; it carries the cone data, so a
remote host only needs the cgv binary). Slots (a local GPU, local CPU, remote hosts) pull
jobs from one queue: GPU slots take the largest remaining job, CPU slots the smallest.
Results go to OUT/<name>.txt (GVs) and OUT/<name>.json (slot, time, status). Finished
jobs are skipped on restart; a failed job is retried on a different slot.

    python dispatch.py slots.json OUT job1.txt job2.txt ...   # or: --jobs list.txt
    python dispatch.py slots.json --deploy                    # copy sources, build on remote hosts

slots.json: [{"name": "gpu0", "cmd": ["/path/to/cgv/cgv_gpu", "-q", "-t", "8", "-g", "0"], "kind": "gpu"},
             {"name": "local", "cmd": ["/path/to/cgv/cgv", "-q", "-t", "16"]},
             {"name": "node1", "ssh": "user@node1", "cmd": ["~/cgv/cgv", "-q", "-t", "32"]}]
Only list GPUs that do not drive a display.
"""
import argparse, json, os, shlex, subprocess, sys, threading, time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)


def job_size(path):
    """Rough cost estimate: the input's max degree relative to its smallest generator degree."""
    with open(path) as f:
        h11, ngen, D = (int(x) for x in f.readline().split()[:3])
    return D, os.path.getsize(path)


def run_slot(slot, path):
    cmd = slot["cmd"]
    if "ssh" in slot:
        remote = "f=$(mktemp); cat > $f; " + " ".join(cmd) + " $f; rc=$?; rm -f $f; exit $rc"
        argv = ["ssh", "-o", "BatchMode=yes", slot["ssh"], remote]
        stdin = open(path, "rb")
    else:
        argv = [os.path.expanduser(c) for c in cmd] + [path]
        stdin = subprocess.DEVNULL
    t0 = time.time()
    p = subprocess.run(argv, stdin=stdin, capture_output=True)
    return p.returncode, p.stdout, p.stderr.decode(errors="replace")[-2000:], time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slots"); ap.add_argument("out", nargs="?"); ap.add_argument("inputs", nargs="*")
    ap.add_argument("--jobs", help="file listing input paths, one per line")
    ap.add_argument("--deploy", action="store_true")
    ap.add_argument("--retries", type=int, default=2)
    a = ap.parse_args()
    slots = json.load(open(a.slots))

    if a.deploy:
        for host in sorted({s["ssh"] for s in slots if "ssh" in s}):
            print(f"deploy {host}", flush=True)
            files = [os.path.join(SRC, f) for f in ("gv.c", "main.c", "cgv.h", "Makefile")]
            subprocess.run(["ssh", host, "mkdir -p ~/cgv"], check=True)
            subprocess.run(["scp", "-q", *files, f"{host}:cgv/"], check=True)
            subprocess.run(["ssh", host, "cd ~/cgv && make cgv >/dev/null && echo built"], check=True)
        return

    paths = list(a.inputs)
    if a.jobs:
        paths += [l.strip() for l in open(a.jobs) if l.strip()]
    os.makedirs(a.out, exist_ok=True)
    name = lambda p: os.path.splitext(os.path.basename(p))[0]
    todo = [p for p in paths if not os.path.exists(os.path.join(a.out, name(p) + ".json"))
            or json.load(open(os.path.join(a.out, name(p) + ".json"))).get("status") != "ok"]
    todo.sort(key=job_size)  # smallest first; GPU slots pop from the end
    print(f"{len(paths)} jobs, {len(todo)} to run, {len(slots)} slots", flush=True)
    lock = threading.Lock()
    fails = {}  # path -> slots that failed it

    def take(slot):
        with lock:
            order = range(len(todo) - 1, -1, -1) if slot.get("kind") == "gpu" else range(len(todo))
            for i in order:
                if slot["name"] not in fails.get(todo[i], ()):
                    return todo.pop(i)
            return None

    def worker(slot):
        while True:
            p = take(slot)
            if p is None:
                return
            rc, out, err, dt = run_slot(slot, p)
            meta = {"input": p, "slot": slot["name"], "seconds": round(dt, 3), "returncode": rc}
            if rc == 0:
                with open(os.path.join(a.out, name(p) + ".txt"), "wb") as f:
                    f.write(out)
                meta.update(status="ok", n_gv=out.count(b"\n"))
            else:
                meta.update(status="failed", stderr=err)
                with lock:
                    fails.setdefault(p, set()).add(slot["name"])
                    if len(fails[p]) <= a.retries and len(fails[p]) < len(slots):
                        todo.append(p); todo.sort(key=job_size)
            json.dump(meta, open(os.path.join(a.out, name(p) + ".json"), "w"))
            print(f"{meta['status']:6s} {name(p):30s} {slot['name']:10s} {dt:9.2f}s", flush=True)

    th = [threading.Thread(target=worker, args=(s,)) for s in slots]
    for t in th: t.start()
    for t in th: t.join()


if __name__ == "__main__":
    main()
