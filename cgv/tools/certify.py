"""Prove that a cgv output is exact, i.e. that its CRT modulus was large enough.

    python certify.py INPUT.txt GVS.txt NPRIMES [--threads N]

NPRIMES: the number of primes the GVs were lifted with (a run's "primes used: N", or
crt_combine's "N primes"); the primes are cgv's fixed sequence below 2^62.

cgv_maj bounds each |GV_C| from above, rounding upward, assuming every GV of lower degree
equals its candidate. If every bound is < M/2, then by induction on degree every candidate
is the true GV (a curve only feeds curves of higher degree). Exit 0 = certified."""
import argparse, math, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("input"); ap.add_argument("gvs"); ap.add_argument("nprimes", type=int)
ap.add_argument("--threads", type=int, default=os.cpu_count())
ap.add_argument("--prime0", type=int, default=0, help="index of the first prime (distributed runs)")
a = ap.parse_args()


def is_prime(n):
    if n % 2 == 0: return False
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for b in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(b, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True


primes, n = [], (1 << 62) - 1
while len(primes) < a.prime0 + a.nprimes:
    if is_prime(n): primes.append(n)
    n -= 2
logM = math.log2(math.prod(primes[a.prime0:]))
env = dict(os.environ, CGV_MAJ_CAND=os.path.abspath(a.gvs))
p = subprocess.run([os.path.join(os.path.dirname(HERE), "cgv_maj"), "-t", str(a.threads), a.input],
                   env=env, capture_output=True, text=True)
if p.returncode:
    sys.exit(f"cgv_maj failed: {p.stderr[-500:]}")
bits = {int(d): float(b) for d, b in (l.split() for l in p.stdout.splitlines())}
top = max(bits.values())
ok = top < logM - 1.01  # |GV| < M/2, with slack for log2's rounding
worst = max(bits, key=bits.get)
print(f"proven max |GV| < 2^{top:.2f} (degree {worst}); modulus 2^{logM:.2f} from {a.nprimes} primes: "
      + ("CERTIFIED" if ok else f"NOT certified; need {math.ceil((top + 1.01) / 61.99)} primes"))
sys.exit(0 if ok else 1)
