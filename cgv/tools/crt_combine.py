"""Combine residue files from distributed cgv runs (CGV_RESIDUES=1, CGV_PRIME0=k) into GVs.

    python crt_combine.py res_a.txt res_b.txt ... > gvs.txt

Each file holds residues of every GV mod its own primes (curves absent = residue 0). The
symmetric CRT lift uses all primes; the lift from all but the last prime must agree with
it for every curve (the same check the single-machine run makes), else exit status 2.
The modulus's size is printed to stderr for certification (tools/certify.py)."""
import sys, math

primes, rows = [], {}
for path in sys.argv[1:]:
    with open(path) as f:
        head = f.readline().split()
        assert head[:2] == ["#", "primes"], f"{path}: not a residue file"
        ps = [int(x) for x in head[2:]]
        off = len(primes)
        primes += ps
        for line in f:
            v = line.split()
            key = " ".join(v[: -len(ps)])
            r = rows.setdefault(key, {})
            for i, x in enumerate(v[-len(ps):]):
                r[off + i] = int(x)
if len(set(primes)) != len(primes):
    sys.exit("duplicate primes across files")
k = len(primes)
M = math.prod(primes); Mk = M // primes[-1]
# CRT basis: e_i = (M/p_i) * ((M/p_i)^-1 mod p_i); and the same without the last prime
E = [(M // p) * pow(M // p, -1, p) for p in primes]
Ek = [(Mk // p) * pow(Mk // p, -1, p) for p in primes[:-1]]
bad = 0
out = sys.stdout
for key, r in rows.items():
    x = sum(E[i] * r.get(i, 0) for i in range(k)) % M
    x = x - M if x > M // 2 else x
    y = sum(Ek[i] * r.get(i, 0) for i in range(k - 1)) % Mk
    y = y - Mk if y > Mk // 2 else y
    if x != y:
        bad += 1
    if x:
        out.write(f"{key} {x}\n")
print(f"crt_combine: {len(rows)} curves, {k} primes, log2 M = {math.log2(M):.2f}, unstable {bad}", file=sys.stderr)
sys.exit(2 if bad else 0)
