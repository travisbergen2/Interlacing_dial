#!/usr/bin/env python3
"""
interlacing_dial.py — Python twin of Room IX, "The Interlacing Dial" (2026-09-24).

The function-field Davenport–Heilbronn twin.  Over F_q[T] with an irreducible cubic
modulus m, every odd Dirichlet character chi (nontrivial on the constants F_q^x) has

    L(u, chi) = 1 + S1 u + S2 u^2,      S1 = sum_a chi(T + a),   S2 = sum_{a,b} chi(T^2 + aT + b),

and Weil's theorem puts both roots on |u| = q^{-1/2}.  Normalise z = sqrt(q) u:

    L~(z) = 1 + s1 z + s2 z^2,   s1 = S1/sqrt q,   s2 = S2/q,   |s2| = 1 (root number),
    functional equation  s1 = s2 conj(s1)  =>  s1 = r e^{i phi/2} with r REAL, phi = arg s2.

The matched real mixture  F = c L~(chi) + conj(c) L~(chi-bar)  with  c/conj(c) = 1/s2
(exactly Davenport–Heilbronn's recipe: (1+i kappa)/(1-i kappa) = root number) is

    F(z) = 4 sin(phi/2) [ cos(phi/2) (1 + z^2) + r z ],

so its zeros are on the unit circle  iff  |r| <= 2 |cos(phi/2)|,
and that inequality is exactly the statement that the zeros of L~(chi) and L~(chi-bar)
INTERLACE on the circle (angles -phi/2 +- psi and +phi/2 -+ psi, cos psi = -r/2).

This script recomputes everything from the field arithmetic (no closed forms are
trusted): finds m and a generator g, tabulates discrete logs, sums the characters,
finds the roots numerically, tests interlacing by sorting angles, tests the mixture's
roots by the quadratic formula, and only THEN checks the closed form against the numerics.
Output: interlacing_dial.json (embedded in index.html as the self-check target) and a
human-readable census.  Nothing here bears on the Riemann Hypothesis.
"""
import cmath, itertools, json, math, os, sys

def prime_factors(n):
    out, p = [], 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: out.append(n)
    return out

class Ring:
    """F_q[T]/(m); elements are coefficient tuples low->high of length d; m monic (m[d] = 1)."""
    def __init__(self, q, m):
        self.q, self.m, self.d = q, m, len(m) - 1
        self.one = tuple([1] + [0] * (self.d - 1))
    def mul(self, x, y):
        q, d, m = self.q, self.d, self.m
        prod = [0] * (2 * d - 1)
        for i, xi in enumerate(x):
            if xi:
                for j, yj in enumerate(y):
                    if yj: prod[i + j] = (prod[i + j] + xi * yj) % q
        for k in range(2 * d - 2, d - 1, -1):
            c = prod[k]
            if c:
                for i in range(d):
                    prod[k - d + i] = (prod[k - d + i] - c * m[i]) % q
                prod[k] = 0
        return tuple(prod[:d])
    def pow(self, x, e):
        r, b = self.one, x
        while e:
            if e & 1: r = self.mul(r, b)
            b = self.mul(b, b); e >>= 1
        return r
    def generator(self):
        N = self.q ** self.d - 1
        pf = prime_factors(N)
        for g in itertools.product(range(self.q), repeat=self.d):
            if not any(g): continue
            if self.pow(g, N) == self.one and all(self.pow(g, N // l) != self.one for l in pf):
                return g
        return None

def find_irreducible(q, d):
    for tail in itertools.product(range(q), repeat=d):
        m = list(tail) + [1]
        if m[0] == 0: continue
        if Ring(q, m).generator() is not None:
            return m
    raise RuntimeError("no irreducible found")

def quad_roots(a, b, c):
    """roots of a z^2 + b z + c (complex coefficients)"""
    disc = cmath.sqrt(b * b - 4 * a * c)
    return [(-b + disc) / (2 * a), (-b - disc) / (2 * a)]

def interlaced(angA, angB):
    pts = sorted([(a % (2 * math.pi), 0) for a in angA] + [(b % (2 * math.pi), 1) for b in angB])
    lab = [l for _, l in pts]
    return all(lab[i] != lab[(i + 1) % len(lab)] for i in range(len(lab)))

def census(q, d=3):
    m = find_irreducible(q, d)
    R = Ring(q, m); N = q ** d - 1; g = R.generator()
    log = {}; x = R.one
    for k in range(N):
        log[x] = k; x = R.mul(x, g)
    assert len(log) == N and x == R.one
    monic1 = [(a, 1, 0) for a in range(q)]                       # T + a
    monic2 = [(b, a, 1) for a in range(q) for b in range(q)]     # T^2 + a T + b
    consts = [(c, 0, 0) for c in range(1, q)]
    def chi(j, f): return cmath.exp(2j * math.pi * j * log[f] / N)
    rows = []; seen = set()
    n_odd = 0; weil_ok = 0
    for j in range(1, N):
        odd = any(abs(chi(j, c) - 1) > 1e-9 for c in consts)
        if not odd: continue
        n_odd += 1
        S1 = sum(chi(j, f) for f in monic1); S2 = sum(chi(j, f) for f in monic2)
        rootsA = quad_roots(S2 / q, S1 / math.sqrt(q), 1)        # roots of L~(z) = 1 + s1 z + s2 z^2
        if all(abs(abs(z) - 1) < 1e-9 for z in rootsA): weil_ok += 1
        jb = N - j
        if j in seen or jb in seen or j == jb: continue      # j == N/2 is the real (self-conjugate) character: no pair
        seen.update({j, jb})
        s1, s2 = S1 / math.sqrt(q), S2 / q
        assert abs(abs(s2) - 1) < 1e-9, "root number must have modulus 1"
        assert abs(s1 - s2 * s1.conjugate()) < 1e-9, "functional equation s1 = s2 conj(s1)"
        phi = cmath.phase(s2)
        r = (s1 * cmath.exp(-1j * phi / 2)).real
        assert abs((s1 * cmath.exp(-1j * phi / 2)).imag) < 1e-9, "s1 e^{-i phi/2} must be real"
        # the matched mixture, built from the definition (no closed form)
        w = s2 - 1
        c = 1j * w.conjugate() if abs(w) > 1e-9 else 1.0
        F0 = 2 * (c * 1).real; F1 = 2 * (c * s1).real; F2 = 2 * (c * s2).real
        assert abs(F0 - F2) < 1e-9, "mixture must be palindromic"
        rootsF = quad_roots(F2, F1, F0)
        on = all(abs(abs(z) - 1) < 1e-7 for z in rootsF)
        rootsB = [z.conjugate() for z in rootsA]                  # roots of L~(chi-bar)
        il = interlaced([cmath.phase(z) for z in rootsA], [cmath.phase(z) for z in rootsB])
        closed = abs(r) <= 2 * abs(math.cos(phi / 2)) + 1e-12
        assert on == il == closed, (q, j, on, il, closed)
        rows.append({"j": j, "jb": jb, "S1": [S1.real, S1.imag], "S2": [S2.real, S2.imag],
                     "r": r, "phi_deg": math.degrees(phi),
                     "rootsA_deg": sorted(math.degrees(cmath.phase(z)) for z in rootsA),
                     "mixture_root_moduli": sorted(abs(z) for z in rootsF), "on": on})
    assert weil_ok == n_odd, "Weil's RH failed for some odd character"
    on = sum(1 for r_ in rows if r_["on"])
    return {"q": q, "d": d, "m_low_to_high": m, "generator": list(g), "N": N,
            "odd_characters": n_odd, "weil_verified": weil_ok,
            "pairs": len(rows), "on": on, "off": len(rows) - on, "rows": rows}

if __name__ == "__main__":
    out = {"generated": "2026-09-24", "census": [], "totals": {}}
    tot_pairs = tot_on = 0
    for q in (3, 5, 7):
        c = census(q)
        out["census"].append(c)
        tot_pairs += c["pairs"]; tot_on += c["on"]
        print(f"q = {q}: m = {c['m_low_to_high']} (low->high), g = {c['generator']}, N = {c['N']}, "
              f"odd characters {c['odd_characters']} (Weil verified {c['weil_verified']}), "
              f"conjugate pairs {c['pairs']}: on circle {c['on']}, off circle {c['off']}")
        if q == 3:
            for r_ in c["rows"]:
                print(f"    (chi_{r_['j']:>2d}, chi_{r_['jb']:>2d})  r = {r_['r']:+.6f}  phi = {r_['phi_deg']:+8.3f} deg  "
                      f"2|cos(phi/2)| = {2*abs(math.cos(math.radians(r_['phi_deg'])/2)):.6f}  "
                      f"L(chi) root angles {r_['rootsA_deg'][0]:+8.2f}, {r_['rootsA_deg'][1]:+8.2f}  "
                      f"mixture |roots| {r_['mixture_root_moduli'][0]:.4f}, {r_['mixture_root_moduli'][1]:.4f}  -> {'ON ' if r_['on'] else 'OFF'}")
    out["totals"] = {"pairs": tot_pairs, "on": tot_on, "off": tot_pairs - tot_on}
    print(f"TOTAL: {tot_pairs} pairs, {tot_on} on the circle, {tot_pairs - tot_on} off; "
          f"criterion |r| <= 2|cos(phi/2)| == interlacing == on-circle held for every pair (asserted).")
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "interlacing_dial.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote interlacing_dial.json")
