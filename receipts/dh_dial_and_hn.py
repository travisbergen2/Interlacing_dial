#!/usr/bin/env python3
"""
dh_dial_and_hn.py  (2026-09-23)

Receipts for the question "what alteration of the Davenport-Heilbronn (D-H) function
would flip its guaranteed off-line zeros to on-line zeros", plus the arithmetic hiding
in the Hadwiger-Nelson plane-coloring bounds.

Sections
  A. The D-H kappa dial: kappa from the root number of the odd character mod 5;
     functional-equation points on the real kappa-axis, Euler-product points at +-i;
     multiplicativity forces kappa^2 = -1; the log-derivative coefficients Lambda_f(n)
     are supported off the prime powers and grow; the mirror orbit at the certified zero.
  B. Function-field twin over F_q[T]: for conjugate characters chi, chi-bar mod an
     irreducible cubic, the real functional-equation-matched mixture (the D-H of the
     finite world) has its zeros on the circle |u| = q^{-1/2} iff the zeros of
     L(u,chi) and L(u,chi-bar) interlace; in degree 2 the criterion is sin(t1) sin(t2) <= 0.
     Higher degree: interlacing => on circle (one direction), tested by census.
  C. Plane-coloring arithmetic: the Moser / de Grey rotation angles are
     2 arcsin(1/(2R)) with cos = 1 - 1/(2 R^2); for R^2 = 3, 4, 16 they are (up to a
     half-turn and mirror) the angles pi/pi-bar of ideals of norm 3, 4, 16 in
     Q(sqrt-11), Q(sqrt-15), Q(sqrt-7).  The hexagonal 7-colouring is
     Z[omega] mod (3+omega) (norm 7); no Eisenstein ideal has norm 5 or 6; square tiles need 9.

Every printed number is recomputed here; asserts are the self-checks.
"""
import cmath, itertools, json, math, sys
import mpmath as mp

OUT = {}
def say(*a):
    print(*a); sys.stdout.flush()

# ----------------------------------------------------------------------------- A
mp.mp.dps = 30
I = mp.mpc(0, 1)
chi = {1: mp.mpc(1), 2: I, 3: -I, 4: mp.mpc(-1)}          # odd character mod 5, chi(2) = i

say("=" * 78)
say("A. THE DAVENPORT-HEILBRONN KAPPA DIAL")
say("=" * 78)

tau = mp.fsum(chi[n] * mp.exp(2 * mp.pi * I * n / 5) for n in range(1, 5))   # Gauss sum
eps = tau / (I * mp.sqrt(5))                                                # root number (odd: a = 1)
assert abs(abs(eps) - 1) < mp.mpf(10) ** -25
phi = mp.arg(eps)
kappa_plus = mp.tan(phi / 2)
kappa_DH = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)             # D-H's closed form
kappa_minus = -1 / kappa_plus
assert abs(kappa_plus - kappa_DH) < mp.mpf(10) ** -25
assert abs((1 + I * kappa_plus) / (1 - I * kappa_plus) - eps) < mp.mpf(10) ** -25
assert abs((1 + I * kappa_minus) / (1 - I * kappa_minus) + eps) < mp.mpf(10) ** -25
say(f"Gauss sum tau            = {mp.nstr(tau, 15)}   |tau|^2 = {mp.nstr(abs(tau)**2, 15)}")
say(f"root number eps          = {mp.nstr(eps, 15)}   arg eps = {mp.nstr(phi*180/mp.pi, 12)} deg")
say(f"kappa_+ = tan(arg eps/2) = {mp.nstr(kappa_plus, 20)}")
say(f"kappa_DH (closed form)   = {mp.nstr(kappa_DH, 20)}   |diff| = {mp.nstr(abs(kappa_plus-kappa_DH), 3)}")
say(f"kappa_- = -1/kappa_+     = {mp.nstr(kappa_minus, 20)}   (the antisymmetric functional equation)")
say("(1+i k)/(1-i k) = eps at kappa_+ and = -eps at kappa_-  [checked to 1e-25]")
OUT["eps"] = [float(mp.re(eps)), float(mp.im(eps))]
OUT["arg_eps_deg"] = float(phi * 180 / mp.pi)
OUT["kappa_plus"] = float(kappa_plus); OUT["kappa_minus"] = float(kappa_minus)

def f_kappa(s, k):
    """D-H family: 5^{-s}[zeta(s,1/5) + k zeta(s,2/5) - k zeta(s,3/5) - zeta(s,4/5)]"""
    return 5 ** (-s) * (mp.zeta(s, mp.mpf(1) / 5) + k * mp.zeta(s, mp.mpf(2) / 5)
                        - k * mp.zeta(s, mp.mpf(3) / 5) - mp.zeta(s, mp.mpf(4) / 5))
def Lam(s, k):
    return (5 / mp.pi) ** (s / 2) * mp.gamma((s + 1) / 2) * f_kappa(s, k)
def Lchi(s, conj=False):
    c = {n: (mp.conj(v) if conj else v) for n, v in chi.items()}
    return 5 ** (-s) * mp.fsum(c[a] * mp.zeta(s, mp.mpf(a) / 5) for a in range(1, 5))

say("\nA2. Functional-equation test  |Lam(s) - Lam(1-s)| (even) and |Lam(s) + Lam(1-s)| (odd)")
pts = [mp.mpc('0.3', '2.0'), mp.mpc('0.7', '17.5'), mp.mpc('1.4', '-3.3')]
fe = {}
for k, label in [(kappa_plus, 'kappa_+'), (kappa_minus, 'kappa_-'), (mp.mpf(0), 'kappa=0'), (mp.mpf(1), 'kappa=1')]:
    ev = max(abs(Lam(s, k) - Lam(1 - s, k)) for s in pts)
    od = max(abs(Lam(s, k) + Lam(1 - s, k)) for s in pts)
    fe[label] = (float(ev), float(od))
    say(f"  {label:8s}: even-defect {mp.nstr(ev, 3):>10s}   odd-defect {mp.nstr(od, 3):>10s}")
assert fe['kappa_+'][0] < 1e-20 and fe['kappa_-'][1] < 1e-20
assert fe['kappa=0'][0] > 1e-3 and fe['kappa=0'][1] > 1e-3
OUT["fe_defects"] = fe
# Euler-product points: kappa = +-i reproduce the single L-functions exactly
s = mp.mpc('0.6', '5.0')
assert abs(f_kappa(s, I) - Lchi(s)) < mp.mpf(10) ** -25 and abs(f_kappa(s, -I) - Lchi(s, True)) < mp.mpf(10) ** -25
say("  kappa = +i gives L(s,chi) exactly, kappa = -i gives L(s,chi-bar) exactly  [checked]")
say("  involution kappa -> -1/kappa: fixes +-i, swaps kappa_+ <-> kappa_-  [algebra: (1-ik)L+(1+ik)Lbar -> difference combination]")

say("\nA3. Multiplicativity: a(2)^2 - a(4) = kappa^2 + 1  -> zero only at kappa = +-i")
for k, label in [(kappa_plus, 'kappa_+'), (kappa_minus, 'kappa_-'), (I, '+i')]:
    say(f"  {label:8s}: a(2)^2 - a(4) = {mp.nstr(k*k + 1, 12)}")

# Lambda_f(n): the coefficients of -f'/f, by the recursion a(n) log n = sum_{d|n} Lambda_f(d) a(n/d)
def coeff_table(k, N):
    tbl = {0: 0, 1: 1, 2: k, 3: -k, 4: -1}
    return [None] + [tbl[n % 5] for n in range(1, N + 1)]
def lambda_f(a, N):
    divs = [[] for _ in range(N + 1)]
    for d in range(2, N + 1):
        for m in range(2 * d, N + 1, d):
            divs[m].append(d)
    L = [0] * (N + 1)
    for n in range(2, N + 1):
        s = a[n] * math.log(n)
        for d in divs[n]:            # proper divisors d with 1 < d < n
            s -= L[d] * a[n // d]
        L[n] = s
    return L
def is_prime_power(n):
    for p in range(2, n + 1):
        if n % p == 0:
            while n % p == 0: n //= p
            return n == 1
    return False
kf = complex(kappa_plus)
N1 = 40
Lf = lambda_f(coeff_table(kf, N1), N1)
Li = lambda_f(coeff_table(1j, N1), N1)      # Euler-product point: must be chi(n) Lambda(n)
say("\n    n   prime-power?   Lambda_f(n) at kappa_+       Lambda_f(n) at kappa=i (should be chi(n)Lambda(n))")
for n in range(2, 25):
    say(f"  {n:3d}   {'yes' if is_prime_power(n) else 'NO ':>4s}         {Lf[n].real:+.6f}                 {Li[n].real:+.6f}{Li[n].imag:+.6f}i")
for n in range(2, N1 + 1):
    if not is_prime_power(n):
        assert abs(Li[n]) < 1e-9, n          # supported on prime powers at the Euler-product point
    else:
        p = min(p for p in range(2, n + 1) if n % p == 0)
        assert abs(abs(Li[n]) - (math.log(p) if n % 5 else 0.0)) < 1e-9, n
say(f"  Lambda_f(6) at kappa_+ = {Lf[6].real:.6f} = (1 + kappa^2) log 6 = {(1 + kf.real**2) * math.log(6):.6f}")
assert abs(Lf[6].real - (1 + kf.real ** 2) * math.log(6)) < 1e-9
N2 = 4000
Lf2 = lambda_f(coeff_table(kf, N2), N2)
rows = []
for cut in (100, 400, 1000, 2000, 4000):
    mx = max(abs(Lf2[n]) / math.log(n) for n in range(2, cut + 1))
    rows.append((cut, mx))
say("  growth of max_{n<=N} |Lambda_f(n)|/log n at kappa_+ (for an Euler product this is <= 1):")
for cut, mx in rows: say(f"     N = {cut:5d}: {mx:12.4f}")
OUT["lambda_f_growth"] = rows
OUT["lambda_f_6"] = Lf[6].real

say("\nA4. Mirror orbit at the certified off-line zero (f_kappa_+ has real coefficients + functional equation)")
mp.mp.dps = 25
s0 = mp.mpc('0.808517182456637', '85.6993484853776')
orbit = {'s0': s0, 'conj s0': mp.conj(s0), '1 - s0': 1 - s0, '1 - conj s0': 1 - mp.conj(s0)}
for name, z in orbit.items():
    say(f"  |f({name:11s})| = {mp.nstr(abs(f_kappa(z, kappa_plus)), 3):>10s}    |L(chi)| = {mp.nstr(abs(Lchi(z)), 8)}   |L(chi-bar)| = {mp.nstr(abs(Lchi(z, True)), 8)}")
say("  -> at a zero of f the two constituents have EQUAL modulus and opposite matched phase; neither vanishes.")
OUT["orbit_fvals"] = {k: float(abs(f_kappa(z, kappa_plus))) for k, z in orbit.items()}

# ----------------------------------------------------------------------------- B
say("\n" + "=" * 78)
say("B. FUNCTION-FIELD TWIN: real matched mixtures of conjugate L-functions over F_q[T]")
say("=" * 78)

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
    """F_q[T]/(m), m monic of degree d given low->high (m[d] = 1)."""
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
            # order exactly N = q^d - 1 is possible only in a field (unit group of a non-field is smaller)
            if self.pow(g, N) == self.one and all(self.pow(g, N // l) != self.one for l in pf):
                return g
        return None          # ring is not a field (m reducible)

def find_irreducible(q, d):
    for tail in itertools.product(range(q), repeat=d):
        m = list(tail) + [1]
        if m[0] == 0: continue
        R = Ring(q, m)
        if R.generator() is not None:
            return m
    raise RuntimeError

def poly_roots(coeffs_low_to_high):
    cs = [complex(c) for c in coeffs_low_to_high]
    while len(cs) > 1 and abs(cs[-1]) < 1e-14: cs.pop()
    if len(cs) <= 1: return []
    mp.mp.dps = 20
    rts = mp.polyroots(list(reversed(cs)), maxsteps=200, extraprec=60)
    return [complex(r) for r in rts]

def interlace(anglesA, anglesB):
    """Do two equal-size multisets of angles on the circle alternate (cyclically)?"""
    pts = sorted([(a % (2 * math.pi), 0) for a in anglesA] + [(b % (2 * math.pi), 1) for b in anglesB])
    labs = [l for _, l in pts]
    return all(labs[i] != labs[(i + 1) % len(labs)] for i in range(len(labs)))

def ff_census(q, d, verbose):
    m = find_irreducible(q, d)
    R = Ring(q, m); N = q ** d - 1; g = R.generator()
    logt, x = {}, R.one
    for k in range(N):
        logt[x] = k; x = R.mul(x, g)
    assert len(logt) == N
    def monics(deg):
        for cs in itertools.product(range(q), repeat=deg):
            yield tuple(list(cs) + [1] + [0] * (d - 1 - deg))
    constants = [tuple([c] + [0] * (d - 1)) for c in range(1, q)]
    def chi_val(j, f): return cmath.exp(2j * math.pi * j * logt[f] / N)
    res = {'q': q, 'd': d, 'm': m, 'N': N, 'pairs': 0, 'on': 0, 'off': 0,
           'interlaced_on': 0, 'interlaced_off': 0, 'nonint_on': 0, 'nonint_off': 0, 'crit_mismatch': 0}
    say(f"\n  q = {q}, modulus m = {m} (low->high), unit group order {N}, generator {g}")
    Ltab = {}
    for j in range(1, N):
        S = [sum(chi_val(j, f) for f in monics(deg)) for deg in range(d)]
        odd = any(abs(chi_val(j, c) - 1) > 1e-9 for c in constants)
        Ltab[j] = (S, odd)
        rts = poly_roots(S)
        if odd:
            assert len(rts) == d - 1
            for r in rts: assert abs(abs(r) - q ** -0.5) < 1e-8, (q, d, j, r)   # Weil's RH, odd chi
        else:
            assert len(rts) == d - 1
            mods = sorted(abs(r) for r in rts)
            assert abs(mods[-1] - 1) < 1e-8 and all(abs(mm - q ** -0.5) < 1e-8 for mm in mods[:-1])  # trivial zero u=1 + RH
    say(f"  Weil's RH checked for all {N-1} nontrivial characters: odd ones have all {d-1} roots on |u| = q^-1/2; even ones have (1-u) times roots on the circle.")
    n = d - 1
    seen = set()
    for j in range(1, N):
        jb = (N - j) % N
        if j in seen or jb in seen or j == jb: continue
        S, odd = Ltab[j]
        if not odd: continue
        seen.update({j, jb})
        Sb = Ltab[jb][0]
        s = [S[k] * q ** (-k / 2) for k in range(d)]        # L~(z) = L(z/sqrt q), z = sqrt(q) u
        sb = [Sb[k] * q ** (-k / 2) for k in range(d)]
        assert abs(abs(s[n]) - 1) < 1e-9
        for k in range(d):                                   # functional equation of L~: s_{n-k} = s_n conj(s_k)
            assert abs(s[n - k] - s[n] * s[k].conjugate()) < 1e-8
        w = s[n] - 1
        c = 1j * w.conjugate() if abs(w) > 1e-9 else 1.0     # even (symmetric) matched mixture
        F = [2 * (c * s[k]).real for k in range(d)]
        # 2Re(c s_k) == c s_k + conj(c) sb_k  because sb_k = conj(s_k): the mixture has real coefficients
        for k in range(d):
            assert abs(sb[k] - s[k].conjugate()) < 1e-9
            assert abs(F[k] - F[n - k]) < 1e-8, "even mixture must be palindromic"
        rts = poly_roots(F)
        on = all(abs(abs(r) - 1) < 1e-7 for r in rts)
        thA = [(-cmath.phase(r)) for r in poly_roots(s)]     # roots of L~(chi) are e^{-i theta}
        rootsA = [cmath.phase(r) for r in poly_roots(s)]
        rootsB = [cmath.phase(r) for r in poly_roots(sb)]
        il = interlace(rootsA, rootsB)
        res['pairs'] += 1
        res['on' if on else 'off'] += 1
        key = ('interlaced_' if il else 'nonint_') + ('on' if on else 'off')
        res[key] += 1
        if n == 2:
            t1, t2 = thA
            crit = math.sin(t1) * math.sin(t2) <= 1e-12
            if crit != on: res['crit_mismatch'] += 1
            if verbose:
                say(f"    pair (chi_{j:>3d}, chi_{jb:>3d}): zero angles of L(chi) {math.degrees(-t1):8.2f}, {math.degrees(-t2):8.2f} deg | "
                    f"sin*sin = {math.sin(t1)*math.sin(t2):+.3f} | interlaced {str(il):5s} | mixture on circle {str(on):5s} | |roots| = "
                    + ", ".join(f"{abs(r):.4f}" for r in rts))
    return res

results_B = []
for (q, d, verbose) in [(3, 3, True), (5, 3, False), (7, 3, False), (3, 4, False), (5, 4, False), (3, 5, False)]:
    r = ff_census(q, d, verbose)
    results_B.append(r)
    say(f"  census q={q} d={d} (L-degree {d-1}): {r['pairs']} conjugate odd pairs | on-circle {r['on']} | off-circle {r['off']} | "
        f"interlaced->on {r['interlaced_on']} | interlaced->off {r['interlaced_off']} | non-interlaced->on {r['nonint_on']} | non-interlaced->off {r['nonint_off']}"
        + (f" | degree-2 criterion mismatches {r['crit_mismatch']}" if d == 3 else ""))
    assert r['interlaced_off'] == 0, "Hermite-Biehler direction violated"
    if d == 3:
        assert r['crit_mismatch'] == 0 and r['nonint_on'] == 0, "degree-2 iff violated"
OUT["ff_census"] = results_B

# ----------------------------------------------------------------------------- C
say("\n" + "=" * 78)
say("C. PLANE-COLOURING ARITHMETIC")
say("=" * 78)
say("C1. Rotation angles 2 arcsin(1/(2R)) move a point at distance R by exactly 1: cos = 1 - 1/(2R^2)")
def rot(R2):
    c = 1 - 1 / (2 * R2); sn = math.sqrt(1 - c * c)
    return complex(c, sn)
w_moser = rot(3); w_4 = rot(4); w_16 = rot(16)
say(f"  R^2 = 3  (Moser spindle):  cos = 5/6  = {w_moser.real:.12f}, e^(i theta) = (5 + sqrt(-11))/6")
say(f"  R^2 = 4  (de Grey):        cos = 7/8  = {w_4.real:.12f}, e^(i theta) = (7 + sqrt(-15))/8")
say(f"  R^2 = 16 (de Grey):        cos = 31/32= {w_16.real:.12f}, e^(i theta) = (31 + 3 sqrt(-7))/32")
# Hecke-angle identifications: e^{i theta} = -conj(pi/pi-bar) for an element pi of the named norm
alpha = complex(0.5, math.sqrt(11) / 2)        # (1+sqrt-11)/2, norm 3
beta = complex(0.5, math.sqrt(15) / 2)         # (1+sqrt-15)/2, norm 4
gamma = complex(0.5, math.sqrt(7) / 2)         # (1+sqrt-7)/2,  norm 2
def spin(z): return z / z.conjugate()
assert abs(abs(alpha) ** 2 - 3) < 1e-12 and abs(abs(beta) ** 2 - 4) < 1e-12 and abs(abs(gamma) ** 2 - 2) < 1e-12
assert abs(w_moser - (-spin(alpha).conjugate())) < 1e-12
assert abs(w_4 - (-spin(beta).conjugate())) < 1e-12
assert abs(w_16 - (-(spin(gamma) ** 4))) < 1e-12
say("  identities checked to 1e-12:")
say("    (5+sqrt-11)/6  = -conj(alpha/alpha-bar),  alpha = (1+sqrt-11)/2, norm 3   [3 splits in Q(sqrt-11)]")
say("    (7+sqrt-15)/8  = -conj(beta/beta-bar),    beta  = (1+sqrt-15)/2, norm 4   [(beta) = p^2, p above 2 in Q(sqrt-15)]")
say("    (31+3sqrt-7)/32= -(gamma/gamma-bar)^4,    gamma = (1+sqrt-7)/2,  norm 2   [2 splits in Q(sqrt-7)]")
say("  i.e. each rotation is, up to a half-turn (a lattice symmetry) and a mirror, the angle pi/pi-bar of an ideal of norm R^2.")
OUT["rotations"] = {"moser_cos": w_moser.real, "degrey_cos": [w_4.real, w_16.real]}

say("\nC2. The hexagonal 7-colouring is Z[omega] mod (3 + omega), an ideal of norm 7")
omega = cmath.exp(2j * math.pi / 3)
def norm_eis(a, b): return a * a - a * b + b * b
norms = sorted({norm_eis(a, b) for a in range(-6, 7) for b in range(-6, 7)})
say(f"  Eisenstein norms that occur: {norms[:12]} ...  (5 and 6 do NOT occur: 2 is inert)")
assert 5 not in norms and 6 not in norms and 7 in norms
# colouring c(a + b omega) = (a - 3b) mod 7 is the quotient map (omega == -3 mod (3+omega), since 3 + omega == 0)
pts = [(a, b) for a in range(-8, 9) for b in range(-8, 9)]
def col(a, b): return (a - 3 * b) % 7
# same colour <=> difference in the ideal (3+omega); minimal nonzero norm of the ideal:
min_same = min(norm_eis(a - a2, b - b2) for (a, b) in pts for (a2, b2) in pts
               if (a, b) != (a2, b2) and col(a, b) == col(a2, b2))
say(f"  minimal squared lattice distance between same-coloured centres = {min_same}  (= N(3+omega))")
assert min_same == 7
# hexagons of circumradius r tile with centre spacing sqrt(3) r; same-colour centres >= sqrt(21) r apart
r_lo, r_hi = 1 / (math.sqrt(21) - 2), 0.5
say(f"  admissible hexagon circumradius r in ({r_lo:.4f}, {r_hi:.4f}): diameter 2r < 1 and same-colour tiles >= sqrt(21) r - 2r > 1")
say(f"  e.g. r = 0.45: diameter 0.90 < 1, same-colour separation >= {math.sqrt(21)*0.45 - 0.9:.4f} > 1  -> a valid 7-colouring")
# general lattice-ideal construction: need sqrt(3 N) r - 2 r > 1 with 2 r < 1  ->  N > 16/3
say(f"  a norm-N Eisenstein ideal works iff sqrt(3N) > 4, i.e. N > {16/3:.3f}: N = 7 is the first norm that occurs -> 7 colours")
# square tiles on Z[i]: side t, same-colour centres sqrt(N) t apart, tile diameter t sqrt 2:  need sqrt(N) t - sqrt2 t > 1, sqrt2 t < 1 -> N > 8 -> N = 9 (3 inert)
gnorms = sorted({a * a + b * b for a in range(-6, 7) for b in range(-6, 7)})
say(f"  square tiles on Z[i]: need N > 8; Gaussian norms {gnorms[:10]} -> N = 9 = 3^2 (3 inert): 9 colours, the classical square-grid count")
say("  lattice unit-distance graphs themselves: Z[i] mod (1+i) -> 2 colours; Z[omega] mod (1-omega) -> 3 colours.")

with open('/agent/workspace/dh_dial/dh_dial_and_hn.json', 'w') as fh:
    json.dump(OUT, fh, indent=1, default=str)
say("\nALL ASSERTS PASSED; json written.")
