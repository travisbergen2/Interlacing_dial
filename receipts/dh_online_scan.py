#!/usr/bin/env python3
"""
dh_online_scan.py  (2026-09-23)

On the critical line the Davenport-Heilbronn even function is, up to a constant phase,
Z_f(t) = Z_chi(t) + Z_chi(-t), where Z_chi(t) = eps^{-1/2} Lambda(1/2 + it, chi) is the
real Hardy-type function of the odd character chi mod 5 (chi(2) = i).  So the D-H function
is the t -> -t symmetrisation of a single complex L-function.  This scan, in a window
around the certified off-line zero 0.808517 + 85.699348 i, lists

  * the zeros of Z_chi(t)   (= ordinates of zeros of L(s,chi) at height +t),
  * the zeros of Z_chi(-t)  (= ordinates of zeros of L(s,chi-bar) at height +t),
  * whether the two lists interlace, and where they fail,
  * the on-line zeros of the even D-H function Z_f = Z_chi(t) + Z_chi(-t),
  * |Z_chi(t)| vs |Z_chi(-t)| at the off-line ordinate (the equal-amplitude locus).

Numerics only (mpmath, 18 digits); nothing here is a theorem.
"""
import json, sys
import mpmath as mp
mp.mp.dps = 18
I = mp.mpc(0, 1)
chi = {1: mp.mpc(1), 2: I, 3: -I, 4: mp.mpc(-1)}
tau = mp.fsum(chi[n] * mp.exp(2 * mp.pi * I * n / 5) for n in range(1, 5))
eps = tau / (I * mp.sqrt(5))
sq = mp.sqrt(eps)

def Lchi(s):
    return 5 ** (-s) * mp.fsum(chi[a] * mp.zeta(s, mp.mpf(a) / 5) for a in range(1, 5))
def Z(t):
    s = mp.mpc(0.5, t)
    v = (5 / mp.pi) ** (s / 2) * mp.gamma((s + 1) / 2) * Lchi(s) / sq
    return v

# reality check of Z (validates the root number)
worst = 0
for t in (3.7, 20.1, 85.7, -61.3):
    v = Z(t); worst = max(worst, abs(mp.im(v)) / abs(v))
print(f"max |Im Z|/|Z| at test points = {mp.nstr(worst, 3)}  (must be ~1e-15)"); sys.stdout.flush()
assert worst < 1e-10

T0 = mp.mpf(sys.argv[1]) if len(sys.argv) > 1 else mp.mpf(78)
T1 = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf(94)
OFF = sys.argv[3] if len(sys.argv) > 3 else '85.6993484853776'     # off-line ordinate to inspect
TAG = sys.argv[4] if len(sys.argv) > 4 else ''
h = mp.mpf('0.02')
grid = []
t = T0
while t <= T1 + h / 2:
    zp = mp.re(Z(t)); zm = mp.re(Z(-t))
    grid.append((t, zp, zm))
    t += h
print(f"grid of {len(grid)} points computed"); sys.stdout.flush()

def refine(fn, a, b):
    fa = fn(a)
    for _ in range(60):
        c = (a + b) / 2; fc = fn(c)
        if fa * fc <= 0: b = c
        else: a, fa = c, fc
    return (a + b) / 2

zp_zeros, zm_zeros, zf_zeros, zg_zeros = [], [], [], []
for (t, zp, zm), (t2, zp2, zm2) in zip(grid, grid[1:]):
    if zp * zp2 < 0: zp_zeros.append(refine(lambda x: mp.re(Z(x)), t, t2))
    if zm * zm2 < 0: zm_zeros.append(refine(lambda x: mp.re(Z(-x)), t, t2))
    if (zp + zm) * (zp2 + zm2) < 0: zf_zeros.append(refine(lambda x: mp.re(Z(x)) + mp.re(Z(-x)), t, t2))
    if (zp - zm) * (zp2 - zm2) < 0: zg_zeros.append(refine(lambda x: mp.re(Z(x)) - mp.re(Z(-x)), t, t2))

merged = sorted([(float(z), 'chi   ') for z in zp_zeros] + [(float(z), 'chibar') for z in zm_zeros])
print("\nordinate      which L-function vanishes there      interlacing")
prev = None; breaks = []
for z, lab in merged:
    flag = ''
    if prev is not None and prev[1] == lab:
        flag = '<-- two in a row: interlacing FAILS between %.3f and %.3f' % (prev[0], z)
        breaks.append((prev[0], z))
    print(f"{z:10.4f}   {lab}                                   {flag}")
    prev = (z, lab)
print(f"\ninterlacing failures in [{float(T0)}, {float(T1)}]: {breaks}")
print(f"on-line zeros of the EVEN D-H function (Z_chi(t)+Z_chi(-t)): {[round(float(z),4) for z in zf_zeros]}")
print(f"on-line zeros of the ODD  D-H function (Z_chi(t)-Z_chi(-t)): {[round(float(z),4) for z in zg_zeros]}")
t_off = mp.mpf(OFF)
a, b = abs(Z(t_off)), abs(Z(-t_off))
print(f"\nat the off-line ordinate t = {t_off}: |Z_chi(t)| = {mp.nstr(a, 6)}, |Z_chi(-t)| = {mp.nstr(b, 6)}, ratio = {mp.nstr(a/b, 6)}")
# where does |Z_chi(t)|/|Z_chi(-t)| cross 1 in the window?
cross = []
for (t, zp, zm), (t2, zp2, zm2) in zip(grid, grid[1:]):
    d1, d2 = abs(zp) - abs(zm), abs(zp2) - abs(zm2)
    if d1 * d2 < 0: cross.append(round(float((t + t2) / 2), 3))
print(f"equal-amplitude crossings |Z_chi(t)| = |Z_chi(-t)| in the window: {cross}")
json.dump({'chi_zeros': [float(z) for z in zp_zeros], 'chibar_zeros': [float(z) for z in zm_zeros],
           'even_DH_online_zeros': [float(z) for z in zf_zeros], 'odd_DH_online_zeros': [float(z) for z in zg_zeros],
           'interlacing_breaks': breaks, 'equal_amplitude_crossings': cross,
           'ratio_at_offline_ordinate': float(a / b)},
          open(f'/agent/workspace/dh_dial/dh_online_scan{TAG}.json', 'w'), indent=1)
print("done")
