# Room IX · The Interlacing Dial

**Live:** https://fractalyouniverse.org/Interlacing_dial/ (also https://travisbergen2.github.io/Interlacing_dial/)

The Davenport–Heilbronn mixture in the finite world. In 1936 Davenport and Heilbronn mixed the two
conjugate odd Dirichlet L-functions mod 5 into one function that keeps the functional equation and loses
the Euler product, and proved its zeros leave the critical line. This room builds the same mixture over
F_q[T] — where the Riemann Hypothesis is Weil's theorem — and shows that there the mixture's fate is
decided by one inequality:

> For an odd character χ mod an irreducible cubic, with L̃(z) = 1 + r·e^{iφ/2} z + e^{iφ} z² (Weil: |r| ≤ 2),
> the matched real mixture c·L̃(χ) + c̄·L̃(χ̄), c/c̄ = 1/e^{iφ} — Davenport–Heilbronn's own recipe — equals
> 4 sin(φ/2)·[cos(φ/2)(1 + z²) + r z]. Its zeros lie on the unit circle **iff |r| ≤ 2|cos(φ/2)|**, **iff the
> zeros of L̃(χ) and L̃(χ̄) interlace on the circle.** (Elementary; proof on the page.)

The page recomputes the entire census in the browser — generator, discrete logarithms, character sums,
roots, interlacing test, mixture roots — for q = 3, 5, 7 (194 conjugate pairs: 3/3, 32/14, 99/43 on/off)
and checks every number against the committed Python twin. It ends with the fence: over ℚ the first half
(constituents on the line) is the generalized Riemann Hypothesis for χ mod 5, unproven; the second half
(the mixture off the line) is a theorem; and what the dial decides in one inequality is, over ℚ, the
Riemann Hypothesis itself. Nothing here bears on RH.

## Files

| file | what |
|---|---|
| `index.html` | the room — self-contained; one external request (Google Fonts); no analytics |
| `receipts/interlacing_dial.py` | the Python twin: recomputes the census from field arithmetic, asserts inequality ⟺ interlacing ⟺ numeric on-circle for every pair, writes `interlacing_dial.json` |
| `receipts/interlacing_dial.json` | the twin's output, embedded in the page as the self-check target |
| `receipts/index.template.html`, `receipts/build_page.py` | the page is `build_page.py` applied to the template with the JSON injected |
| `receipts/dh_dial_and_hn.py` (+ `.log`) | the ℚ side (2026-09-23): κ₊ = tan(½ arg ε) from the Gauss sum, the κ-dial, multiplicativity ⟺ κ² = −1, the mirror orbit at the certified off-line zero, the degree-2/3/4 function-field census (194 + 314 pairs), and the plane-colouring arithmetic of the same day |
| `receipts/dh_online_scan.py` (+ `_78_94`, `_114`, `_166_176.log`) | the on-line scans: zeros of L(χ), L(χ̄) and of the even/odd Davenport–Heilbronn functions around the four off-line zeros below t = 200 — the seam specimens quoted on the page |

Rebuild: `pip install mpmath && python3 receipts/interlacing_dial.py && python3 receipts/build_page.py`.

## What is and is not claimed

* [T] The degree-2 theorem above (two lines). The one-way half, interlacing ⟹ zeros on the circle, holds in every degree (Hermite–Biehler on the unit circle); the converse fails from degree 3 (census in `dh_dial_and_hn.py`: 0 interlaced pairs off the circle in 314, but 61 non-interlaced pairs on it at q = 5).
* [T, literature] Over ℚ: Davenport–Heilbronn 1936; Saias–Weingartner 2009 (a positive proportion of zeros in every strip inside ½ < σ ≤ 1 + η unless the series is P(s)·L(s,χ)); Kaczorowski–Kulas 2007 (degree one + functional equation + RH ⟹ Euler product); Bombieri–Hejhal 1995; Karatsuba 1990; Lagarias 1999.
* [Numerics] The four seam specimens (mpmath, 18 digits) — four specimens, no theorem; seams are not sufficient.
* Nothing bears on the Riemann Hypothesis. Instruments, not proofs.

Built 2026-09-24 by Travis Bergen with the Riemann agent. Companion rooms: [The Angle Atlas](https://fractalyouniverse.org/Angle_atlas/) (VII), [The Product's Edge](https://fractalyouniverse.org/Product_edge/) (VIII).
