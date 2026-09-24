#!/usr/bin/env python3
"""
build_page.py — inject the Python twin's census into the page template.

    python3 receipts/interlacing_dial.py      # writes receipts/interlacing_dial.json
    python3 receipts/build_page.py            # writes index.html from receipts/index.template.html

The page recomputes the whole census in the browser; the injected object is only the
self-check target (per pair: [j, r, phi_deg, on]).  Nothing else is copied in.
"""
import json, math, os, hashlib

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
data = json.load(open(os.path.join(here, "interlacing_dial.json")))

compact = {
    "generated": data["generated"],
    "kappa_plus": 0.28407904384041229603,          # dh_dial_and_hn.py (2026-09-23): tan(arg eps / 2) = D-H radical to 2e-31
    "arg_eps_deg": 31.7174744115,                  # same receipt
    "census": [],
}
for c in data["census"]:
    compact["census"].append({
        "q": c["q"], "m": c["m_low_to_high"], "g": c["generator"], "N": c["N"],
        "odd": c["odd_characters"], "pairs": c["pairs"], "on": c["on"], "off": c["off"],
        "rows": [[r["j"], round(r["r"], 9), round(r["phi_deg"], 6), 1 if r["on"] else 0] for r in c["rows"]],
    })
compact["totals"] = data["totals"]
blob = json.dumps(compact, separators=(",", ":"))

tpl = open(os.path.join(here, "index.template.html"), encoding="utf-8").read()
assert tpl.count("__PY_CENSUS__") == 1
html = tpl.replace("__PY_CENSUS__", blob)
# Greek letter + combining macron (chi-bar) is not rendered by the page fonts in every browser (tofu boxes were
# observed in headless Chromium).  In the HTML text, draw the bar with CSS instead; the canvas code draws its own bars.
head, script = html.split("<script>", 1)
head = head.replace("\u03c7\u0304", '<span class="ol">\u03c7</span>')
assert "\u03c7\u0304" not in head
html = head + "<script>" + script
out = os.path.join(root, "index.html")
open(out, "w", encoding="utf-8").write(html)
print(f"wrote {out}: {len(html):,} bytes, census blob {len(blob):,} bytes, sha256 {hashlib.sha256(html.encode()).hexdigest()}")
print("totals:", compact["totals"])
