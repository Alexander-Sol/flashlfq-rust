# -*- coding: utf-8 -*-
import base64, os
D = os.environ["CLAUDE_JOB_DIR"] + r"\tmp"

def b64(p):
    with open(p, "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()

fig_detect  = b64(D + r"\decoy_hist.png")
fig_refined = b64(D + r"\decoy_hist_refined.png")
fig_pos     = b64(D + r"\decoy_positional.png")

HTML = f"""<title>Decoy comb FDR — first experiment</title>
<style>
:root {{
  --bg:#f6f8fb; --surface:#ffffff; --ink:#161a21; --muted:#59616e;
  --hair:#e2e6ee; --target:#2f6fd0; --decoy:#cf3f56; --warn:#b8791a;
  --mono:ui-monospace,"SF Mono","Cascadia Mono","JetBrains Mono",Menlo,Consolas,monospace;
  --sans:system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
}}
@media (prefers-color-scheme:dark) {{
  :root {{ --bg:#0e1116; --surface:#161b22; --ink:#e7ebf2; --muted:#98a2b2;
    --hair:#242b35; --target:#5b95e6; --decoy:#e8697e; --warn:#d69a3f; }}
}}
:root[data-theme="dark"] {{ --bg:#0e1116; --surface:#161b22; --ink:#e7ebf2; --muted:#98a2b2;
  --hair:#242b35; --target:#5b95e6; --decoy:#e8697e; --warn:#d69a3f; }}
:root[data-theme="light"] {{ --bg:#f6f8fb; --surface:#ffffff; --ink:#161a21; --muted:#59616e;
  --hair:#e2e6ee; --target:#2f6fd0; --decoy:#cf3f56; --warn:#b8791a; }}

* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:var(--sans);
  line-height:1.6; -webkit-font-smoothing:antialiased; }}
.wrap {{ max-width:820px; margin:0 auto; padding:64px 24px 96px; }}
.eyebrow {{ font-family:var(--mono); font-size:12px; letter-spacing:.16em; text-transform:uppercase;
  color:var(--muted); margin:0 0 14px; }}
h1 {{ font-size:34px; line-height:1.15; letter-spacing:-.02em; margin:0 0 12px; text-wrap:balance; font-weight:650; }}
.lede {{ font-size:18px; color:var(--muted); margin:0; max-width:64ch; text-wrap:pretty; }}
h2 {{ font-size:14px; font-family:var(--mono); letter-spacing:.12em; text-transform:uppercase;
  color:var(--muted); margin:52px 0 16px; padding-bottom:8px; border-bottom:1px solid var(--hair); }}
p {{ max-width:66ch; }}
.legend {{ display:flex; gap:20px; margin:18px 0 0; font-family:var(--mono); font-size:13px; }}
.legend span {{ display:inline-flex; align-items:center; gap:7px; }}
.dot {{ width:11px; height:11px; border-radius:3px; display:inline-block; }}
.dot.t {{ background:var(--target); }} .dot.d {{ background:var(--decoy); }}

.verdict {{ margin:28px 0 8px; background:var(--surface); border:1px solid var(--hair);
  border-left:3px solid var(--decoy); border-radius:10px; padding:20px 22px; }}
.verdict .k {{ font-family:var(--mono); font-size:12px; letter-spacing:.14em; text-transform:uppercase;
  color:var(--decoy); margin:0 0 8px; }}
.verdict p {{ margin:0; font-size:16px; }}

.caveat {{ margin:0 0 30px; background:color-mix(in srgb, var(--warn) 8%, var(--surface));
  border:1px solid color-mix(in srgb, var(--warn) 35%, var(--hair)); border-left:3px solid var(--warn);
  border-radius:10px; padding:16px 20px; }}
.caveat .k {{ font-family:var(--mono); font-size:12px; letter-spacing:.14em; text-transform:uppercase;
  color:var(--warn); margin:0 0 8px; }}
.caveat p {{ margin:0; font-size:14.5px; max-width:none; }}

figure {{ margin:22px 0 8px; }}
.plate {{ background:#ffffff; border:1px solid var(--hair); border-radius:10px; padding:10px;
  overflow-x:auto; box-shadow:0 1px 2px rgba(20,26,40,.04); }}
.plate img {{ display:block; width:100%; height:auto; }}
figcaption {{ font-family:var(--mono); font-size:12.5px; color:var(--muted); margin-top:10px; max-width:70ch; }}

.tbl {{ overflow-x:auto; margin:16px 0; }}
table {{ border-collapse:collapse; width:100%; font-family:var(--mono); font-size:13px;
  font-variant-numeric:tabular-nums; }}
th,td {{ text-align:right; padding:8px 12px; border-bottom:1px solid var(--hair); white-space:nowrap; }}
th:first-child,td:first-child {{ text-align:left; }}
thead th {{ color:var(--muted); font-weight:600; letter-spacing:.04em; border-bottom:1.5px solid var(--hair); }}
tbody tr:last-child td {{ border-bottom:none; }}
.t-col {{ color:var(--target); }} .d-col {{ color:var(--decoy); }}
.hi {{ background:color-mix(in srgb, var(--decoy) 9%, transparent); }}

ul {{ max-width:66ch; padding-left:20px; }} li {{ margin:6px 0; }}
code {{ font-family:var(--mono); font-size:.9em; background:color-mix(in srgb,var(--ink) 7%,transparent);
  padding:1px 5px; border-radius:4px; }}
.foot {{ margin-top:56px; padding-top:18px; border-top:1px solid var(--hair);
  font-family:var(--mono); font-size:12px; color:var(--muted); line-height:1.9; }}
.strong {{ color:var(--ink); font-weight:600; }}
</style>

<div class="wrap">
  <div class="caveat">
    <p class="k">Caveat · flawed pathway</p>
    <p>Every result below comes from the <span class="strong">detector-decoy / refiner-real</span>
    pathway: the decoy model was applied only at feature <em>detection</em>, while <em>refinement</em>
    re-anchored the mono/charge and scored against the <em>real</em> averagine. That pathway
    <span class="strong">launders the decoy</span> — it re-fits the decoy hypothesis back onto real
    signal before scoring — so the target-vs-decoy separations here reflect that flaw and understate
    what a correct <em>end-to-end</em> decoy (the decoy model carried through refinement too) would
    show. Read these as diagnostics of the flawed pathway, not as the achievable FDR.</p>
  </div>
  <p class="eyebrow">MS1 feature detection · target–decoy FDR</p>
  <h1>Can a chlorinated-averagine decoy give us known negatives?</h1>
  <p class="lede">To build a target–decoy FDR we need known negatives: features the detector accepts
  that cannot be real analytes. We tried three ways to manufacture them inside the identical detection
  tree — a decoy isotope comb whose <em>shape</em> is wrong, one whose <em>spacing</em> is wrong, and
  one whose spacing is <em>globally impossible</em> — and asked whether decoy detections separate from
  real (averagine) targets.</p>

  <div class="verdict">
    <p class="k">Verdict — shape fails, position separates weakly</p>
    <p>A wrong-<em>shape</em> comb (on-lattice) does <strong>not</strong> separate at all — it fits as
    well or better than targets. Wrong-<em>position</em> combs (off-lattice, mixed-charge)
    <strong>do</strong> separate, but only on isotope-ladder length, not on fit score, and not cleanly.
    The usable signal is envelope completeness, not "is there an isotope-spaced peak here."</p>
  </div>

  <p class="eyebrow" style="margin:40px 0 0">Experiment 1 · wrong shape</p>

  <h2>What was compared</h2>
  <p>One detection pass with the averagine comb (<span class="strong">target</span>, 96,385 features)
  and one with the chlorinated-averagine comb (<span class="strong">decoy</span>, 60,996 features) on
  the CA/Lumos 10-min run at 90% coverage. Same seeds, same tolerances, same everything but the comb
  weights. The decoy's ³⁷Cl A+2 ladder pushes its envelope mode 4–6 ¹³C units off the monoisotope —
  a shape no tryptic peptide produces.</p>
  <div class="legend">
    <span><i class="dot t"></i>target — averagine comb</span>
    <span><i class="dot d"></i>decoy — chlorinated-averagine comb</span>
  </div>

  <h2>Detection score</h2>
  <figure>
    <div class="plate"><img src="{fig_detect}" alt="Target vs decoy detection score histograms"></div>
    <figcaption>Left: raw matched-filter score (Σ w·g·I) — decoy runs <em>higher</em> than target
    because the score scales with intensity and the wider comb claims more peaks. Right: score ÷
    summed intensity, a crude shape/RT-fit proxy — the two distributions sit on top of each other.</figcaption>
  </figure>
  <div class="tbl">
    <table>
      <thead><tr><th>fit proxy ≥</th><th class="t-col">target</th><th class="d-col">decoy</th><th>decoy / target</th></tr></thead>
      <tbody>
        <tr><td>0.30</td><td>85,458</td><td>57,185</td><td>0.669</td></tr>
        <tr><td>0.50</td><td>64,650</td><td>46,108</td><td>0.713</td></tr>
        <tr><td>0.70</td><td>26,692</td><td>20,937</td><td>0.784</td></tr>
        <tr><td>0.80</td><td>13,779</td><td>9,784</td><td>0.710</td></tr>
      </tbody>
    </table>
  </div>
  <p>A valid decoy would thin out as fit quality rises — the ratio should fall toward zero. It stays
  flat near 0.7–0.8. No separation.</p>

  <h2>Refined envelope-fit (the strict metric)</h2>
  <p>The detection score is lenient — it never penalizes unexplained peaks or missing teeth. The
  refined <code>Decon Score</code> (envelope-fit cosine over the union grid) does. This is where
  separation should appear if anywhere.</p>
  <figure>
    <div class="plate"><img src="{fig_refined}" alt="Target vs decoy envelope-fit score histogram"></div>
    <figcaption>Envelope-fit cosine after shift-decon refinement. The decoy distribution is shifted
    slightly <em>right</em> of the target — it fits marginally better, not worse.</figcaption>
  </figure>
  <div class="tbl">
    <table>
      <thead><tr><th>Decon Score ≥</th><th class="t-col">target</th><th class="d-col">decoy</th><th>decoy / target</th></tr></thead>
      <tbody>
        <tr><td>0.70</td><td>25,112</td><td>18,935</td><td>0.754</td></tr>
        <tr><td>0.80</td><td>20,369</td><td>16,316</td><td>0.801</td></tr>
        <tr><td>0.90</td><td>14,480</td><td>12,506</td><td>0.864</td></tr>
        <tr class="hi"><td>0.99</td><td>3,821</td><td>3,462</td><td>0.906</td></tr>
      </tbody>
    </table>
  </div>
  <p>The ratio doesn't collapse at high confidence — it <strong>climbs</strong>, to 0.91 at Decon
  Score ≥ 0.99. The most confident decoy fits are nearly as numerous as the most confident targets.
  Precisely backwards from a working null.</p>

  <h2>Why it fails</h2>
  <p>The decoy only changed the <em>detection</em> comb — but the data underneath is still real
  peptide signal, and two things let it through:</p>
  <ul>
    <li><span class="strong">The teeth stay on the ¹³C lattice.</span> Composition barely moves the
    spacing (~3–5 mDa/step, inside 10 ppm tolerance — measured earlier). So the decoy comb still lands
    on real, abundant isotope-structured peaks.</li>
    <li><span class="strong">Refinement re-optimizes freely.</span> The envelope-fit is scored against
    the real averagine template after shift-decon re-anchors the window — so both target and decoy
    seeds converge onto whatever real envelope best explains the data. The decoy label never
    propagates to the score.</li>
  </ul>
  <p>Root cause in one line: <span class="strong">a decoy that only reshapes the comb still queries
  real signal.</span> To manufacture negatives you have to query where real signal <em>cannot</em>
  be — off the isotope lattice.</p>

  <p class="eyebrow" style="margin:52px 0 0">Experiment 2 · wrong position</p>
  <h2>Off-lattice &amp; mixed-charge lattices</h2>
  <p>So we moved the teeth off the real isotope lattice, keeping the normal averagine weights.
  <span class="strong">Off-lattice ×0.5</span> and <span class="strong">×0.9</span> scale the ¹³C
  spacing by a non-physical factor (querying the gaps between real peaks); <span class="strong">
  mixed-charge</span> gives every inter-tooth step a different <code>1/z′</code> spacing, so each step
  is plausible but the sequence matches no real envelope. Judged at detection (positions are fixed
  there; refinement would launder them back).</p>
  <p>The fit proxy is still useless — it runs <em>higher</em> for decoys, because it only measures how
  well the peaks a comb <em>did</em> land on fit, and a 2-tooth coincidence on two real peaks looks
  perfect. The real signal is <span class="strong">how many teeth the comb can chain</span>: a real
  feature builds a long isotope ladder; a wrong-spacing comb mostly finds the bare 2-tooth minimum.</p>
  <figure>
    <div class="plate"><img src="{fig_pos}" alt="Positional decoy detection distributions"></div>
    <figcaption>Left: detection fit proxy — decoys sit at or above target (non-discriminating).
    Right: teeth found per feature — targets carry a visibly heavier tail of long ladders.</figcaption>
  </figure>
  <div class="tbl">
    <table>
      <thead><tr><th>per-feature rate</th><th class="t-col">target</th><th>off ×0.5</th><th>off ×0.9</th><th>mixed</th></tr></thead>
      <tbody>
        <tr><td>exactly 2 teeth</td><td>43.0%</td><td>57.8%</td><td class="hi">71.4%</td><td>63.5%</td></tr>
        <tr><td>≥ 4 teeth</td><td>31.9%</td><td>19.2%</td><td class="hi">11.4%</td><td>14.4%</td></tr>
        <tr><td>reached 90% coverage</td><td>yes</td><td>yes (89.7)</td><td class="hi">no — 83.3%</td><td>yes</td></tr>
      </tbody>
    </table>
  </div>
  <p>Rates are per-feature, so the coverage confound (decoys need ~1.6× more features) cancels. Targets
  are <span class="strong">1.7–2.8× richer in long ladders</span>. The ordering is itself telling:
  <span class="strong">×0.9 separates best</span>, not ×0.5 — an almost-right spacing lets the first
  one or two teeth match within tolerance but drifts teeth 3+ out of tolerance (cumulative drift caps
  the ladder), while a badly-wrong ×0.5 drops teeth into gaps that dense co-eluting data still fills.
  ×0.9 also can't explain 90% of the signal at all. Mixed-charge lands between the two.</p>
  <p>But none is a <em>clean</em> null: every positional decoy still yields 150k+ accepted features,
  many with 3–4 coincidental teeth, because MS1 data is dense enough to assemble short ladders at
  almost any spacing. The 2-tooth acceptance floor does most of the "detecting."</p>

  <h2>Where this leaves us</h2>
  <ul>
    <li><span class="strong">Shape decoys are dead</span> — real signal is on-lattice and refinement
    re-optimizes to it.</li>
    <li><span class="strong">Position decoys are a soft null</span> — real, but only under a score that
    rewards <em>envelope completeness</em> (ladder length, gap-penalized fit), not the lenient
    detection fit proxy. Detection's 2-tooth floor is too permissive to be the FDR score.</li>
  </ul>
  <p>Two moves make the soft null usable: (1) score with a <span class="strong">gap-penalizing
  envelope-fit</span> — one that charges for predicted-but-absent teeth (the union-grid cosine already
  in the refiner, applied at fixed decoy positions rather than re-anchored); and (2) a
  <span class="strong">paired, same-seed comparison</span> — score every seed under target and decoy
  combs without the claim/coverage loop, so counts share a denominator and FDR = decoy/target is exact.
  That combination is the honest test of whether ×0.9 or mixed-charge can carry an FDR.</p>

  <div class="foot">
    branch <span class="strong">worktree-decoy-comb-fdr</span> · commits 801105b, bdcbe03 · 211 lib tests pass<br>
    CombWeightModel::Decoy · LatticeMode {{Scaled, MixedCharge}} · COMB_MODEL / DECOY_LATTICE / DECOY_SPACING_SCALE<br>
    CA/Lumos 04-17-23_CA_Tryp_HCD_10min · cov90 · detect-only · data-driven σ 1.89 s
  </div>
</div>
"""

out = D + r"\decoy_report.html"
with open(out, "w", encoding="utf-8") as f:
    f.write(HTML)
print("wrote", out, "-", os.path.getsize(out), "bytes")
