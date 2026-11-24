"""Add the 'Fair comparison' tab after Learning loop: baseline vs full system (+ training-only ablation),
3 repeats each on the 65 held-out questions, and the end-to-end DFT check.
usage: add_fair_tab.py <page.html> <compare_v2.json>
"""
import json
import re
import statistics as st
import sys
from collections import Counter

path, cmp_path = sys.argv[1:3]
page = open(path).read()
assert 'id="panel-fair"' not in page
D = json.load(open(cmp_path))
ids = iter(range(4000, 4999))


def e():
    return f'data-e="e{next(ids)}"'


def td(c, num=True, cls=""):
    k = " ".join(x for x in ["num" if num else "", cls] if x)
    return f'<td {e()}{f" class={chr(34)}{k}{chr(34)}" if k else ""}>{c}</td>'


def th(h, num=True):
    return f'<th {e()}{" class=num" if num else ""}>{h}</th>'


A = D["arms"]
ARMS = [("A", "Baseline", "nemotron-3-nano (base)", "off", "off"),
        ("C", "Training only (ablation)", "nemotron-3-nano-learn-v2", "off", "off"),
        ("B", "Full system", "nemotron-3-nano-learn-v2", "on", "on")]
REFS = [("REF_base_critic", "Earlier: base model, critic on", "nemotron-3-nano (base)", "on", "off"),
        ("REF_ft1_critic", "Earlier: round-1 model, critic on", "nemotron-3-nano-learn-v1", "on", "off")]


def stat(a, k, fmt=lambda x: str(x), pct=False):
    s = a[k]
    v = s["values"]
    f = (lambda x: f"{100 * x:.0f}%") if pct else fmt
    if len(v) == 1:
        return f"{f(v[0])}"
    return f"<strong>{f(s['mean']) if pct else round(s['mean'], 1)}</strong><small>{' · '.join(f(x) for x in v)}</small>"


def mod(a, m):
    x, n = a["first_attempt_accept_by_module"].get(m, [0, 0])
    return f"{100 * x / n:.0f}%<small>{x} / {n}</small>" if n else "–"


# ---- main table -------------------------------------------------------------
heads = ["Arm", "Model", "Plan critic", "RAG", "Correct first try (of 65)", "After one hinted retry", "No plan (budget ran out)",
         "Wrong plan (caught by oracle)", "Replies that were not tool calls", "Attempts with a content mistake",
         "Target accepted 1st attempt", "Fragments accepted 1st attempt", "Attempts per campaign", "Planning time, 65 tasks (min)"]
thead = "".join(th(h, i >= 4) for i, h in enumerate(heads))


def row(key, name, model, critic, rag, ref=False):
    a = A[key]
    wrong = {"values": [x for x in a["wrong_plan_caught_by_oracle"]["values"]], "mean": a["wrong_plan_caught_by_oracle"]["mean"]}
    cells = [td(f"<strong>{key if not ref else '–'}</strong><small>{name}</small>", False), td(f"<code>{model}</code>", False),
             td(critic, False), td(rag, False),
             td(stat(a, "first_try_correct")), td(stat(a, "after_retry_correct")), td(stat(a, "no_plan")),
             td(stat({"w": wrong}, "w")), td(stat(a, "format_failure_rate", pct=True)), td(stat(a, "content_mistake_rate", pct=True)),
             td(mod(a, "target")), td(mod(a, "fragments")), td(stat(a, "attempts_per_campaign")), td(stat(a, "plan_wall_min"))]
    return f'<tr{" class=ref" if ref else ""}>' + "".join(cells) + "</tr>"


rows = "".join(row(*r) for r in ARMS)
refrows = "".join(row(*r, ref=True) for r in REFS)

# ---- paired test ------------------------------------------------------------
P = D["paired"]
pt = []
for key, lab in [("CA", "Training only (C) vs baseline (A)"), ("BA", "Full system (B) vs baseline (A)"), ("BC", "Full system (B) vs training only (C)")]:
    w, l, p = P[key]
    x, y = key[0], key[1]
    diff = A[x]["first_try_correct"]["mean"] - A[y]["first_try_correct"]["mean"]
    pt.append("<tr>" + td(lab, False) + td(f"{diff:+.1f}") + td(str(w)) + td(str(l)) + td(str(65 - w - l)) + td(f"{p:.2f}") + "</tr>")
pthead = "".join(th(h, i > 0) for i, h in enumerate(["Comparison", "Difference in mean correct", "Tasks solved more often",
                                                     "Tasks solved less often", "Tasks tied", "Sign test p (two-sided)"]))
cnt = D["per_task"]
cons = "".join("<tr>" + td(f"<strong>{k}</strong> {n}", False) + td(str(sum(v == 3 for v in cnt[k].values())))
               + td(str(sum(0 < v < 3 for v in cnt[k].values()))) + td(str(sum(v == 0 for v in cnt[k].values()))) + "</tr>"
               for k, n, *_ in ARMS)
conshead = "".join(th(h, i > 0) for i, h in enumerate(["Arm", "Correct in all 3 repeats", "Correct in 1 or 2", "Never correct"]))

# ---- failure classes and bond classes ---------------------------------------
fc = "".join("<tr>" + td(f"<strong>{k}</strong> {n}", False) + "".join(td(str(A[k]["failure_classes_total"].get(c, 0)))
             for c in ["no_plan", "wrong_bond", "wrong_parent"]) + td(str(sum(A[k]["failure_classes_total"].values()))) + "</tr>"
             for k, n, *_ in ARMS)
fchead = "".join(th(h, i > 0) for i, h in enumerate(["Arm (3 × 65 = 195 campaigns)", "No plan in budget", "Wrong bond", "Wrong parent", "All misses"]))
classes = sorted(A["A"]["by_bond_class"], key=lambda c: -A["A"]["by_bond_class"][c][1])
bc = "".join("<tr>" + td(c, False) + td(str(A["A"]["by_bond_class"][c][1] // 3)) + "".join(
    td(f"{A[k]['by_bond_class'][c][0]} / {A[k]['by_bond_class'][c][1]}") for k in ["A", "C", "B"]) + "</tr>" for c in classes)
bchead = "".join(th(h, i > 0) for i, h in enumerate(["Bond class", "Tasks", "A baseline", "C training only", "B full system"]))

# ---- end to end -------------------------------------------------------------
EA, EB = D["e2e"]["A"], D["e2e"]["B"]
rb = {r["task"]: r for r in EB["rows"]}


def cell(r):
    if r["plan_correct"]:
        return f'<span class="ok">correct</span> · {r["bde_kcal"]:.1f} · {r["results_verdict"]}'
    if r["bde_kcal"] is not None:
        return f'<span class="bad">{r["error_class"].replace("_", " ")}</span> · {r["bde_kcal"]:.1f} · {r["results_verdict"]}'
    return f'<span class="bad">{(r["error_class"] or "?").replace("_", " ")}</span> · not computed'


e2e = "".join("<tr>" + td(r["task"].split(":")[0] + f'<small>{r["task"].split(":")[1].replace("_", " ")}</small>', False)
              + td(cell(r), False) + td(cell(rb[r["task"]]), False) + "</tr>" for r in EA["rows"])
e2ehead = "".join(th(h, False) for h in ["Question (molecule, bond)", "A baseline: plan · BDE kcal/mol · results check",
                                         "B full system: plan · BDE kcal/mol · results check"])


def e2esum(E):
    wrong_pass = sum(1 for r in E["rows"] if not r["plan_correct"] and r["bde_kcal"] is not None and r["results_verdict"] == "PASS")
    return [E["plan_correct"], E["reports"], E["completed_correct"], wrong_pass, f'{E["species_ok"]} / {E["species_tasks"]}',
            E["db_hits"], E["wall_min"]["plan"], round(E["wall_min"]["compute"] + E["wall_min"]["assemble"], 1)]


e2es = "".join("<tr>" + td(f"<strong>{k}</strong> {n}", False) + "".join(td(str(x)) for x in e2esum(E)) + "</tr>"
               for k, n, E in [("A", "Baseline", EA), ("B", "Full system", EB)])
e2eshead = "".join(th(h, i > 0) for i, h in enumerate(
    ["Arm (17 questions)", "Correct plans", "Reports written", "Correct plan + computed BDE", "Wrong plan computed and passed the results check",
     "Species computed OK", "Species reused from DB", "Planning (min)", "DFT + assembly (min)"]))

a, b, c = (A[k]["first_try_correct"]["mean"] for k in "ABC")
SECTION = f"""<section class="rplan">
    <h2 id="fair-comparison" {e()}>Fair comparison: baseline vs the full system, three repeats each</h2>
    <p class="lede" {e()}>Does the whole learning loop make the planner better? nemotron-3-nano with no plan critic and no retrieval (A) against the round-2 fine-tuned model with the plan critic and retrieval from the database of computed bonds (B), on the same 65 held-out questions, three full repeats each, plus an end-to-end DFT check. A training-only arm (C) separates the effect of the fine-tuning from the critic and retrieval. Run on Frontier 2026-11-23/24.</p>
    <div class="tiles">
      <div class="tile"><div class="big" {e()}>{a:.1f} · {c:.1f} · {b:.1f}</div><div class="cap" {e()}>mean correct first try (of 65) for A baseline · C training only · B full system</div></div>
      <div class="tile"><div class="big" {e()}>±2–3</div><div class="cap" {e()}>run-to-run spread (standard deviation over 3 repeats) of the same arm: larger than any difference between arms</div></div>
      <div class="tile"><div class="big" {e()}>p ≥ 0.58</div><div class="cap" {e()}>paired sign test per question, every pair of arms: no significant difference</div></div>
      <div class="tile"><div class="big" {e()}>12 vs 11</div><div class="cap" {e()}>correct plans that went on to a computed BDE in the end-to-end DFT check (B vs A, 17 questions)</div></div>
    </div>
    <div class="findings">
      <div class="finding bad"><h3 {e()}>No measurable gain from the full system</h3><p {e()}>B solved {b:.1f} of 65 on average against {a:.1f} for A. The difference is within the spread of repeats of one arm, and per question B did better on {P['BA'][0]} and worse on {P['BA'][1]}.</p></div>
      <div class="finding bad"><h3 {e()}>Round-2 training alone did not help either</h3><p {e()}>C (fine-tuned v2, critic and retrieval off) averaged {c:.1f}, the same as A. Training did change <em>how</em> the planner fails: fewer replies that are not tool calls (52% → 42% of attempts) but more content mistakes (10% → 19%), so more plans reach the oracle and more of them are wrong.</p></div>
      <div class="finding"><h3 {e()}>The earlier 40 → 50 gain was mostly the plan critic</h3><p {e()}>The Learning loop tab compared single runs with the critic on: base 40, round-1 model 50. With the critic off and three repeats, the base model alone reaches {a:.1f}. The critic was holding the base model back (24 budget timeouts in that run against 12–16 here); fine-tuning made the model more tolerant of the critic rather than better at planning.</p></div>
      <div class="finding"><h3 {e()}>What the critic and retrieval change</h3><p {e()}>With the same fine-tuned model, turning them on (C → B) halves first-attempt acceptance of the target module (57% → 31%) and gives more no-plan campaigns ({A['C']['failure_classes_total']['no_plan']} → {A['B']['failure_classes_total']['no_plan']} of 195), but fewer wrong plans ({sum(v for k,v in A['C']['failure_classes_total'].items() if k!='no_plan')} → {sum(v for k,v in A['B']['failure_classes_total'].items() if k!='no_plan')}). The two effects roughly cancel. Planning time is set by the model, not by them: about 220 minutes per 65 questions for both fine-tuned arms against 164 for the base model.</p></div>
      <div class="finding bad"><h3 {e()}>The results check cannot see a wrong plan</h3><p {e()}>In the end-to-end run, 3 wrong plans were computed to completion and all 3 passed the post-calculation QC check (A: styrene, wrong ring position; B: pyridine wrong position, vinyl chloride wrong parent, 64 kcal/mol). Only a check against the question, like the oracle, catches them.</p></div>
      <div class="finding good"><h3 {e()}>The computations themselves are reproducible</h3><p {e()}>Where both arms chose the same species and method, the BDE agreed to 0.1 kcal/mol (hydrazine 72.2, benzene 117.6, butadiene 120.0, acetamide 116.6, bromoethane 95.2). Every species task that ran completed; no SCF failed to converge.</p></div>
    </div>
    <div class="summary">
      <h3 {e()}>What was held fixed, and what changed</h3>
      <div class="tw"><table style="min-width:760px"><thead><tr>{th('Held fixed for every arm', False)}{th('Value', False)}</tr></thead><tbody>
        <tr>{td('Questions', False)}{td('65 held-out tasks (learn/tasks.json, split eval), the same three parts: 0–21, 22–43, 44–64', False)}</tr>
        <tr>{td('Repeats', False)}{td('3 full repeats per arm (9 runs), arms interleaved on the same queues: one repeat on service/develop, two on batch/normal', False)}</tr>
        <tr>{td('Serving', False)}{td('1 Frontier node per run, 8 Ollama servers (one per GCD), num_parallel 1, context 8192, same quantised base weights (LoRA loaded as an ADAPTER on nemotron-3-nano:latest)', False)}</tr>
        <tr>{td('Planning protocol', False)}{td('plan board v2, 30-minute budget per campaign, Flowcept tracing on, no lessons', False)}</tr>
        <tr>{td('Grading', False)}{td('learn/oracle.py against the RDKit-computed answer; a wrong plan sent back once with a hint', False)}</tr>
        <tr>{td('Database', False)}{td('each run gets its own fresh copy of the frozen db/rag_eval_v1.sqlite (103 species, 90 BDEs, sha256 c339a259…); with RAG off it is not read for retrieval', False)}</tr>
        <tr>{td('End-to-end check', False)}{td('17 of the 65 (every fourth), plan → xtb + DFT → assembly + results check, first-try plan only (no oracle hint, as in real use)', False)}</tr>
      </tbody></table></div>
      <div class="tw"><table style="min-width:760px;margin-top:12px"><thead><tr>{th('Changed', False)}{th('A baseline', False)}{th('C training only', False)}{th('B full system', False)}</tr></thead><tbody>
        <tr>{td('Model', False)}{td('nemotron-3-nano', False)}{td('+ LoRA round 2', False)}{td('+ LoRA round 2', False)}</tr>
        <tr>{td('Plan critic (LLM)', False)}{td('off: validator only', False)}{td('off', False)}{td('on', False)}</tr>
        <tr>{td('Retrieval of computed bonds (RAG)', False)}{td('off', False)}{td('off', False)}{td('on', False)}</tr>
      </tbody></table></div>
    </div>
    <div class="summary">
      <h3 {e()}>Results on the 65 held-out questions (planning only)</h3>
      <div class="tw"><table style="min-width:1500px"><thead><tr>{thead}</tr></thead><tbody>
        {rows}
        <tr><td colspan="14" {e()} style="font-size:12px;padding-top:14px"><strong>Reference, not part of the controlled comparison:</strong> the single runs from the Learning loop tab (critic on, RAG off, one repeat)</td></tr>
        {refrows}
      </tbody></table></div>
      <p {e()} style="font-size:13px">Bold = mean of 3 repeats; small = each repeat. Module acceptance and failure counts are pooled over the 3 repeats. Planning time is the sum of the three parts' planning stages. "Not tool calls" and "content mistake" are fractions of all planner attempts.</p>
    </div>
    <div class="summary">
      <h3 {e()}>Is any difference real? Paired per question</h3>
      <p {e()} style="font-size:13px">For each question, the number of repeats (0–3) in which the arm's first plan was correct, compared between two arms.</p>
      <div class="tw"><table style="min-width:760px"><thead><tr>{pthead}</tr></thead><tbody>{"".join(pt)}</tbody></table></div>
      <div class="tw"><table style="min-width:560px;margin-top:12px"><thead><tr>{conshead}</tr></thead><tbody>{cons}</tbody></table></div>
    </div>
    <div class="summary">
      <h3 {e()}>How the misses break down</h3>
      <div class="tw"><table style="min-width:640px"><thead><tr>{fchead}</tr></thead><tbody>{fc}</tbody></table></div>
      <div class="tw"><table style="min-width:640px;margin-top:12px"><thead><tr>{bchead}</tr></thead><tbody>{bc}</tbody></table></div>
      <p {e()} style="font-size:13px">Bond-class cells: correct first tries / campaigns, pooled over 3 repeats. C–Cl was held out of all training data.</p>
    </div>
    <div class="summary">
      <h3 {e()}>End-to-end DFT check (17 held-out questions)</h3>
      <div class="tw"><table style="min-width:1000px"><thead><tr>{e2eshead}</tr></thead><tbody>{e2es}</tbody></table></div>
      <div class="tw"><table style="min-width:900px;margin-top:12px"><thead><tr>{e2ehead}</tr></thead><tbody>{e2e}</tbody></table></div>
      <p {e()} style="font-size:13px">BDE = electronic energy difference, xtb geometries with a B3LYP/6-31G* single point (the method both arms' planners chose), no zero-point correction, so values sit above experiment. The results check is the post-calculation QC verdict in each report. One run per arm: 12 against 11 correct is not a meaningful difference.</p>
    </div>
    <div class="summary">
      <h3 {e()}>Caveats</h3>
      <ul class="next" style="list-style:disc">
        <li {e()}>A against B changes three things at once (model, critic, retrieval); C isolates the model. The effects of the critic and of retrieval are not separated from each other.</li>
        <li {e()}>Training data came from campaigns run with the critic on and retrieval off, so B runs the fine-tuned model in a setting it was not trained in, and C in one closer to the baseline's.</li>
        <li {e()}>The retrieval database holds earlier results for 4 of the eval molecules; with RAG on, B can see them. It did not translate into more correct plans.</li>
        <li {e()}>Three repeats detect a difference of roughly 5 questions in 65; smaller effects would need more repeats.</li>
        <li {e()}>Round 2 added one collection run to round 1's data (519 SFT examples and 346 DPO pairs against 455 and 305). It was not re-evaluated with the critic on, so round 2 against round 1 is not compared here.</li>
      </ul>
    </div>
    <div class="summary">
      <h3 {e()}>Training round 2</h3>
      <p {e()}>Episodes from the 8 round-1 collection runs plus the rerun of the failed one (2,011 attempts) → 519 SFT examples, 346 DPO pairs. SFT and DPO ran as two jobs of under 2 h on batch/normal (the second starting from the first's adapter, <code>train_lora --init-adapter</code>): SFT validation loss 0.583 (round 1: 0.587); DPO reward accuracy on held-out pairs 0.935 (round 1: 0.964). Converted with llama.cpp and registered as <code>nemotron-3-nano-learn-v2</code> in a third compute-node job.</p>
    </div>
    <div class="summary">
      <h3 {e()}>Files</h3>
      <p {e()} style="font-size:13px">Scripts: <code>frontier-native/experiments/learn/</code> (run_arm.sh, arm_chain.sh, train_v2_sft.sbatch, train_v2_dpo.sbatch, register_v2.sbatch); the critic switch <code>--no-plan-critic</code> / <code>COMPCHEM_PLAN_CRITIC=off</code>; summary <code>src/learn/compare_arms.py</code>. Results: <code>/lustre/orion/stf053/scratch/akhilesh/compchem/learn/eval/compare_v2.json</code>; runs in <code>learn/runs/</code> (driver_armA_*, driver_armB_*, driver_armC_* logs name each run directory); end-to-end runs <code>{EA['run']}</code> (A) and <code>{EB['run']}</code> (B); adapters in <code>learn/train/adapter_v2_sft</code>, <code>adapter_v2</code>.</p>
    </div>
  </section>"""

bar = '<button type="button" role="tab" id="tab-log"'
assert page.count(bar) == 1
page = page.replace(bar, '<button type="button" role="tab" id="tab-fair" aria-controls="panel-fair" aria-selected="false">Fair comparison</button>' + bar, 1)
anchor = '<div class="tabpanel" role="tabpanel" id="panel-log"'
assert page.count(anchor) == 1
page = page.replace(anchor, f'<div class="tabpanel" role="tabpanel" id="panel-fair" aria-labelledby="tab-fair" hidden>\n  {SECTION}\n</div>\n\n' + anchor, 1)
m = re.search(r'<li data-e="e\d+"><a href="#learning-loop">Learning loop</a>[^<]*</li>', page)
assert m
page = page.replace(m.group(0), m.group(0) + f'\n          <li {e()}><a href="#fair-comparison">Fair comparison</a> Baseline against the full system (fine-tuned + critic + retrieval), three repeats each, with a training-only ablation and an end-to-end DFT check.</li>', 1)
css = ".finding.bad { border-left-color: var(--fail); }"
assert css in page
page = page.replace(css, css + "\ntr.ref td { opacity: .75; }\ntd .ok { color: var(--pass); font-weight: 600; }\ntd .bad { color: var(--fail); font-weight: 600; }", 1)
dup = [x for x, n in Counter(re.findall(r'data-e="(e\d+)"', page)).items() if n > 1]
assert not dup, dup
open(path, "w").write(page)
print("ok", len(page), "ids to", next(ids) - 1)
