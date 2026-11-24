"""Add 'How training changed the way the model fails' to the Fair comparison tab:
arm A (base) vs arm C (fine-tuned v2), both critic off / RAG off, from the first-try episodes.
usage: add_failure_breakdown.py <page.html> <episodes_armA.jsonl> <episodes_armC.jsonl>
"""
import json
import re
import sys
from collections import Counter, defaultdict

path, ea, ec = sys.argv[1:4]
page = open(path).read()
assert 'id="failure-breakdown"' not in page
ids = iter(range(4400, 4799))


def e():
    return f'data-e="e{next(ids)}"'


def td(c, num=True, cls=""):
    k = " ".join(x for x in ["num" if num else "", cls] if x)
    return f'<td {e()}{f" class={chr(34)}{k}{chr(34)}" if k else ""}>{c}</td>'


def th(h, num=True):
    return f'<th {e()}{" class=num" if num else ""}>{h}</th>'


E = {"A": [json.loads(l) for l in open(ea)], "C": [json.loads(l) for l in open(ec)]}
OK = ("accepted", "oracle_wrong", "critic_revise")


# 1. outcome of each attempt per module
def mod_rates(es, m):
    x = [e for e in es if e["module"] == m]
    n = len(x)
    f = sum(e["blame"] == "planner_format" for e in x)
    c = sum(e["blame"] == "planner_content" for e in x)
    a = sum(e["blame"] in OK for e in x)
    return n, f, c, a


def pc(k, n):
    return f"{100 * k / n:.0f}%<small>{k} / {n}</small>"


r1 = []
for m in ["target", "fragments", "method", "rationale"]:
    nA, fA, cA, aA = mod_rates(E["A"], m)
    nC, fC, cC, aC = mod_rates(E["C"], m)
    r1.append("<tr>" + td(f"<strong>{m}</strong>", False) + td(pc(fA, nA)) + td(pc(fC, nC), cls="hl" if fC < fA else "")
              + td(pc(cA, nA)) + td(pc(cC, nC)) + td(pc(aA, nA)) + td(pc(aC, nC)) + "</tr>")
h1 = "".join(th(h, i > 0) for i, h in enumerate(["Module", "Not a tool call: A", "C", "Content rejected: A", "C", "Accepted: A", "C"]))


# 2. kinds of format failure
def fkind(e):
    s = str(e["reply"].get("error", ""))
    t = s.split("plain text", 1)[-1]
    if "{" in t and ("smiles" in t.lower() or "parent" in t or '"name"' in t or "arguments" in t):
        return "json"
    return "prose"


fk = {a: Counter(fkind(e) for e in E[a] if e["blame"] == "planner_format") for a in "AC"}
r2 = []
for k, lab in [("prose", "Prose with no JSON (an answer, a number such as <code>\\boxed{54.7}</code>, or reasoning)"),
               ("json", 'The right JSON written as text, e.g. <code>{"name": "publish_method", "arguments": {"xc": "b3lyp", ...}}</code>')]:
    a, c = fk["A"][k], fk["C"][k]
    r2.append("<tr>" + td(lab, False) + td(str(a)) + td(str(c)) + td(f"{100 * (c - a) / a:+.0f}%") + "</tr>")
h2 = "".join(th(h, i > 0) for i, h in enumerate(["Reply that was not a tool call", "A", "C", "Change"]))


# 3. kinds of content mistake
def ckind(p):
    p = str(p)
    if p.startswith("parent: SMILES"):
        return "Parent SMILES is a formula or half-SMILES (<code>C6H6</code>, <code>C2H4</code>, <code>HC(=O)NH2</code>)"
    if "cannot have multiplicity" in p:
        return "Fragment with an impossible spin (e.g. an acyl radical given as an even-electron doublet)"
    if "SMILES" in p and "not valid" in p:
        return "Invalid radical SMILES (vinyl, dichloromethyl radicals …)"
    if "but the question asks" in p:
        return "Valid split of the wrong bond (e.g. allylic C–H when the question asks for methyl)"
    if "add up to" in p:
        return "Fragments do not add up to the parent"
    if "JSON array" in p:
        return "Fragments not a JSON array of species"
    return "Other"


ck = {a: Counter(ckind((e.get("problems") or ["?"])[0]) for e in E[a] if e["blame"] == "planner_content") for a in "AC"}
order = sorted(set(ck["A"]) | set(ck["C"]), key=lambda k: (k == "Other", -(ck["C"][k] + ck["A"][k])))
r3 = ["<tr>" + td(k, False) + td(str(ck["A"][k])) + td(str(ck["C"][k])) + "</tr>" for k in order]
r3.append("<tr>" + td("<strong>All content rejections</strong>", False) + td(f"<strong>{sum(ck['A'].values())}</strong>")
          + td(f"<strong>{sum(ck['C'].values())}</strong>") + "</tr>")
h3 = "".join(th(h, i > 0) for i, h in enumerate(["Why the board's validator rejected the content", "A", "C"]))


# 4. campaign outcomes: where no-plan campaigns stall, and wrong plans that got through
def campaigns(es):
    camp = defaultdict(list)
    for x in es:
        camp[(x["run"], x["task_index"])].append(x)
    stuck, wrong = Counter(), Counter()
    for k, xs in camp.items():
        o = xs[0]["oracle"]
        if o["error_class"] == "no_plan":
            acc = {x["module"] for x in xs if x["blame"] in OK}
            first = [m for m in ["target", "fragments", "method", "rationale"] if m not in acc]
            stuck[first[0] if first else "?"] += 1
        elif not o["correct"]:
            wrong[o["error_class"]] += 1
    return len(camp), stuck, wrong


cA, sA, wA = campaigns(E["A"])
cC, sC, wC = campaigns(E["C"])
r4 = []
for lab, a, c in [("No plan in budget: stalled at target", sA["target"], sC["target"]),
                  ("No plan in budget: stalled at fragments", sA["fragments"], sC["fragments"]),
                  ("No plan in budget: stalled at method or rationale", sA["method"] + sA["rationale"], sC["method"] + sC["rationale"]),
                  ("Wrong plan accepted: wrong bond", wA["wrong_bond"], wC["wrong_bond"]),
                  ("Wrong plan accepted: wrong parent", wA["wrong_parent"], wC["wrong_parent"])]:
    r4.append("<tr>" + td(lab, False) + td(str(a)) + td(str(c)) + "</tr>")
h4 = "".join(th(h, i > 0) for i, h in enumerate([f"Campaign outcome ({cA} campaigns per arm)", "A", "C"]))

nA, nC = len(E["A"]), len(E["C"])
SEC = f"""<div class="summary">
      <h3 id="failure-breakdown" {e()}>How training changed the way the model fails</h3>
      <p {e()} style="font-size:13px">Arm A (base nemotron-3-nano) against arm C (round-2 fine-tuned), both with the plan critic and retrieval off, so training is the only difference. Every first-try planner attempt over 3 repeats × 65 questions: {nA:,} attempts in A, {nC:,} in C, from the Flowcept episodes (<code>learn/eval/episodes_armA.jsonl</code>, <code>episodes_armC.jsonl</code>).</p>
      <div class="findings" style="margin:8px 0 12px">
        <div class="finding good"><h3 {e()}>More replies are tool calls</h3><p {e()}>Non-tool replies fell from 47% to 30% of target attempts and from 71% to 61% of fragments attempts. Most of the gain is prose answers that ignored the task ({fk['A']['prose']} → {fk['C']['prose']}); the commonest failure, the right JSON written as text, barely moved ({fk['A']['json']} → {fk['C']['json']}).</p></div>
        <div class="finding bad"><h3 {e()}>…but more of them are wrong</h3><p {e()}>Content rejections rose from {sum(ck['A'].values())} to {sum(ck['C'].values())}: formulas given as SMILES, impossible radical spins, invalid radical SMILES and valid splits of the wrong bond. The model now commits to an answer before getting the chemistry right, so acceptance is unchanged.</p></div>
        <div class="finding bad"><h3 {e()}>New: wrong parents get through</h3><p {e()}>Wrong plans accepted by the board went from {wA['wrong_parent']} to {wC['wrong_parent']} with the wrong parent (styrene meta C–H planned on benzene; hydrazine as <code>[N]N</code>). Ring positions stay unsolved in both: pyridine 2-/3-/4-position C–H is most of the wrong-bond plans.</p></div>
        <div class="finding"><h3 {e()}>Why, and what to try</h3><p {e()}>DPO paired correct replies against format failures, so the model learned "call the tool"; nothing rewarded a valid SMILES or spin. Two cheap fixes target the biggest buckets without retraining: accept JSON written as text ({fk['C']['json']} attempts in C) and turn a formula into a SMILES with RDKit before validating ({ck['C'][order[0]] if order else 0}).</p></div>
      </div>
      <div class="tw"><table style="min-width:760px"><thead><tr>{h1}</tr></thead><tbody>{"".join(r1)}</tbody></table></div>
      <div class="tw"><table style="min-width:760px;margin-top:12px"><thead><tr>{h2}</tr></thead><tbody>{"".join(r2)}</tbody></table></div>
      <div class="tw"><table style="min-width:760px;margin-top:12px"><thead><tr>{h3}</tr></thead><tbody>{"".join(r3)}</tbody></table></div>
      <div class="tw"><table style="min-width:640px;margin-top:12px"><thead><tr>{h4}</tr></thead><tbody>{"".join(r4)}</tbody></table></div>
      <p {e()} style="font-size:13px">Percentages are of that module's attempts. "Accepted" includes attempts the board accepted that the oracle later showed wrong. A no-plan campaign is counted at the first module that never got accepted. Content-rejection reasons are the validator's first problem for each attempt.</p>
    </div>
    """

m = re.search(r'<div class="summary">\s*<h3 data-e="e\d+">Caveats</h3>', page)
assert m
page = page[:m.start()] + SEC + page[m.start():]
toc_anchor = '<a href="#fair-comparison">Fair comparison</a>'
assert toc_anchor in page
css = "td .bad { color: var(--fail); font-weight: 600; }"
if css in page and "td.hl" not in page:
    page = page.replace(css, css + "\ntd.hl strong, td.hl { color: var(--pass); }", 1)
dup = [x for x, n in Counter(re.findall(r'data-e="(e\d+)"', page)).items() if n > 1]
assert not dup, dup
open(path, "w").write(page)
print("ok", len(page), "ids to", next(ids) - 1, fk, wA, wC, sA, sC)
