"""Build the 'How the errors came down' section: failure reasons per run, by stage."""
import html
import json
import sys

# From the page's "Why campaigns fail: per run" table (campaigns, then counts per reason),
# plus NE2 (8 of 8 correct, section 11).
REASONS = [("smiles", "Invalid SMILES", "var(--z-plan)"), ("spin", "Impossible spin", "var(--z-critic)"),
           ("format", "LLM output format", "var(--z-exec)"), ("lost", "Worker lost", "var(--z-result)"),
           ("wrong", "Wrong answer passed", "var(--z-wait)")]
STAGES = [
    ("1 · Water: results database and retrieval", "nemotron · 2026-10-14 to 10-19", [
        ("Baseline", "water · no database", 56, [21, 1, 1, 0, 0]),
        ("A", "water · empty DB", 56, [19, 0, 2, 0, 0]),
        ("B", "water · warm DB, retrieval", 56, [11, 0, 3, 0, 0]),
        ("D", "water · warm DB, no retrieval", 56, [19, 0, 2, 0, 1]),
        ("C", "water · plain mode", 56, [11, 0, 1, 0, 1]),
        ("V", "water · review check", 8, [1, 0, 0, 0, 0])]),
    ("2 · A new molecule, no plan checks", "nemotron · 10-19", [
        ("R1", "methanol · retrieval", 56, [23, 18, 1, 8, 1]),
        ("R2", "methanol · no retrieval", 56, [27, 19, 2, 1, 1])]),
    ("3 · A larger planner", "Llama 3.1 405B · 10-19", [
        ("L405", "methanol", 8, [0, 0, 6, 0, 1])]),
    ("4 · Plan checks, radical rules, stricter results critic", "10-20", [
        ("L405v2", "methanol · 405B", 8, [0, 0, 0, 0, 0]),
        ("F1", "methanol · nemotron", 56, [0, 2, 3, 0, 0]),
        ("F2", "methanol · 405B", 56, [1, 0, 0, 0, 1]),
        ("F3", "toluene · 405B", 8, [1, 3, 1, 0, 0])]),
    ("5 · Plus the LLM plan critic", "10-20", [
        ("G3", "methanol · nemotron", 24, [1, 0, 3, 0, 0]),
        ("G1", "methanol · 405B", 16, [1, 0, 1, 0, 0]),
        ("G2", "toluene · 405B", 8, [0, 3, 0, 0, 1])]),
    ("6 · Plan board", "405B · 10-20", [
        ("B1", "toluene", 8, [0, 0, 0, 0, 3]),
        ("B2", "methanol", 8, [0, 0, 2, 0, 0])]),
    ("7 · Plan board v2: sample formats and site labels", "405B · 10-21", [
        ("B1v2", "toluene", 8, [0, 0, 0, 0, 0]),
        ("B2v2", "methanol", 8, [0, 0, 0, 0, 0])]),
    ("8 · New molecules on board v2", "405B · 10-21", [
        ("NE1", "8 molecules · retrieval", 8, [0, 0, 0, 0, 0]),
        ("NE2", "8 molecules · no retrieval", 8, [0, 0, 0, 0, 0])]),
]


def chart():
    X0, X1 = 270, 830
    px = (X1 - X0) / 100
    y, out = 34, []
    for name, sub, runs in STAGES:
        n_all = sum(r[2] for r in runs)
        bad_all = sum(sum(r[3]) for r in runs)
        out.append(f'<text class="tg" x="16" y="{y + 14}">{html.escape(name)}</text>'
                   f'<text class="tk" x="{X1 + 8}" y="{y + 14}">{html.escape(sub)} · stage: {bad_all}/{n_all} not correct ({100 * bad_all / n_all:.0f}%)</text>')
        y += 24
        for label, desc, n, counts in runs:
            out.append(f'<text class="tl" x="16" y="{y + 14}"><tspan class="tlb">{label}</tspan> · {html.escape(desc)} · {n}</text>')
            x = X0
            for (key, rname, color), c in zip(REASONS, counts):
                if not c:
                    continue
                w = 100 * c / n * px
                out.append(f'<rect x="{x:.1f}" y="{y}" width="{max(0.8, w - 1.5):.1f}" height="18" rx="2" fill="{color}">'
                           f'<title>{label}: {rname.lower()} {c} of {n} campaigns ({100 * c / n:.0f}%)</title></rect>')
                if w > 22:
                    ink = {"smiles": "plan", "spin": "critic", "format": "exec", "lost": "result", "wrong": "claim_wait"}[key]
                    out.append(f'<text class="tv" x="{x + w / 2 - 0.75:.1f}" y="{y + 13}" text-anchor="middle" style="fill:var(--z-ink-{ink})">{c}</text>')
                x += w
            bad = sum(counts)
            if bad < n:
                out.append(f'<rect class="okseg" x="{x:.1f}" y="{y + 0.5}" width="{X1 - x - 1:.1f}" height="17" rx="2">'
                           f'<title>{label}: {n - bad} of {n} campaigns correct</title></rect>')
            out.append(f'<text class="te" x="{X1 + 8}" y="{y + 13}">{bad}/{n} not correct · {100 * bad / n:.0f}%</text>')
            y += 26
        y += 6
    h = y + 4
    axis = []
    for t in range(0, 101, 25):
        tx = X0 + t * px
        axis.append(f'<line class="grid" x1="{tx:.1f}" x2="{tx:.1f}" y1="26" y2="{h - 8}"/>'
                    f'<text class="tk" x="{tx:.1f}" y="20" text-anchor="middle">{t}%</text>')
    axis.append(f'<text class="tk" x="{X0 - 8}" y="20" text-anchor="end">share of the run\'s campaigns</text>')
    return f'<svg viewBox="0 0 1160 {h}" role="img">' + "".join(axis) + "".join(out) + "</svg>"


def stage_table():
    rows = []
    for name, sub, runs in STAGES:
        n = sum(r[2] for r in runs)
        c = [sum(r[3][i] for r in runs) for i in range(5)]
        bad = sum(c)
        cells = "".join(f"<td class=num>{v or '–'}</td>" for v in c)
        rows.append(f"<tr><td>{html.escape(name)}<small>{html.escape(sub)} · {', '.join(r[0] for r in runs)}</small></td>"
                    f"<td class=num>{n}</td>{cells}<td class=num><strong>{bad}</strong> <small>{100 * bad / n:.0f}%</small></td></tr>")
    head = ("<tr><th>Stage</th><th class=num>Campaigns</th>" + "".join(f"<th class=num>{r[1]}</th>" for r in REASONS)
            + "<th class=num>Not correct</th></tr>")
    return "<table><thead>" + head + "</thead><tbody>" + "".join(rows) + "</tbody></table>"


if __name__ == "__main__":
    json.dump({"chart": chart(), "table": stage_table()}, open(sys.argv[1], "w"))
    for name, sub, runs in STAGES:
        n = sum(r[2] for r in runs); bad = sum(sum(r[3]) for r in runs)
        print(f"{name}: {bad}/{n} = {100 * bad / n:.1f}%")
    print("total", sum(sum(sum(r[3]) for r in s[2]) for s in STAGES), sum(sum(r[2] for r in s[2]) for s in STAGES))
