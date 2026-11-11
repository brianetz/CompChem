"""Build the 'Where the time goes' artifact section from stage_times.json."""
import html
import json
import sys

D = json.load(open(sys.argv[1]))
R = {r["label"]: r for r in D}

GROUPS = [
    ("nemotron-3-nano · one-shot planner, before the plan checks",
     [("Baseline", "water · no database"), ("A", "water · staged, empty DB"), ("B", "water · warm DB, retrieval on"),
      ("D", "water · warm DB, retrieval off"), ("C", "water · plain mode, empty DB"), ("V", "water · review check"),
      ("R1", "methanol · retrieval on"), ("R2", "methanol · retrieval off")]),
    ("nemotron-3-nano · with plan checks",
     [("F1", "methanol · critic off"), ("G3", "methanol · + LLM plan critic")]),
    ("Llama 3.1 405B · one-shot planner",
     [("L405", "methanol · before plan checks"), ("L405v2", "methanol · with plan checks"),
      ("F2", "methanol · critic off"), ("F3", "toluene · critic off"),
      ("G1", "methanol · + LLM plan critic"), ("G2", "toluene · + LLM plan critic")]),
    ("Llama 3.1 405B · plan board (planner and critic review every module)",
     [("B1", "toluene"), ("B2", "methanol"), ("B1v2", "toluene · formats + site labels"),
      ("B2v2", "methanol · formats + site labels"), ("NE1", "8 new molecules · retrieval on"),
      ("NE2", "8 new molecules · retrieval off")]),
]
ZONES = [("startup", "Startup", "var(--z-start)"), ("plan", "Planning", "var(--z-plan)"),
         ("critic", "Plan critic", "var(--z-critic)"), ("exec", "Execution", "var(--z-exec)"),
         ("result", "Result generation", "var(--z-result)"), ("claim_wait", "Waiting on claims", "var(--z-wait)")]
PER_ZONES = [z for z in ZONES if z[0] in ("plan", "critic", "exec", "result")]


def fmt_min(s):
    return f"{s / 60:.1f}"


def chart(kind):
    """kind='wall' (run wall clock vs time limit) or 'per' (average per campaign)."""
    X0, X1 = 270, 880
    vmax = 125 if kind == "wall" else 30           # minutes
    px = (X1 - X0) / vmax
    rows, y = [], 34
    ticks = range(0, vmax + 1, 20 if kind == "wall" else 5)
    out = []
    for gname, members in GROUPS:
        out.append(f'<text class="tg" x="16" y="{y + 14}">{html.escape(gname)}</text>')
        y += 24
        for label, desc in members:
            r = R[label]
            src = r["wall"] if kind == "wall" else r["per_campaign"]
            zones = ZONES if kind == "wall" else PER_ZONES
            n = r["campaigns"]
            out.append(f'<text class="tl" x="16" y="{y + 14}"><tspan class="tlb">{label}</tspan> · {html.escape(desc)} · {n}</text>')
            x = X0
            used = 0.0
            for key, name, color in zones:
                sec = src.get(key) or 0.0
                if sec <= 0:
                    continue
                w = max(1.0, sec / 60 * px)
                tip = f"{label}: {name.lower()} {fmt_min(sec)} min" + ("" if kind == "wall" else " per campaign")
                out.append(f'<rect x="{x:.1f}" y="{y}" width="{max(0.5, w - 1.5):.1f}" height="18" rx="2" fill="{color}"><title>{tip}</title></rect>')
                if w > 34:
                    out.append(f'<text class="tv" x="{x + w / 2 - 0.75:.1f}" y="{y + 13}" text-anchor="middle" style="fill:var(--z-ink-{key})">{fmt_min(sec)}</text>')
                x += w
                used += sec
            if kind == "wall":
                lim = (r["limit_min"] or 0) * 60
                lx = X0 + r["limit_min"] * px
                if used < lim:
                    out.append(f'<rect class="unused" x="{x:.1f}" y="{y + 0.5}" width="{lx - x - 1:.1f}" height="17" rx="2"><title>{label}: {fmt_min(lim - used)} min of the {r["limit_min"]}-min limit unused</title></rect>')
                over = used > lim
                out.append(f'<line class="{"limx" if over else "lim"}" x1="{lx:.1f}" x2="{lx:.1f}" y1="{y - 2}" y2="{y + 20}"/>')
                pct = round(100 * used / lim) if lim else 0
                note = (f"{fmt_min(used)} min · 2 blocks, first hit the {r['limit_min']}-min limit" if over
                        else f"{fmt_min(used)} min · {pct}% of {r['limit_min']} min")
                tx = max(x, lx) + 8 if over else lx + 8
                if over:
                    out.append(f'<text class="te" x="{tx:.1f}" y="{y + 13}">{note}</text>')
                else:
                    out.append(f'<text class="te" x="{tx:.1f}" y="{y + 13}">{note}</text>')
            else:
                out.append(f'<text class="te" x="{x + 6:.1f}" y="{y + 13}">{fmt_min(used)} min</text>')
            y += 26
        y += 6
    h = y + 6
    axis = []
    for t in ticks:
        tx = X0 + t * px
        axis.append(f'<line class="grid" x1="{tx:.1f}" x2="{tx:.1f}" y1="26" y2="{h - 10}"/>')
        axis.append(f'<text class="tk" x="{tx:.1f}" y="20" text-anchor="middle">{t}</text>')
    unit = "minutes of wall-clock time" if kind == "wall" else "minutes per campaign"
    axis.append(f'<text class="tk" x="{X0 - 8}" y="20" text-anchor="end">{unit}</text>')
    width = 1160 if kind == "wall" else 960
    return f'<svg viewBox="0 0 {width} {h}" role="img">' + "".join(axis) + "".join(out) + "</svg>"


def table():
    head = ("<tr><th>Run</th><th>Model · method</th><th class=num>Limit</th><th class=num>Startup</th>"
            "<th class=num>Planning</th><th class=num>Plan critic</th><th class=num>Execution</th>"
            "<th class=num>Result gen.</th><th class=num>Claim wait</th><th class=num>Total</th>"
            "<th class=num>Decision share</th><th class=num>Per campaign: plan · critic · exec</th></tr>")
    body = []
    for gname, members in GROUPS:
        for label, desc in members:
            r = R[label]; w = r["wall"]; p = r["per_campaign"]
            used = sum(w[k] for k, _, _ in ZONES)
            dec = (w["plan"] + w["critic"]) / used * 100 if used else 0
            cell = lambda s: f"{s / 60:.1f}" if s else "–"
            body.append(
                f"<tr><td>{label}<small>{r['run']}</small></td><td>{html.escape(r['model'])} · {html.escape(desc)}</td>"
                f"<td class=num>{r['limit_min']}</td><td class=num>{cell(w['startup'])}</td><td class=num>{cell(w['plan'])}</td>"
                f"<td class=num>{cell(w['critic'])}</td><td class=num>{cell(w['exec'])}</td><td class=num>{cell(w['result'])}</td>"
                f"<td class=num>{cell(w['claim_wait'])}</td><td class=num>{used / 60:.1f}</td><td class=num>{dec:.0f}%</td>"
                f"<td class=num>{p['plan'] / 60:.1f} · {p['critic'] / 60:.1f} · {p['exec']:.1f} s</td></tr>")
    return "<table><thead>" + head + "</thead><tbody>" + "".join(body) + "</tbody></table>"


def stats():
    s = {}
    for r in D:
        w = r["wall"]; used = sum(w[k] for k, _, _ in ZONES)
        s[r["label"]] = {"used_min": used / 60, "dec_pct": (w["plan"] + w["critic"]) / used * 100,
                         "exec_pct": w["exec"] / used * 100, "result_pct": w["result"] / used * 100,
                         "limit_pct": used / (r["limit_min"] * 60) * 100,
                         "per_total_min": sum(r["per_campaign"][k] for k in ("plan", "critic", "exec", "result")) / 60,
                         "critic_share": r["per_campaign"]["critic"] / max(1e-9, r["per_campaign"]["plan"] + r["per_campaign"]["critic"]) * 100}
    return s


if __name__ == "__main__":
    out = {"wall": chart("wall"), "per": chart("per"), "table": table()}
    json.dump(out, open(sys.argv[2], "w"))
    for k, v in stats().items():
        print(k, {a: round(b, 1) for a, b in v.items()})
