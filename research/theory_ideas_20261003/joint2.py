"""J2 extension sweep (PROTOCOL_J2.md); reuses joint.evaluate."""
import itertools, json
from pathlib import Path
import joint, sweeps_abc as m
ROOT = Path(__file__).resolve().parent
J, sumw = m.load()
grid = [("onsager", q, l, True, t0, S) for q, l, t0, S in itertools.product((8.0, 10.0, 12.0), (1.05, 1.2, 1.4), (8.0, 10.0, 12.0), (240, 360, 560, 960))]
grid += [("plain", q, 0.0, False, t0, S) for q, t0, S in itertools.product((10.0, 12.0), (12.0, 20.0, 30.0), (560, 960, 1560))]
pilot = joint.evaluate(J, sumw, grid, 256, 50001)
sel = {}
for E in m.E_LIST:
    for fam in ("plain", "onsager"):
        b = min((r for r in pilot if r["mode"] == fam), key=lambda r: r["score"][E])
        sel[f"{fam}_E{E}"] = (b["mode"], b["q"], b["par"], b["ramp"], b["T0"], b["S"])
hold = joint.evaluate(J, sumw, sorted(set(sel.values())), 1024, 50002)
h = {(r["mode"], r["q"], r["par"], r["ramp"], r["T0"], r["S"]): r for r in hold}
prev = json.loads((ROOT / "J_results.json").read_text())["summary"]
summary = {}
for E in m.E_LIST:
    row = {fam: dict(cfg=sel[f"{fam}_E{E}"], p=h[sel[f"{fam}_E{E}"]]["p"], t_ms=h[sel[f"{fam}_E{E}"]]["t_ms"], tts_ms=h[sel[f"{fam}_E{E}"]]["tts"][E]) for fam in ("plain", "onsager")}
    row["J_onsager_tts_ms"] = prev[str(E)]["onsager"]["tts_ms"]; row["J_plain_tts_ms"] = prev[str(E)]["plain"]["tts_ms"]
    summary[E] = row; print("SUMMARY", E, json.dumps(row), flush=True)
(ROOT / "J2_results.json").write_text(json.dumps(dict(pilot=pilot, selection={k: list(v) for k, v in sel.items()}, holdout=hold, summary=summary), indent=1) + "\n")
