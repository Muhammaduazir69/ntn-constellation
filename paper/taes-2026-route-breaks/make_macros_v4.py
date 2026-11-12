#!/usr/bin/env python3
"""Write ../macros_p7.tex from out/: closed-form law, real-shell check and the ns-3 campaign.
Every number in the manuscript comes from here."""
import json, os
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
M = {}
TAG = {"starlink": "Starlink", "kuiper": "Kuiper", "oneweb": "Oneweb"}
PTAG = {"lon-nyc": "Lon", "fra-sin": "Fra", "syd-lax": "Syd"}


def f(x, d=1):
    return f"{x:,.{d}f}" if abs(x) >= 1e4 else f"{x:.{d}f}"

from p7util import dedupe


# ---- closed form
L = json.load(open(os.path.join(OUT, "law.json")))
G = json.load(open(os.path.join(OUT, "grid_allF.json")))
for s, t in TAG.items():
    d = L[s]
    M[f"GridLinks{t}"] = f(d["grid_links"], 0)
    M[f"GridDmax{t}"], M[f"GridDmin{t}"] = f(d["grid_dmax_km"], 0), f(d["grid_dmin_km"], 0)
    M[f"GridAllF{t}"] = f(G[s], 0)
    M[f"Rcap{t}"] = f(d["rcap_km"], 0)
    M[f"AfDeg{t}"], M[f"AfDegS{t}"] = f(d["af_degree"], 2), f(d["af_degree_sampled"], 2)
    M[f"AfEv{t}"], M[f"AfEvS{t}"] = f(d["af_events"], 0), f(d["af_events_sampled"], 0)
    M[f"AfEpochs{t}"] = f(d["af_epochs"], 0)
    M[f"AfGap{t}"] = f(d["af_epoch_gap_ms"], 1)
    # longest pass with the link dropped at 23 deg (2 deg hysteresis below the 25 deg mask), over all sites
    from law_tables import pass_max  # noqa: E402
    import isl_core as _C  # noqa: E402
    P_, S_, F_, i_, h_ = d["P"], d["S"], d["F"], d["inc"], d["alt_km"] * 1e3
    sh_ = _C.walker(P_, S_, F_, i_, h_, s == "oneweb")
    pm = max(pass_max(sh_, lat, 23.0) for lat in (51.51, 40.71, 50.11, 1.35, -33.87, 34.05))
    M[f"Pass{t}"] = f(pm, 0)
    M[f"Floor{t}"] = f(3600 / pm, 1)
    M[f"Period{t}"] = f(d["period_s"] / 60, 1)

PC = json.load(open(os.path.join(OUT, "pair_check.json")))
M["PairN"], M["PairErrNm"] = f(PC["pairs"], 0), f(PC["max_err_m"] * 1e9, 0)
R = json.load(open(os.path.join(OUT, "real_shell.json")))
M["RealSats"] = f(R["shell_sats"], 0); M["RealPairs"] = f(R["pairs"], 0)
M["RealErrMed"], M["RealErrNf"], M["RealErrMax"] = (f(R[k], 1) for k in
                                                    ("range_err_km_median", "range_err_km_p95", "range_err_km_max"))
M["RealEpMed"], M["RealEpNf"] = f(R["epoch_err_s_median"], 0), f(R["epoch_err_s_p95"], 0)

# ---- ns-3 campaign v4 (out/ns3_routes4.jsonl)
RNG = np.random.default_rng(7)
rows = [json.loads(l) for l in open(os.path.join(OUT, "ns3_routes4.jsonl"))] if os.path.exists(os.path.join(OUT, "ns3_routes4.jsonl")) else []
rows = [r for r in rows if "error" not in r]
for r in rows:
    r["cfgd"] = dict(r["cfg"])
M["NsRuns"] = f(sum(r["part"] == "main" for r in rows), 0)
M["NsRunsAll"] = f(len(rows), 0)
POL = [("SnapSlow", dict(policy="snap", period=1.0)), ("SnapFast", dict(policy="snap", period=0.1)),
       ("Od", dict(policy="ondemand")), ("Asr", dict(policy="asr", window=30.0)),
       ("AsrLong", dict(policy="asr", window=120.0)), ("Pshort", dict(policy="life", eps=0.0)),
       ("LifeOne", dict(policy="life", eps=0.01)), ("LifeFive", dict(policy="life", eps=0.05)),
       ("LifeTen", dict(policy="life", eps=0.1)), ("Hold", dict(policy="life", eps=0.05, hold=1)),
       ("Mrv", dict(policy="mrv"))]


def match(r, spec, part="main"):
    d = r["cfgd"]
    if r["part"] != part:
        return False
    if spec.get("hold", 0) != d.get("hold", 0):
        return False
    return all(abs(d.get(k, -9) - v) < 1e-9 if isinstance(v, float) else d.get(k) == v for k, v in spec.items() if k != "hold")


def cells(spec, part="main"):
    out = {}
    for r in rows:
        if match(r, spec, part):
            out.setdefault((r["shell"], r["pair"]), []).append(r)
    return out


def strat(spec, key, scale=1.0, part="main", B=2000):
    """Mean stratified over (shell, pair) cells, with a cluster bootstrap over start instants in each cell."""
    C = cells(spec, part)
    if not C:
        return None
    vals = [np.array([r[key] * scale for r in v]) for v in C.values()]
    m = np.mean([v.mean() for v in vals])
    bs = [np.mean([RNG.choice(v, len(v)).mean() for v in vals]) for _ in range(B)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return m, (hi - lo) / 2


def ratio(spec_a, spec_b, key, part="main", B=2000):
    """Stratified mean of b over stratified mean of a, paired by run, with a cluster bootstrap."""
    A, Bc = cells(spec_a, part), cells(spec_b, part)
    ks = [k for k in A if k in Bc]
    pairs = []
    for k in ks:
        da = {r["start"]: r[key] for r in A[k]}; db = {r["start"]: r[key] for r in Bc[k]}
        pairs.append(np.array([(da[s], db[s]) for s in da if s in db]))
    if not pairs:
        return None
    est = np.mean([p[:, 1].mean() for p in pairs]) / np.mean([p[:, 0].mean() for p in pairs])
    bs = []
    for _ in range(B):
        res = [p[RNG.integers(0, len(p), len(p))] for p in pairs]
        bs.append(np.mean([q[:, 1].mean() for q in res]) / np.mean([q[:, 0].mean() for q in res]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return est, lo, hi


for tag, spec in POL:
    for key, name, sc, dec in [("loss", "Loss", 100, 2), ("changes_per_h", "Chg", 1, 0), ("forced_per_h", "Forced", 1, 0),
                               ("dijkstra_per_h", "Dij", 1, 0), ("predictions_per_h", "Pred", 1, 0),
                               ("stretch_mean", "Str", 100, 1), ("stretch_p99", "StrNn", 100, 1)]:
        s = strat(spec, key, sc)
        if s:
            M[f"{name}{tag}"], M[f"{name}{tag}Ci"] = f(s[0], dec), f(s[1], dec)
    for sh, t in TAG.items():
        for key, name, sc, dec in [("stretch_mean", "Str", 100, 1), ("stretch_p99", "StrNn", 100, 1), ("changes_per_h", "Chg", 1, 0)]:
            v = [r[key] * sc for (s_, p_), rr in cells(spec).items() if s_ == sh for r in rr]
            if v:
                M[f"{name}{tag}{t}"] = f(np.mean(v), dec)
    M[f"Lost{tag}"] = str(int(sum(r["sent"] - r["received"] for rr in cells(spec).values() for r in rr)))
M["DelayLonNycStarlink"] = f(np.mean([r["delay_ms_mean"] for r in cells(dict(policy="snap", period=1.0)).get(("starlink", "lon-nyc"), [])] or [0]), 1)
for a_tag, a_spec, b_tag, b_spec, key in [("LifeFive", dict(policy="life", eps=0.05), "Asr", dict(policy="asr", window=30.0), "changes_per_h"),
                                          ("LifeFive", dict(policy="life", eps=0.05), "Od", dict(policy="ondemand"), "changes_per_h"),
                                          ("LifeFive", dict(policy="life", eps=0.05), "SnapSlow", dict(policy="snap", period=1.0), "changes_per_h"),
                                          ("LifeFive", dict(policy="life", eps=0.05), "Asr", dict(policy="asr", window=30.0), "dijkstra_per_h"),
                                          ("LifeFive", dict(policy="life", eps=0.05), "Pshort", dict(policy="life", eps=0.0), "changes_per_h"),
                                          ("LifeFive", dict(policy="life", eps=0.05), "AsrLong", dict(policy="asr", window=120.0), "changes_per_h"),
                                          ("LifeTen", dict(policy="life", eps=0.10), "Pshort", dict(policy="life", eps=0.0), "changes_per_h")]:
    r_ = ratio(a_spec, b_spec, key)
    if r_:
        nm = f"Ratio{b_tag}{'Chg' if key == 'changes_per_h' else 'Dij'}" + ("" if a_tag == "LifeFive" else f"Vs{a_tag}")
        M[nm], M[nm + "Lo"], M[nm + "Hi"] = f(r_[0], 2), f(r_[1], 2), f(r_[2], 2)

# Proposition 2 with acquisition and in-flight terms, on every snapshot and on-demand run
pred, meas = [], []
for r in rows:
    d = r["cfgd"]
    if r["part"] not in ("main", "tacq") or d.get("policy") not in ("snap", "ondemand"):
        continue
    lam = r["breaks_per_h"] / 3600.0
    wait = d["period"] / 2 if d["policy"] == "snap" else 0.05
    pred.append(lam * (wait + d.get("tAcq", 0.0) + 0.5 * r["delay_ms_mean"] / 1e3)); meas.append(r["loss"])  # W_b ~ W/2
if pred:
    pred, meas = np.array(pred), np.array(meas)
    ok = pred > 0
    rel = np.abs(meas[ok] - pred[ok]) / pred[ok]
    M["StaleN"] = f(int(ok.sum()), 0)
    M["StaleRelErr"], M["StaleRelErrNf"] = f(100 * np.median(rel), 0), f(100 * np.percentile(rel, 95), 0)
    json.dump(dict(pred=pred.tolist(), meas=meas.tolist()), open(os.path.join(OUT, "stale_check.json"), "w"))
M["IslDownMain"] = f(sum(r["isl_down_events"] for r in rows if r["part"] == "main"), 0)

# stress: acquisition time
for ta, tt in ((0.0, "Zero"), (3.0, "Three")):
    for tag, spec in [("SnapSlow", dict(policy="snap", period=1.0)), ("Od", dict(policy="ondemand")),
                      ("Asr", dict(policy="asr", window=30.0)), ("LifeFive", dict(policy="life", eps=0.05))]:
        rr = [r for r in rows if r["part"] == "tacq" and abs(r["cfgd"].get("tAcq", -1) - ta) < 1e-9 and match(dict(r, part="main"), spec)]
        if rr:
            M[f"LossAcq{tt}{tag}"] = f(100 * np.mean([r["loss"] for r in rr]), 2)
# stress: guard and ephemeris error
for gd, gt in ((0.1, "Tenth"), (0.5, "Half"), (2.0, "Two")):
    for sg, st_ in ((0.0, "Zero"), (5.0, "Five"), (20.0, "Twenty")):
        if gd == 0.5 and sg == 0.0:
            rr = [r for r in rows if r["part"] == "main" and r["pair"] == "lon-nyc" and r["shell"] in ("starlink", "oneweb")
                  and match(r, dict(policy="life", eps=0.05)) and r["start"] in (0, 21600, 43200)]
        else:
            rr = [r for r in rows if r["part"] == "pred" and r["cfgd"].get("policy") == "life"
                  and abs(r["cfgd"].get("guard", -1) - gd) < 1e-9 and abs(r["cfgd"].get("predSigmaKm", -1) - sg) < 1e-9]
        if rr:
            M[f"LossG{gt}S{st_}"] = f(100 * np.mean([r["loss"] for r in rr]), 3)
            M[f"BrkG{gt}S{st_}"] = f(np.mean([r["breaks_per_h"] for r in rr]), 1)
            M[f"ChgG{gt}S{st_}"] = f(np.mean([r["changes_per_h"] for r in rr]), 0)
for sg, st_ in ((5.0, "Five"), (20.0, "Twenty")):
    rr = [r for r in rows if r["part"] == "pred" and r["cfgd"].get("policy") == "asr" and abs(r["cfgd"].get("predSigmaKm", -1) - sg) < 1e-9]
    if rr:
        M[f"LossAsrS{st_}"] = f(100 * np.mean([r["loss"] for r in rr]), 3)
# stress: cross-plane latitude cutoff on OneWeb
for tag, spec in [("SnapSlow", dict(policy="snap", period=1.0)), ("Od", dict(policy="ondemand")),
                  ("Asr", dict(policy="asr", window=30.0)), ("LifeFive", dict(policy="life", eps=0.05))]:
    rr = [r for r in rows if r["part"] == "cross" and match(dict(r, part="main"), spec)]
    if rr:
        M[f"CrossLoss{tag}"] = f(100 * np.mean([r["loss"] for r in rr]), 2)
        M[f"CrossBrk{tag}"] = f(np.mean([r["breaks_per_h"] for r in rr]), 0)
        M[f"CrossChg{tag}"] = f(np.mean([r["changes_per_h"] for r in rr]), 0)
        M[f"CrossIslDown{tag}"] = f(np.mean([r["isl_down_events"] for r in rr]), 0)
# offline optimum of (P1)
if os.path.exists(os.path.join(OUT, "offline_opt.json")):
    O = json.load(open(os.path.join(OUT, "offline_opt.json")))
    for sh, t in TAG.items():
        g_ = [v["greedy_per_h"] for k, v in O.items() if k.startswith(sh) and v["greedy_per_h"] is not None]
        o_ = [v["opt_per_h"] for k, v in O.items() if k.startswith(sh) and v["opt_per_h"] is not None]
        if g_ and o_:
            M[f"Greedy{t}"], M[f"Opt{t}"] = f(np.mean(g_), 1), f(np.mean(o_), 1)
            M[f"GreedyGap{t}"] = f(100 * (np.mean(g_) / np.mean(o_) - 1), 0)

with open(os.path.join(HERE, "..", "macros_p7.tex"), "w") as fh:
    fh.write("% generated by sim/make_macros.py; do not edit\n")
    for k, v in sorted(M.items()):
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print(len(M), "macros")
