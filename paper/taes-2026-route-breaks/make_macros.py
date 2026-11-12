#!/usr/bin/env python3
"""Write ../macros_p7.tex from out/: closed-form law, real-shell check and the ns-3 campaign out/p7_v8.jsonl (one binary).
Every number in the manuscript comes from here."""
import json, os
import numpy as np
from scipy import stats as st
from p7util import *

M = {}
TAG = {"starlink": "Starlink", "kuiper": "Kuiper", "oneweb": "Oneweb"}
CT = {("starlink", "lon-nyc"): "StL", ("starlink", "syd-lax"): "StS", ("kuiper", "lon-nyc"): "KuL",
      ("kuiper", "syd-lax"): "KuS", ("oneweb", "lon-nyc"): "OnL", ("oneweb", "syd-lax"): "OnS"}


def f(x, d=1):
    return f"{x:,.{d}f}" if abs(x) >= 1e4 else f"{x:.{d}f}"


# ---- closed form
L = json.load(open(os.path.join(OUT, "law.json")))
G = json.load(open(os.path.join(OUT, "grid_allF.json")))
from law_tables import pass_max
import isl_core as _C
FLOOR = {}
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
    sh_ = _C.walker(d["P"], d["S"], d["F"], d["inc"], d["alt_km"] * 1e3, s == "oneweb")
    pm = max(pass_max(sh_, lat, 23.0) for lat in (51.51, 40.71, -33.87, 34.05))
    M[f"Pass{t}"] = f(pm, 0)
    FLOOR[s] = 3600 / pm
    M[f"Floor{t}"] = f(FLOOR[s], 1)
    M[f"Period{t}"] = f(d["period_s"] / 60, 1)
PC = json.load(open(os.path.join(OUT, "pair_check.json")))
M["PairN"], M["PairErrNm"] = f(PC["pairs"], 0), f(PC["max_err_m"] * 1e9, 0)
RS = json.load(open(os.path.join(OUT, "real_shell.json")))
M["RealSats"] = f(RS["shell_sats"], 0); M["RealPairs"] = f(RS["pairs"], 0)
M["RealErrMed"], M["RealErrNf"], M["RealErrMax"] = (f(RS[k], 1) for k in ("range_err_km_median", "range_err_km_p95", "range_err_km_max"))
M["RealEpMed"], M["RealEpNf"] = f(RS["epoch_err_s_median"], 0), f(RS["epoch_err_s_p95"], 0)
M["RealEvents"] = f(RS["events_true"], 0)
# Lemma 1 against the library propagator (m10)
errs = {}
for l in open(os.path.join(OUT, "law_check.jsonl")):
    r = json.loads(l); d = np.array(r["range"]) / 1e3
    x = L[r["shell"]]
    sh = _C.walker(x["P"], x["S"], x["F"], x["inc"], x["alt_km"] * 1e3, r["shell"] == "oneweb")
    a, b = r["a"], r["b"]
    M_, N_ = _C.MN(sh["Om"][b] - sh["Om"][a], sh["u0"][b] - sh["u0"][a], sh["i"])
    dx, dn = sh["r"] * np.sqrt(2 - 2 * (M_ - abs(N_))) / 1e3, sh["r"] * np.sqrt(max(2 - 2 * (M_ + abs(N_)), 0)) / 1e3
    errs.setdefault(r["shell"], []).append(max(abs(d.max() - dx), abs(d.min() - dn)))
for s, t in TAG.items():
    M[f"LawPropErr{t}"] = f(max(errs[s]), 1)
M["LawPropErrMax"] = f(max(max(v) for v in errs.values()), 1)
TC = json.load(open(os.path.join(OUT, "track_check.json")))
for s, t in TAG.items():
    M[f"TrackMed{t}"], M[f"TrackNf{t}"], M[f"TrackMax{t}"] = f(TC[s]["med_abs"], 1), f(TC[s]["p95_abs"], 1), f(TC[s]["max_abs"], 1)
M["TrackNfMax"] = f(max(TC[s]["p95_abs"] for s in TAG), 1)
M["TrackMedMax"] = f(max(TC[s]["med_abs"] for s in TAG), 1)
TN = json.load(open(os.path.join(OUT, "track_noise.json")))
for s, t in (("starlink", "Starlink"), ("oneweb", "Oneweb")):
    for tag, tt in (("three", "Three"), ("lsq", "Lsq")):
        for sg, st_ in ((0.01, "Hund"), (0.05, "Five")):
            e = TN[f"{s}_{tag}_{sg}"]
            M[f"Tn{t}{tt}{st_}Med"], M[f"Tn{t}{tt}{st_}Nf"] = f(e["med"], 1), f(e["p95"], 0)
RSH = json.load(open(os.path.join(OUT, "real_shells.json")))
for s, t in (("starlink", "Starlink"), ("oneweb", "Oneweb")):
    M[f"ShellSats{t}"], M[f"ShellPlanes{t}"] = f(RSH[s]["sats"], 0), f(RSH[s]["planes"], 0)
    M[f"ShellAlt{t}"] = f(RSH[s]["alt_km"], 0)

# ---- campaign
R = load()
sha = {r["bin_sha"] for r in R}.pop()
M["BinSha"] = sha[:12]
M["NsRunsAll"] = f(len(R), 0)
M["NsRunsMain"] = f(sum(r["part"] == "main" for r in R), 0)
M["NsRunsData"] = f(sum(not r["cfgd"].get("noData") for r in R), 0)
M["NsRunsCtrl"] = f(sum(bool(r["cfgd"].get("noData")) for r in R), 0)
M["NsRunsStress"] = f(sum(r["part"] in ("tacq", "bbm", "pred", "cross") for r in R), 0)
MAINX = ("main", "extra")
E24 = lambda e: dict(policy="elevb", eps=e, trigger=24.0)
POL = [("SnapSlow", dict(policy="snap", period=1.0)), ("Od", dict(policy="ondemand")), ("Elev", dict(policy="elev")),
       ("ElevbFive", dict(policy="elevb", eps=0.05)), ("ElevbTen", dict(policy="elevb", eps=0.10)), ("ElevtFive", E24(0.05)), ("ElevtTen", E24(0.10)),
       ("TrkFive", dict(policy="track", eps=0.05)), ("TrkTen", dict(policy="track", eps=0.10)),
       ("Asr", dict(policy="asr", window=30.0)), ("AsrLong", dict(policy="asr", window=120.0)), ("AsrVlong", dict(policy="asr", window=300.0)),
       ("AsrkLong", dict(policy="asrk", window=120.0)), ("Alpr", dict(policy="alpr")), ("AlprR", dict(policy="alpr", proactive=0)),
       ("Isasr", dict(policy="isasr", gamma=1.0)), ("Mrv", dict(policy="mrv")), ("Pshort", dict(policy="life", eps=0.0)),
       ("LifeFive", dict(policy="life", eps=0.05)), ("LifeTen", dict(policy="life", eps=0.10)), ("HoldFive", dict(policy="life", eps=0.05, hold=1)),
       ("HoldTwo", dict(policy="life", eps=0.02, hold=1)), ("HoldTen", dict(policy="life", eps=0.10, hold=1))]
SP = dict(POL)
for tag, spec in POL:
    for key, name, sc, dec in [("loss", "Loss", 100, 2), ("changes_per_h", "Chg", 1, 1), ("acq_per_h", "Acq", 1, 1), ("isl_only_per_h", "IslOnly", 1, 2),
                               ("forced_per_h", "Forced", 1, 1), ("dijkstra_per_h", "Dij", 1, 0), ("stretch_mean", "Str", 100, 1),
                               ("stretch_p99", "StrNn", 100, 1)]:
        s = strat(R, spec, key, sc, "main" if key == "loss" else MAINX)
        if s:
            M[f"{name}{tag}"], M[f"{name}{tag}Ci"] = f(s[0], dec), f(s[1], dec)
    M[f"Lost{tag}"] = str(int(sum(r["sent"] - r["received"] for r in R if match(r, spec, "main"))))
    for c, ct in CT.items():
        a, b, l_ = (by_start(R, spec, c, "main" if k == "loss" else MAINX, k, x) for k, x in (("changes_per_h", 1), ("stretch_mean", 100), ("loss", 100)))
        if a:
            M[f"Chg{tag}{ct}"], M[f"Str{tag}{ct}"], M[f"Loss{tag}{ct}"] = f(np.mean(list(a.values())), 1), f(np.mean(list(b.values())), 1), f(np.mean(list(l_.values())) if l_ else float("nan"), 2)
    for sh, t in TAG.items():
        v = [r["changes_per_h"] for r in R if r["shell"] == sh and match(r, spec, MAINX)]
        if v:
            M[f"Chg{tag}{t}"] = f(np.mean(v), 1)
M["NStarts"] = str(len({r["start"] for r in R if r["part"] == "main"}))
M["NStartsOneweb"] = str(len(by_start(R, SP["LifeFive"], ("oneweb", "lon-nyc"), MAINX)))
M["IslDownMain"] = f(sum(r["isl_down_events"] for r in R if r["part"] in MAINX), 0)
M["NMainX"] = f(sum(r["part"] in MAINX for r in R), 0)
# changes that replace ground links; ends
life_rows = [r for r in R if match(r, SP["LifeFive"], MAINX)]
allrows = [r for r in R if r["part"] in MAINX and r["changes_per_h"] > 0]
M["ShareIslOnly"] = f(100 * sum(r["isl_only_per_h"] for r in allrows) / sum(r["changes_per_h"] for r in allrows), 1)
M["ShareBothEnds"] = f(100 * sum(r["ingress_per_h"] + r["egress_per_h"] - r["acq_per_h"] for r in allrows) / sum(r["changes_per_h"] for r in allrows), 0)
od = [r for r in R if match(r, SP["Od"], MAINX)]
M["BreaksInOd"], M["BreaksOutOd"] = f(np.mean([r["breaks_in_per_h"] for r in od]), 1), f(np.mean([r["breaks_out_per_h"] for r in od]), 1)
M["ResInLife"], M["ResOutLife"] = f(np.mean([r["res_in_mean_s"] for r in life_rows]), 0), f(np.mean([r["res_out_mean_s"] for r in life_rows]), 0)
M["ResInOd"], M["ResOutOd"] = f(np.mean([r["res_in_mean_s"] for r in od]), 0), f(np.mean([r["res_out_mean_s"] for r in od]), 0)

# seeds: starts change geometry
gsl = {}
for r in od:
    gsl.setdefault((r["shell"], r["pair"]), []).append(r["changes_per_h"])
M["OdChgSpreadOnL"] = f(np.ptp(gsl[("oneweb", "lon-nyc")]), 1)
v = [r["changes_per_h"] for r in R if match(r, SP["LifeFive"], "main") and r["shell"] == "starlink" and r["pair"] == "lon-nyc"]
M["LifeStLSpread"] = f(np.ptp(v), 1); M["LifeStLMean"] = f(np.mean(v), 1)
v = list(by_start(R, SP["LifeFive"], ("oneweb", "lon-nyc"), MAINX, "stretch_mean", 100).values())
M["LifeOnLStrLo"], M["LifeOnLStrHi"] = f(min(v), 1), f(max(v), 1)

from frontier_stats import frontier_macros
frontier_macros(R, M, f, SP, CT)

# control plane against data plane
bad = 0
for r in R:
    if r["part"] == "check":
        sp = {k: v for k, v in r["cfgd"].items() if k not in ("noData", "tAcq")}
        v = by_start(R, sp, (r["shell"], r["pair"]), "main").get(r["start"])
        bad += (v is None) or abs(v - r["changes_per_h"]) > 1e-6
M["CheckBad"] = str(bad)
# trigger sweep
for trig, tt in ((23.3, "TwentyThreeA"), (23.5, "TwentyThreeB"), (24.0, "TwentyFour"), (25.0, "TwentyFive"), (27.0, "TwentySeven"), (30.0, "Thirty")):
    sp = dict(policy="elevb", eps=0.05) if trig == 25.0 else dict(policy="elevb", eps=0.05, trigger=trig)
    a = strat(R, sp, "acq_per_h", 1, ("main", "trig")); b = strat(R, sp, "stretch_mean", 100, ("main", "trig")); c_ = strat(R, sp, "forced_per_h", 1, ("main", "trig"))
    if a:
        M[f"ChgTrig{tt}"], M[f"StrTrig{tt}"], M[f"ForcedTrig{tt}"] = f(a[0], 1), f(b[0], 1), f(c_[0], 2)
a = strat(R, dict(policy="elevb", eps=0.05, trigger=24.0, floorGrid=1.0), "acq_per_h", 1, "trig"); b = strat(R, dict(policy="elevb", eps=0.05, trigger=24.0, floorGrid=1.0), "stretch_mean", 100, "trig")
if a:
    M["ChgGridOne"], M["StrGridOne"] = f(a[0], 1), f(b[0], 1)

# M5: realized stretch of the selection budget, enforced ablation
lr = [r for r in R if match(r, SP["LifeFive"], MAINX)]
M["LifeStrRunMax"] = f(100 * max(r["stretch_mean"] for r in lr), 1)
M["LifeStrPnnRunMax"] = f(100 * max(r["stretch_p99"] for r in lr), 1)
worst = max(lr, key=lambda r: r["stretch_p99"])
M["LifeWorstShell"] = {"starlink": "Starlink", "kuiper": "Kuiper", "oneweb": "OneWeb"}[worst["shell"]]
M["LifeOverFive"] = f(100 * np.mean([r["over5"] for r in lr]), 0)
M["LifeOverTen"] = f(100 * np.mean([r["over10"] for r in lr]), 0)
hr = [r for r in R if match(r, SP["HoldFive"], "main")]
if hr:
    M["HoldStrRunMax"] = f(100 * max(r["stretch_mean"] for r in hr), 1)
    M["HoldStrPnnRunMax"] = f(100 * max(r["stretch_p99"] for r in hr), 1)
    M["HoldOverFive"] = f(100 * np.mean([r["over5"] for r in hr]), 0)
    M["HoldPerH"] = f(np.mean([r["hold_per_h"] for r in hr]), 0)
# M6: tail, per cell
tl = []
for c in CELLS:
    a = by_start(R, SP["LifeFive"], c, "main", "stretch_p99", 100); b = by_start(R, SP["AsrLong"], c, "main", "stretch_p99", 100)
    tl.append(np.mean(list(a.values())) - np.mean(list(b.values())))
M["TailLower"] = str(sum(x < 0 for x in tl)); M["TailHigher"] = str(sum(x > 0 for x in tl))

# M15: offline optimum
O = json.load(open(os.path.join(OUT, "offline_opt.json")))
ratios = []
gs, os_ = [], []
for c in CELLS:
    g_, o_ = [], []
    for k, v in O.items():
        if k.startswith(c[0] + "_" + c[1]) and v.get("opt_per_h"):
            g_.append(v["greedy_per_h"]); o_.append(v["opt_per_h"])
    ratios.append(np.mean(g_) / np.mean(o_))
M["DpGapMax"] = f(100 * (max(ratios) - 1), 1)
M["DpN"] = str(sum(1 for v in O.values() if v.get("opt_per_h")))
M["DpEq"] = str(sum(1 for v in O.values() if v.get("opt_per_h") and abs(v["greedy_per_h"] - v["opt_per_h"]) < 1e-9))

# M4: Proposition 2 out of sample. Lambda and route delay from control-plane traces at held-out starts (no packets); ns-3 loss at the main starts.
pred, meas, lab = [], [], []
PRD = {}
for name, sp, wait in (("Od", SP["Od"], 0.05), ("SnapSlow", SP["SnapSlow"], 0.5)):
    for c in CELLS:
        t = [r for r in R if r["part"] == "trace" and r["shell"] == c[0] and r["pair"] == c[1] and match(r, sp)]
        if not t:
            continue
        lam = np.mean([r["breaks_per_h"] for r in t]) / 3600; W = np.mean([r["route_ms_mean"] for r in t]) / 1e3
        PRD[(name, c)] = (lam, W, wait)
        m = list(by_start(R, sp, c, "main", "loss", 100).values())
        if m:
            p = 100 * lam * (wait + 1.0 + 0.5 * W)
            pred.append(p); meas.append(np.mean(m)); lab.append(name)
pred, meas = np.array(pred), np.array(meas)
rel = np.abs(meas - pred) / pred
M["OosN"] = str(len(pred)); M["OosErrMed"], M["OosErrMax"] = f(100 * np.median(rel), 1), f(100 * rel.max(), 1)
M["TraceStarts"] = "12"; M["NTrace"] = f(sum(r["part"] == "trace" for r in R), 0)
for c, ct in CT.items():
    if ("Od", c) in PRD:
        M[f"LamOd{ct}"] = f(3600 * PRD[("Od", c)][0], 1)
json.dump(dict(pred=pred.tolist(), meas=meas.tolist(), lab=lab), open(os.path.join(OUT, "stale_check.json"), "w"))
# same-run check over all reactive runs: accounting
rp, rm = [], []
for r in R:
    d = r["cfgd"]
    if r["part"] not in ("main", "tacq") or d.get("policy") not in ("snap", "ondemand") or d.get("bbm"):
        continue
    wait = d["period"] / 2 if d["policy"] == "snap" else 0.05
    x = r["breaks_per_h"] / 3600 * (wait + d.get("tAcq", 0.0) + 0.5 * r["delay_ms_mean"] / 1e3)
    if x > 0:
        rp.append(x); rm.append(r["loss"])
if rp:
    e = np.abs(np.array(rm) - np.array(rp)) / np.array(rp)
    M["CheckN"] = str(sum(r["part"] == "check" for r in R))
    M["SameN"], M["SameErrMed"], M["SameErrNf"] = str(len(rp)), f(100 * np.median(e), 1), f(100 * np.percentile(e, 95), 1)

# M12: acquisition time, two terminals and one terminal (Starlink London-New York)
CELL0 = ("starlink", "lon-nyc")
TAS = ((0.0, "Zero"), (0.2, "Fifth"), (1.0, "One"), (3.0, "Three"))
AT = {}
for ta, tt in TAS:
    for tag in ("SnapSlow", "Od", "ElevtFive", "TrkFive", "LifeFive"):
        sp = SP[tag]
        rr = [r for r in R if r["shell"] == CELL0[0] and r["pair"] == CELL0[1] and match(r, sp) and abs(r["cfgd"].get("tAcq", -1) - ta) < 1e-9
              and ((r["part"] == "tacq") or (r["part"] == "main"))]
        if rr:
            M[f"LossAcq{tt}{tag}"] = f(100 * np.mean([r["loss"] for r in rr]), 2); AT[(ta, tag, 0)] = rr
        rb = [r for r in R if r["part"] == "bbm" and r["shell"] == CELL0[0] and r["pair"] == CELL0[1] and match(r, sp) and abs(r["cfgd"].get("tAcq", -1) - ta) < 1e-9]
        if rb:
            M[f"LossBbm{tt}{tag}"] = f(100 * np.mean([r["loss"] for r in rb]), 2); AT[(ta, tag, 1)] = rb
# one terminal, t_a = 1 s, all cells
BB = [("Od", SP["Od"]), ("ElevtFive", SP["ElevtFive"]), ("TrkFive", SP["TrkFive"]), ("AsrLong", SP["AsrLong"]), ("AsrThree", dict(policy="asr", window=180.0)), ("LifeFive", SP["LifeFive"]), ("LifeTen", SP["LifeTen"])]
for tag, sp in BB:
    pool = [r for r in R if r["part"] == "bbm" and match(r, sp) and abs(r["cfgd"].get("tAcq", -1) - 1.0) < 1e-9]
    if pool:
        s = strat(pool, sp, "loss", 100, "bbm")
        M[f"LossBbm{tag}"], M[f"LossBbm{tag}Ci"] = f(s[0], 2), f(s[1], 2)
        for c, ct in CT.items():
            v = [r["loss"] for r in pool if (r["shell"], r["pair"]) == c]
            if v:
                M[f"LossBbm{tag}{ct}"] = f(100 * np.mean(v), 2)
        # Prop 2(4): loss = acquisitions x t_a, from the two-terminal runs of the same starts
        e = []
        for r in pool:
            a = [q["acq_per_h"] for q in R if q["part"] == "main" and q["shell"] == r["shell"] and q["pair"] == r["pair"] and q["start"] == r["start"] and match(q, sp)]
            if a and a[0] > 0 and tag != "Od":
                e.append(abs(100 * r["loss"] - 100 * a[0] / 3600) / (100 * a[0] / 3600))
        if e:
            M[f"BbmErr{tag}"] = f(100 * np.median(e), 1)
ee = []
ee0 = []
for tag, sp in BB:
    if tag == "Od":
        continue
    for r in [q for q in R if q["part"] == "bbm" and match(q, sp) and abs(q["cfgd"].get("tAcq", -1) - 1.0) < 1e-9]:
        a = [q for q in R if q["part"] == "main" and q["shell"] == r["shell"] and q["pair"] == r["pair"] and q["start"] == r["start"] and match(q, sp)]
        if a and a[0]["acq_per_h"] > 0:
            q = a[0]
            both = q["ingress_per_h"] + q["egress_per_h"] - q["acq_per_h"]
            x = (q["acq_per_h"] * 1.0 + both * q["route_ms_mean"] / 1e3) / 3600
            ee.append(abs(r["loss"] - x) / x)
            ee0.append(abs(r["loss"] - q["acq_per_h"] / 3600) / (q["acq_per_h"] / 3600))
if ee:
    M["BbmN"], M["BbmErrMed"], M["BbmErrNf"], M["BbmErrNoW"] = str(len(ee)), f(100 * np.median(ee), 1), f(100 * np.percentile(ee, 95), 1), f(100 * np.median(ee0), 1)
# ephemeris error (per shell, m3)
for sh, t in (("starlink", "Starlink"), ("oneweb", "Oneweb")):
    for gd, gt in ((0.5, "Half"), (2.0, "Two")):
        for sg, s_ in ((5.0, "Five"), (20.0, "Twenty")):
            rr = [r for r in R if r["part"] == "pred" and r["shell"] == sh and r["cfgd"].get("policy") == "life" and abs(r["cfgd"].get("guard", -1) - gd) < 1e-9 and abs(r["cfgd"].get("predSigmaKm", -1) - sg) < 1e-9]
            if rr:
                M[f"LossG{gt}S{s_}{t}"] = f(100 * np.mean([r["loss"] for r in rr]), 3)
    for sg, s_ in ((5.0, "Five"), (20.0, "Twenty")):
        for pol, pt in (("asr", "Asr"), ("asrk", "Asrk")):
            rr = [r for r in R if r["part"] == "pred" and r["shell"] == sh and r["cfgd"].get("policy") == pol and abs(r["cfgd"].get("predSigmaKm", -1) - sg) < 1e-9]
            if rr:
                M[f"Loss{pt}S{s_}{t}"] = f(100 * np.mean([r["loss"] for r in rr]), 3)
# cross-plane cutoff
for tag, sp in [("Od", SP["Od"]), ("Elev", SP["Elev"]), ("ElevtFive", SP["ElevtFive"]), ("TrkFive", SP["TrkFive"]), ("AsrLong", SP["AsrLong"]),
                ("AsrkLong", SP["AsrkLong"]), ("Alpr", SP["Alpr"]), ("LifeFive", SP["LifeFive"])]:
    rr = [r for r in R if r["part"] == "cross" and match(r, dict(sp, crossLatDeg=70.0))]
    if rr:
        M[f"CrossLoss{tag}"] = f(100 * np.mean([r["loss"] for r in rr]), 2)
        M[f"CrossChg{tag}"] = f(np.mean([r["changes_per_h"] for r in rr]), 1)
        M[f"CrossIslDown{tag}"] = f(np.mean([r["isl_down_events"] for r in rr]), 0)
# real element sets
rr = [r for r in R if r["part"] == "real"]
if rr:
    M["RealRuns"] = str(len(rr)); M["RealStarts"] = str(len({r["start"] for r in rr}))
    M["RealPropagator"] = {r["propagator"] for r in rr}.pop()
    M["RealIslDownMax"] = f(max(r["isl_down_events"] for r in rr), 0)
    M["RealBreaksIsl"] = f(sum(r["breaks_isl_per_h"] for r in rr), 1)
    dp = [r for r in rr if not r["cfgd"].get("noData") and r["cfgd"].get("policy") == "life"]
    M["RealLostLife"] = str(int(sum(r["sent"] - r["received"] for r in dp))); M["RealLifeDataRuns"] = str(len(dp))
    # all real runs of the search lose nothing by forced breaks
    M["RealForcedLife"] = f(max([r["forced_per_h"] for r in rr if r["cfgd"].get("policy") == "life"]), 2)
    od_ = [r for r in rr if match(r, SP["Od"]) and not r["cfgd"].get("noData")]
    e_ = [abs(r["loss"] - r["breaks_per_h"] / 3600 * (0.05 + 1.0 + 0.5 * r["delay_ms_mean"] / 1e3)) / max(r["loss"], 1e-9) for r in od_ if r["loss"] > 0]
    if e_:
        M["RealOdErrMed"] = f(100 * np.median(e_), 1)
# noise
for tag, cfg in (("Thr", dict(trackArc=30.0, trackSamples=3)), ("Lsq", dict(trackArc=60.0, trackSamples=30))):
    for sg, st_ in ((0.01, "Hund"), (0.05, "Five")):
        for ee, et in ((0.05, "E"), (0.10, "T")):
            sp = dict(policy="track", eps=ee, elevNoiseDeg=sg, **cfg)
            for c, ct in CT.items():
                a, b, fr = (by_start(R, sp, c, "noise", k, x) for k, x in (("acq_per_h", 1), ("stretch_mean", 100), ("forced_per_h", 1)))
                if a:
                    M[f"Ns{tag}{st_}{et}{ct}A"], M[f"Ns{tag}{st_}{et}{ct}S"], M[f"Ns{tag}{st_}{et}{ct}F"] = f(np.mean(list(a.values())), 1), f(np.mean(list(b.values())), 1), f(np.mean(list(fr.values())), 2)
            pool = [r for r in R if r["part"] == "noise" and match(r, sp) and r["cfgd"].get("noData")]
            if pool:
                M[f"NsForced{tag}{st_}{et}"] = f(np.mean([r["forced_per_h"] for r in pool]), 2)
dp = [r for r in R if r["part"] == "noise" and not r["cfgd"].get("noData")]
if dp:
    M["NoiseDataRuns"] = str(len(dp)); M["NoiseDataLoss"] = f(100 * np.mean([r["loss"] for r in dp]), 3)
# floor (m9)
for sh, t in TAG.items():
    M[f"ChgMrv{t}Ratio"] = f(np.mean([r["changes_per_h"] for r in R if match(r, SP["Mrv"], "main") and r["shell"] == sh]) / FLOOR[sh], 2) if any(match(r, SP["Mrv"], "main") for r in R) else "?"
below = 0
nn = 0
for r in R:
    if r["part"] in MAINX and r["forced_per_h"] == 0:
        nn += 1; below += r["changes_per_h"] < FLOOR[r["shell"]] * 0.75 * 0 + FLOOR[r["shell"]] - 1e-9
M["BelowFloor"], M["FloorRuns"] = str(below), str(nn)
# merged fallback share (m6) and the price of ISASR
for w, t in ((60.0, "Sixty"), (120.0, "Long"), (180.0, "Three"), (300.0, "Vlong")):
    v = [r["fallback_share"] for r in R if r["part"] in ("main", "front") and match(r, dict(policy="asr", window=w)) and r["fallback_share"] == r["fallback_share"]]
    if v:
        M[f"Fallback{t}"] = f(100 * np.mean(v), 0)
M["IsasrDijRatio"] = f(float(M["DijIsasr"].replace(",", "")) / float(M["DijLifeFive"].replace(",", "")), 0)
# cost
for tag in ("LifeFive",):
    rr = [r for r in R if match(r, SP[tag], "main") and r["search_cpu_ms"] == r["search_cpu_ms"]]
    M["SearchMsMed"] = f(np.median([r["search_cpu_ms"] for r in rr]), 1)
    M["SearchMsMax"] = f(max(r["search_cpu_ms"] for r in rr), 1)
    M["DijPerSearch"] = f(np.mean([r["dijkstra_per_h"] / max(r["computations_per_h"], 1e-9) for r in rr]), 1)
M["CpuPerRun"] = f(np.median([r["cpu_s"] for r in R if not r["cfgd"].get("noData")]), 0)

with open(os.path.join(HERE, "..", "macros_p7.tex"), "w") as fh:
    fh.write("% generated by sim/make_macros.py; do not edit\n")
    for k, v in sorted(M.items()):
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print(len(M), "macros")
for k in sorted(M):
    if k.startswith(("R", "Rng", "Oos", "Dp", "Front", "Rival", "Trk")) and not k.startswith(("Real", "Rcap")):
        print(k, M[k], end="; ")
print()
