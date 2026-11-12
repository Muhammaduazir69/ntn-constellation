#!/usr/bin/env python3
"""Macros for the TAES correspondence (Nov 2026). Every number in the text comes from here.
Data: out/p7_v8.jsonl (binary r2) for configurations the r3 change cannot reach, and out/p7_v9.jsonl (binary r3) for ISASR,
the real-shell elevation sweep and the held-out tests; r3 reproduced 148 stratified r2 rows bit for bit (repro rows)."""
import json, os
import numpy as np
from scipy import stats as st
import p7util
from p7util import CELLS, OUT, HERE, strat, by_start
import frontier_stats as F
import run_p7 as R7

M = {}


def f(x, d=1):
    return f"{x:,.{d}f}" if abs(x) >= 1e4 else f"{x:.{d}f}"


def load(path, keep=None):
    out = []
    for l in open(os.path.join(OUT, path)):
        r = json.loads(l)
        if "error" in r:
            continue
        r["cfgd"] = dict(r["cfg"])
        if keep is None or keep(r):
            out.append(r)
    return out


_m = p7util.match
def match2(r, spec, part=None):
    d = r["cfgd"]
    if "isasrIslAct" in d and "isasrIslAct" not in spec:
        return False
    if "isasrIslAct" in spec and d.get("isasrIslAct") != spec["isasrIslAct"]:
        return False
    return _m(r, {k: v for k, v in spec.items() if k != "isasrIslAct"}, part)
p7util.match = match2; F.match = match2

v8 = load("p7_v8.jsonl", lambda r: r["cfgd"].get("policy") != "isasr")
v9 = load("p7_v9.jsonl")
# ---- one binary: reproduction
idx = {(r["part"], r["shell"], r["pair"], r["start"], json.dumps(r["cfg"])): r for r in load("p7_v8.jsonl")}
KEYS = ("changes_per_h", "acq_per_h", "stretch_mean", "stretch_p99", "forced_per_h", "dijkstra_per_h", "sent", "received")
same = tot = 0
for r in v9:
    if r["part"].startswith("repro:"):
        o = idx.get((r["part"][6:], r["shell"], r["pair"], r["start"], json.dumps(r["cfg"])))
        if o:
            tot += 1; same += all(o.get(k) == r.get(k) for k in KEYS)
M["ReproN"], M["ReproSame"] = str(tot), str(same)
M["ShaOld"], M["ShaNew"] = {r["bin_sha"][:12] for r in v8}.pop(), {r["bin_sha"][:12] for r in v9}.pop()
S6 = set(R7.S6); SF = set(R7.SF)
new = []
for r in v9:
    if r["part"].startswith("repro:"):
        continue
    p = r["part"]
    if p in ("isasr", "relev"):
        r["part"] = "real" if r["shell"].startswith("real_") else ("fresh" if r["start"] in SF and r["shell"] in ("starlink", "kuiper") else ("extra" if r["shell"] == "oneweb" and r["start"] not in S6 else "main"))
    new.append(r)
R = v8 + new
M["NRuns"] = f(len([r for r in R if r["part"] != "fresh"]), 0)
M["NRunsNew"] = f(len(new), 0)
M["NData"] = f(sum(not r["cfgd"].get("noData") for r in R), 0)

# ---- link law and ground breaks
L = json.load(open(os.path.join(OUT, "law.json"))); G = json.load(open(os.path.join(OUT, "grid_allF.json")))
for s, t in (("starlink", "Starlink"), ("kuiper", "Kuiper"), ("oneweb", "Oneweb")):
    M[f"GridDmax{t}"], M[f"GridDmin{t}"], M[f"GridAllF{t}"] = f(L[s]["grid_dmax_km"], 0), f(L[s]["grid_dmin_km"], 0), f(G[s], 0)
    M[f"Rcap{t}"] = f(L[s]["rcap_km"], 0)
    M[f"AfEv{t}"], M[f"AfEpochs{t}"] = f(L[s]["af_events"], 0), f(L[s]["af_epochs"], 0)
PC = json.load(open(os.path.join(OUT, "pair_check.json"))); M["PairN"], M["PairErrNm"] = f(PC["pairs"], 0), f(PC["max_err_m"] * 1e9, 0)
pk = [r for r in R if r["part"] == "main" and not r["cfgd"].get("noData")]
M["NPacketMain"] = f(len(pk), 0); M["IslDownPacket"] = f(sum(r["isl_down_events"] for r in pk), 0)
M["IslDownAllMain"] = f(sum(r["isl_down_events"] for r in R if r["part"] in ("main", "extra", "front", "trig")), 0)
M["NAllMain"] = f(sum(1 for r in R if r["part"] in ("main", "extra", "front", "trig")), 0)
rr = [r for r in R if r["part"] == "real"]
M["RealIslDownMax"] = f(max(r["isl_down_events"] for r in rr), 0); M["RealBreaksIsl"] = f(sum(r["breaks_isl_per_h"] for r in rr), 1)
hb = [r for r in new if r["part"] in ("hbench", "hpacket", "htrace")]
M["HeldIslDown"] = f(sum(r["isl_down_events"] for r in hb if not r["shell"].startswith("real")), 0)
M["HeldBreaksIsl"] = f(sum(r["breaks_isl_per_h"] for r in hb), 1)
allc = [r for r in R if r["part"] in ("main", "extra") and r["changes_per_h"] > 0]
M["ShareIslOnly"] = f(100 * sum(r["isl_only_per_h"] for r in allc if r["cfgd"]["policy"] != "isasr") / sum(r["changes_per_h"] for r in allc if r["cfgd"]["policy"] != "isasr"), 1)
od = [r for r in R if match2(r, dict(policy="ondemand"), "main")]
M["BreaksInOd"], M["BreaksOutOd"] = f(np.mean([r["breaks_in_per_h"] for r in od]), 1), f(np.mean([r["breaks_out_per_h"] for r in od]), 1)
for tag, sp in (("Od", dict(policy="ondemand")), ("Snap", dict(policy="snap", period=1.0))):
    s_ = strat(R, sp, "loss", 100, "main"); M[f"Loss{tag}"] = f(s_[0], 2)
    s_ = strat(R, sp, "changes_per_h", 1, "main"); M[f"Chg{tag}"] = f(s_[0], 1)

# ---- loss predictor
def pred_cell(trace, wait):
    lam = np.mean([r["breaks_per_h"] for r in trace]) / 3600; W = np.mean([r["route_ms_mean"] for r in trace]) / 1e3
    return 100 * lam * (wait + 1.0 + 0.5 * W)
PP = []
for name, sp, wait in (("Od", dict(policy="ondemand"), 0.05), ("Snap", dict(policy="snap", period=1.0), 0.5)):
    for c in CELLS:
        t = [r for r in R if r["part"] == "trace" and (r["shell"], r["pair"]) == c and match2(r, sp)]
        m = [100 * r["loss"] for r in R if r["part"] == "main" and (r["shell"], r["pair"]) == c and match2(r, sp)]
        if t and m:
            PP.append(dict(kind="design", c=c, pol=name, pred=pred_cell(t, wait), meas=float(np.mean(m))))
    for c in [(s, p) for s, p in __import__("run_p7_r3").HELD]:
        t = [r for r in new if r["part"] == "htrace" and (r["shell"], r["pair"]) == c and match2(r, sp)]
        m = [100 * r["loss"] for r in new if r["part"] == "hpacket" and (r["shell"], r["pair"]) == c and match2(r, sp)]
        if t and m:
            PP.append(dict(kind="held", c=c, pol=name, pred=pred_cell(t, wait), meas=float(np.mean(m)), measlo=min(m), meashi=max(m)))
for kind, tag in (("design", "Des"), ("held", "Held")):
    e = [abs(p["meas"] - p["pred"]) / p["pred"] for p in PP if p["kind"] == kind]
    if e:
        M[f"Pred{tag}N"], M[f"Pred{tag}Med"], M[f"Pred{tag}Max"] = str(len(e)), f(100 * np.median(e), 1), f(100 * max(e), 1)
json.dump(PP, open(os.path.join(OUT, "pred_check.json"), "w"), indent=0, default=str)
# OneWeb: a cell mean is not a run; per-run prediction band from 40 trace starts
for name, sp, wait in (("Od", dict(policy="ondemand"), 0.05), ("Snap", dict(policy="snap", period=1.0), 0.5)):
    inside = tot = 0; widths = []
    for pr in ("lon-nyc", "syd-lax"):
        t = [r for r in new if r["part"] == "owtrace" and r["pair"] == pr and match2(r, sp)]
        if not t:
            continue
        per = [100 * r["breaks_per_h"] / 3600 * (wait + 1.0 + 0.5 * r["route_ms_mean"] / 1e3) for r in t]
        lo, hi = np.percentile(per, [2.5, 97.5]); widths.append((lo, hi))
        for r in R:
            if r["part"] == "main" and r["shell"] == "oneweb" and r["pair"] == pr and match2(r, sp):
                tot += 1; inside += lo <= 100 * r["loss"] <= hi
        M[f"OwBand{name}{'L' if pr == 'lon-nyc' else 'S'}Lo"], M[f"OwBand{name}{'L' if pr == 'lon-nyc' else 'S'}Hi"] = f(lo, 2), f(hi, 2)
        cm = np.mean(per); M[f"OwMean{name}{'L' if pr == 'lon-nyc' else 'S'}"] = f(cm, 2)
    if tot:
        M[f"OwInside{name}"], M[f"OwTot{name}"] = str(inside), str(tot)
M["OwTraceStarts"] = str(len(R7.np.random.default_rng(20261112).uniform(0, 1, 40)))
# one-terminal consistency
ee = []
for r in [q for q in R if q["part"] == "bbm" and abs(q["cfgd"].get("tAcq", -1) - 1.0) < 1e-9 and q["cfgd"]["policy"] not in ("ondemand", "snap")]:
    a = [q for q in R if q["part"] == "main" and q["shell"] == r["shell"] and q["pair"] == r["pair"] and q["start"] == r["start"] and match2(q, {k: v for k, v in r["cfgd"].items() if k not in ("bbm", "tAcq")})]
    if a and a[0]["acq_per_h"] > 0:
        q = a[0]; both = q["ingress_per_h"] + q["egress_per_h"] - q["acq_per_h"]
        x = (q["acq_per_h"] + both * q["route_ms_mean"] / 1e3) / 3600
        ee.append(abs(r["loss"] - x) / x)
M["BbmN"], M["BbmErrMed"], M["BbmErrNf"] = str(len(ee)), f(100 * np.median(ee), 1), f(100 * np.percentile(ee, 95), 1)
for tag, sp in (("Life", dict(policy="life", eps=0.05)), ("Asr", dict(policy="asr", window=120.0)), ("Od", dict(policy="ondemand"))):
    pool = [r for r in R if r["part"] == "bbm" and match2(r, sp) and abs(r["cfgd"].get("tAcq", -1) - 1.0) < 1e-9]
    s_ = strat(pool, sp, "loss", 100, "bbm"); M[f"LossBbm{tag}"] = f(s_[0], 2)

# ---- benchmark: best loss-free setting per cell, acquisitions at equal delay, jackknife CI
C = F.Cells(R)
HC = F.Cells([dict(r, part="main") for r in new if r["part"] == "hbench"])
G_ = (0.1, 0.3, 1.0, 3.0, 1000.0)
FAM = {"Elev": F.FAM["Elev"], "Merged": F.FAM["Merged"], "IsasrA": [dict(policy="isasr", gamma=g, isasrIslAct=1) for g in G_],
       "IsasrV": [dict(policy="isasr", gamma=g) for g in (0.3, 1.0, 3.0)], "Alpr": F.FAM["Alpr"], "Trk": F.FAM["Trk"]}
HFAM = {"Elev": [dict(policy="elev", trigger=24.0), F.E(0.02, 24.0, 1.0), F.E(0.05, 24.0), F.E(0.10, 24.0)],
        "Merged": [dict(policy="asr", window=w) for w in (60.0, 120.0, 180.0)],
        "IsasrA": [dict(policy="isasr", gamma=g, isasrIslAct=1) for g in (1.0, 1000.0)], "IsasrV": [dict(policy="isasr", gamma=1.0)],
        "Trk": [dict(policy="track", eps=e) for e in (0.05, 0.10)]}
HLIFE = [dict(policy="life", eps=e) for e in (0.0, 0.02, 0.05, 0.10, 0.15)]
CT = dict(zip(CELLS + F.REALC, ("StL", "StS", "KuL", "KuS", "OnL", "OnS", "RsL", "RsS", "RoL", "RoS")))
HELDC = __import__("run_p7_r3").HELD
HT = dict(zip(HELDC, ("HkM", "HkS", "HsM", "HsS", "HrL", "HrS")))
res = {}
def best(Cx, c, sets, life, metric, axis):
    b = None
    for sp in sets:
        o = F.compare(Cx, c, sp, life, metric, axis)
        if o is None or o["forced"] >= 0.1:
            continue
        if o["cl"] == "dom":
            return dict(o, sp=sp)
        if o["cl"] == "on" and (b is None or b["cl"] != "on" or o["ratio"] < b["ratio"]):
            b = dict(o, sp=sp)
        elif b is None:
            b = dict(o, sp=sp)
    return b
TESTABLE = {"oneweb", "real_starlink", "real_oneweb", "real_starlink16"}  # starts vary; lattice starts repeat with the plane spacing


def cellfmt(o, ci=True):
    if o is None:
        return "?"
    if o["cl"] != "on":
        return {"dom": "D", "beyond": "$>$", "below": "$<$"}[o["cl"]]
    s = f(o["ratio"], 2)
    if ci and o.get("testable") and o["lo"] == o["lo"] and (o["lo"] > 1 or o["hi"] < 1):
        s = r"\textbf{" + s + "}"
    return s
for axis, at in (("stretch_mean", "Mn"), ("stretch_p99", "Pn")):
    for fam, sets in FAM.items():
        for c, ct in CT.items():
            o = best(C, c, sets, F.LIFE, "acq_per_h", axis); res[(fam, at, c)] = o
            if o:
                o["testable"] = c[0] in TESTABLE
            M[f"B{fam}{at}{ct}"] = cellfmt(o)
        for c, ct in HT.items():
            if fam in HFAM:
                o = best(HC, c, HFAM[fam], HLIFE, "acq_per_h", axis); res[(fam, at, c)] = o
                if o:
                    o["testable"] = c[0] in TESTABLE
                M[f"B{fam}{at}{ct}"] = cellfmt(o)
            else:
                M[f"B{fam}{at}{ct}"] = "--"
def rng(fam, at, cells, tag):
    v = [res[(fam, at, c)]["ratio"] for c in cells if res.get((fam, at, c)) and res[(fam, at, c)]["cl"] == "on"]
    if v:
        M[f"R{fam}{at}{tag}Lo"], M[f"R{fam}{at}{tag}Hi"] = f(min(v), 2), f(max(v), 2)
    M[f"D{fam}{at}{tag}"] = str(sum(1 for c in cells if res.get((fam, at, c)) and res[(fam, at, c)]["cl"] == "dom"))
for fam in FAM:
    for at in ("Mn", "Pn"):
        rng(fam, at, CELLS, "Lat"); rng(fam, at, F.REALC, "Real"); rng(fam, at, HELDC, "Held")
# gamma optimum per real cell (grid edge, m4)
for c in F.REALC:
    best_g = None
    for g in G_:
        o = F.compare(C, c, dict(policy="isasr", gamma=g, isasrIslAct=1), F.LIFE, "acq_per_h", "stretch_mean")
        if o and o["cl"] == "on" and (best_g is None or o["ratio"] < best_g[1]):
            best_g = (g, o["ratio"])
    if best_g:
        M[f"IsasrBestGamma{CT[c]}"] = f(best_g[0], 1)
# cross-plane latitude cutoff (OneWeb, 70 deg), r2 rows
cr = [r for r in R if r["part"] == "cross"]
M["CrossIslDown"] = f(np.mean([r["isl_down_events"] for r in cr]), 0)
for tag, sp in (("Od", dict(policy="ondemand")), ("Elevt", dict(policy="elevb", eps=0.05, trigger=24.0)), ("Life", dict(policy="life", eps=0.05)),
                ("Asr", dict(policy="asr", window=120.0)), ("Alpr", dict(policy="alpr"))):
    v = [r for r in cr if match2(r, dict(sp, crossLatDeg=70.0))]
    M[f"CrossBrk{tag}"] = f(max(r["breaks_isl_per_h"] for r in v), 0); M[f"CrossLoss{tag}"] = f(100 * np.mean([r["loss"] for r in v]), 2)
M["CrossBrkMax"] = f(max(r["breaks_isl_per_h"] for r in cr), 0)
# ISASR price
for tag, sp in (("IsasrA", dict(policy="isasr", gamma=1000.0, isasrIslAct=1)), ("IsasrAOne", dict(policy="isasr", gamma=1.0, isasrIslAct=1)), ("IsasrV", dict(policy="isasr", gamma=1.0))):
    a = [r for r in R if match2(r, sp) and r["part"] == "main"]
    M[f"Dij{tag}"] = f(np.mean([r["dijkstra_per_h"] for r in a]), 0)
    M[f"IslOnly{tag}"] = f(np.mean([r["isl_only_per_h"] for r in a]), 1)
    M[f"Chg{tag}"] = f(np.mean([r["changes_per_h"] for r in a]), 1); M[f"Acq{tag}"] = f(np.mean([r["acq_per_h"] for r in a]), 1)
l5 = [r for r in R if match2(r, dict(policy="life", eps=0.05)) and r["part"] == "main"]
M["DijLife"] = f(np.mean([r["dijkstra_per_h"] for r in l5]), 0); M["SearchMs"] = f(np.median([r["search_cpu_ms"] for r in l5]), 1)
M["DijRatio"] = f(float(M["DijIsasrA"].replace(",", "")) / float(M["DijLife"].replace(",", "")), 0)
M["DataLossIsasr"] = f(100 * np.mean([r["loss"] for r in R if match2(r, dict(policy="isasr", gamma=1.0, isasrIslAct=1)) and not r["cfgd"].get("noData")]), 2)
# route changes of the cost-free variant
cv = []
for c in CELLS:
    o = best(C, c, FAM["IsasrV"], F.LIFE, "changes_per_h", "stretch_mean")
    if o and o["cl"] == "on":
        cv.append(o["ratio"])
M["IsasrVChgLo"], M["IsasrVChgHi"] = f(min(cv), 2), f(max(cv), 2)
# counts
M["NStarts"] = "6"; M["NStartsOneweb"] = "20"; M["RealStarts"] = "12"; M["HeldStarts"] = "6"; M["HeldTraceStarts"] = "12"
M["NElevSettings"] = str(len(F.FAM["Elev"]))
# noise
for tag, cfg in (("Lsq", dict(trackArc=60.0, trackSamples=30)), ("Thr", dict(trackArc=30.0, trackSamples=3))):
    v = []
    for e in (0.05, 0.10):
        for c in CELLS + F.REALC:
            o = F.compare(C, c, dict(policy="track", eps=e, elevNoiseDeg=0.05, **cfg), F.LIFE, "acq_per_h", "stretch_mean")
            if o and o["cl"] == "on":
                v.append(o["ratio"])
    if v:
        M[f"Noise{tag}Med"], M[f"Noise{tag}Lo"], M[f"Noise{tag}Hi"] = f(float(np.median(v)), 2), f(min(v), 2), f(max(v), 2)
with open(os.path.join(HERE, "..", "macros_p7.tex"), "w") as fh:
    fh.write("% generated by sim/make_macros_corr.py; do not edit\n")
    for k, v in sorted(M.items()):
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print(len(M), "macros")
for k in sorted(M):
    if k.startswith(("R", "D", "Pred", "Ow", "Repro", "Held", "Bbm", "Dij", "Noise", "Isasr")) and not k.startswith("Rcap"):
        print(k, M[k], end="; ")
print()
