#!/usr/bin/env python3
"""Round-3 stage 1: rival/search ratios at equal delay, best loss-free setting per cell, with jackknife 95 % intervals.
Rows: r2 binary for the search and untouched rivals (p7_v8), r3 binary for ISASR and the real-shell elevation sweep (p7_v9);
r3 reproduces r2 on untouched configurations (repro rows). Writes out/stage1_r3.json and prints a table."""
import json, os, sys
import numpy as np
from p7util import match, CELLS, OUT
import frontier_stats as F

def rows(path, keep=None):
    out = []
    for l in open(os.path.join(OUT, path)):
        r = json.loads(l)
        if "error" in r:
            continue
        r["cfgd"] = dict(r["cfg"])
        if keep is None or keep(r):
            out.append(r)
    return out

v8 = rows("p7_v8.jsonl", lambda r: r["cfgd"].get("policy") != "isasr")
v9 = rows("p7_v9.jsonl")
rep = [r for r in v9 if r["part"].startswith("repro:")]
# repro check
idx = {(r["part"], r["shell"], r["pair"], r["start"], json.dumps(r["cfg"])): r for r in rows("p7_v8.jsonl")}
KEYS = ("changes_per_h", "acq_per_h", "stretch_mean", "stretch_p99", "forced_per_h", "dijkstra_per_h", "sent", "received")
same = diff = 0
bad = []
for r in rep:
    o = idx.get((r["part"][6:], r["shell"], r["pair"], r["start"], json.dumps(r["cfg"])))
    if o is None:
        continue
    ok = all((o.get(k) == r.get(k)) or (isinstance(o.get(k), float) and abs(o[k] - r[k]) <= 1e-9 * max(1, abs(o[k]))) for k in KEYS)
    same += ok; diff += (not ok)
    if not ok:
        bad.append((r["part"], r["shell"], r["pair"], r["start"], r["cfg"]))
print("repro identical", same, "differ", diff, bad[:5])
new = [r for r in v9 if not r["part"].startswith("repro:")]
for r in new:
    r["part"] = {"isasr": "main", "relev": "real"}.get(r["part"], r["part"])
    if r["shell"].startswith("real_"):
        r["part"] = "real"
    elif r["shell"] == "oneweb" and r["start"] not in __import__("run_p7").S6:
        r["part"] = "extra"
fresh = [r for r in new if r["shell"] in ("starlink", "kuiper") and r["start"] in __import__("run_p7").SF]
for r in fresh:
    r["part"] = "fresh"
R = v8 + new
C = F.Cells(R)
G = (0.3, 1.0, 3.0, 1000.0)
FAM = {
    "Elevation rank": F.FAM["Elev"],
    "Merged snapshot": F.FAM["Merged"],
    "ISASR, Algorithm 3": [dict(policy="isasr", gamma=g, isasrIslAct=1) for g in G],
    "ISASR, ISL cost-free (r2)": [dict(policy="isasr", gamma=g) for g in G[:3]],
    "ALPR": F.FAM["Alpr"],
    "Tracked search, noise-free": F.FAM["Trk"],
}
# match() rejects unknown knobs; isasrIslAct must be admitted
import p7util
_m = p7util.match
def match2(r, spec, part=None):
    d = r["cfgd"]
    if "isasrIslAct" in d and "isasrIslAct" not in spec:
        return False
    s2 = {k: v for k, v in spec.items() if k != "isasrIslAct"}
    if "isasrIslAct" in spec and d.get("isasrIslAct") != spec["isasrIslAct"]:
        return False
    return _m(r, s2, part)
F.match = match2
p7util.match = match2
REALC = F.REALC
allc = CELLS + REALC
out = {}
for metric in ("acq_per_h", "changes_per_h"):
    for axis in ("stretch_mean", "stretch_p99"):
        for fam, sets in FAM.items():
            for c in allc:
                best = None
                for sp in sets:
                    o = F.compare(C, c, sp, F.LIFE, metric, axis)
                    if o is None or o["forced"] >= 0.1:
                        continue
                    if o["cl"] == "dom":
                        best = dict(o, sp=sp); break
                    if o["cl"] == "on" and (best is None or best["cl"] != "on" or o["ratio"] < best["ratio"]):
                        best = dict(o, sp=sp)
                    elif best is None:
                        best = dict(o, sp=sp)
                out[(metric, axis, fam, c)] = best
def fmt(o):
    if o is None:
        return "    --    "
    if o["cl"] == "on":
        return f"{o['ratio']:.2f}[{o['lo']:.2f},{o['hi']:.2f}]" if o["lo"] == o["lo"] else f"{o['ratio']:.2f}[  n/a  ]"
    return {"dom": "D", "beyond": ">", "below": "<"}[o["cl"]].center(10)
names = ["St-L", "St-S", "Ku-L", "Ku-S", "On-L", "On-S", "rSt-L", "rSt-S", "rOw-L", "rOw-S"]
for metric in ("acq_per_h", "changes_per_h"):
    for axis in ("stretch_mean", "stretch_p99"):
        print(f"\n== {metric} at equal {axis}; best loss-free setting per cell; jackknife 95 % CI")
        print(" " * 28 + "".join(n.center(17) for n in names))
        for fam in FAM:
            print(f"{fam:28s}" + "".join(fmt(out[(metric, axis, fam, c)]).center(17) for c in allc))
json.dump({"|".join([m, a, f, c[0], c[1]]): ({k: (v if not isinstance(v, (np.floating,)) else float(v)) for k, v in o.items()} if o else None)
           for (m, a, f, c), o in out.items()}, open(os.path.join(OUT, "stage1_r3.json"), "w"), indent=0, default=str)
# ISASR Algorithm 3 cost: route changes, ISL-only, Dijkstra, against the search eps 0.05
print()
for g in G:
    sp = dict(policy="isasr", gamma=g, isasrIslAct=1)
    a = [r for r in R if match2(r, sp) and r["shell"] in [c[0] for c in CELLS]]
    if a:
        print(f"ISASR Alg.3 gamma {g}: changes/h {np.mean([r['changes_per_h'] for r in a]):.1f}, ISL-only {np.mean([r['isl_only_per_h'] for r in a]):.1f}, acq {np.mean([r['acq_per_h'] for r in a]):.1f}, stretch {100*np.mean([r['stretch_mean'] for r in a]):.1f} %, Dijkstra/h {np.mean([r['dijkstra_per_h'] for r in a]):.0f}, cpu_s {np.mean([r['cpu_s'] for r in a]):.1f}, forced {np.mean([r['forced_per_h'] for r in a]):.2f}")
l = [r for r in R if match2(r, dict(policy="life", eps=0.05)) and r["part"] == "main"]
print(f"search eps 0.05: changes/h {np.mean([r['changes_per_h'] for r in l]):.1f}, acq {np.mean([r['acq_per_h'] for r in l]):.1f}, Dijkstra/h {np.mean([r['dijkstra_per_h'] for r in l]):.0f}, cpu_s {np.mean([r['cpu_s'] for r in l]):.1f}, search ms {np.median([r['search_cpu_ms'] for r in l]):.1f}")
dp = [r for r in R if match2(r, dict(policy="isasr", gamma=1.0, isasrIslAct=1)) and not r["cfgd"].get("noData")]
if dp:
    print("ISASR Alg.3 gamma 1 data-plane runs", len(dp), "loss %", 100 * np.mean([r["loss"] for r in dp]))
