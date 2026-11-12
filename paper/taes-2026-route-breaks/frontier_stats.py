"""Comparison of rivals with the lifetime search at equal delay, on cell means over the same starts.
Rule (M2): every loss-free setting of a family enters. Metric: ground-link acquisitions per hour (route changes given as well).
Delay axes (M3): mean stretch and 99th percentile of stretch. Rivals are placed on the search's selection-budget frontier and on the enforced
frontier. Off-frontier rival points are classified; none is dropped. Inference (M5): start-level tests only where starts are independent in
effect (OneWeb lattice, real shells); Starlink and Kuiper lattices repeat with the plane spacing, so they are descriptive."""
import json, os
import numpy as np
from scipy import stats as st
from p7util import match, CELLS, OUT, holm

E = lambda e, t=None, g=None: dict(policy="elevb", eps=e, **({"trigger": t} if t else {}), **({"floorGrid": g} if g else {}))
LIFE = [dict(policy="life", eps=e) for e in (0., .01, .02, .03, .05, .07, .10, .15, .20)]
HOLD = [dict(policy="life", eps=e, hold=1) for e in (.02, .03, .05, .07, .10, .15)]
FAM = {
    "Elev": [E(e, t) for t in (23.3, 23.5, 24.0) for e in (.02, .05, .10)] + [E(.02), E(.05), E(.10), E(.05, 27.0), E(.05, 30.0)] +
            [E(e, t, 1.0) for t in (23.5, 24.0, 25.0) for e in (.02, .05, .10)] + [dict(policy="elev"), dict(policy="elev", trigger=24.0)],
    "Merged": [dict(policy="asr", window=w) for w in (30., 60., 120., 180.)] + [dict(policy="asrk", window=w) for w in (60., 120., 180.)],
    "Isasr": [dict(policy="isasr", gamma=g) for g in (.3, 1., 3.)],
    "Alpr": [dict(policy="alpr")],
    "Trk": [dict(policy="track", eps=e) for e in (.02, .05, .07, .10, .15)] + [dict(policy="track", eps=.05, trigger=24.0), dict(policy="track", eps=.10, trigger=24.0)],
}
PARTS = ("main", "front", "extra", "trig", "real")
NS = {"Ns" + a + b: dict(policy="track", elevNoiseDeg=s_, **c) for a, c in (("Thr", dict(trackArc=30.0, trackSamples=3)), ("Lsq", dict(trackArc=60.0, trackSamples=30))) for b, s_ in (("Hund", .01), ("Five", .05))}
REALC = [("real_starlink", "lon-nyc"), ("real_starlink", "syd-lax"), ("real_oneweb", "lon-nyc"), ("real_oneweb", "syd-lax")]
REALT = dict(zip(REALC, ("RsL", "RsS", "RoL", "RoS")))


class Cells:
    def __init__(self, R):
        self.idx = {}
        self.memo = {}
        for r in R:
            if r["part"] in PARTS or r["part"] == "noise":
                self.idx.setdefault((r["shell"], r["pair"]), {}).setdefault(r["start"], []).append(r)

    def row(self, c, s, sp):
        k = (c, s, json.dumps(sp, sort_keys=True))
        if k in self.memo:
            return self.memo[k]
        v = None
        for r in self.idx.get(c, {}).get(s, []):
            if match(r, sp):
                v = r
                break
        self.memo[k] = v
        return v

    def starts(self, c, sp):
        return sorted(s for s in self.idx.get(c, {}) if self.row(c, s, sp) is not None)


def fr(C, c, F, ss, metric, axis):
    pts = []
    for sp in F:
        rr = [C.row(c, s, sp) for s in ss]
        if all(rr):
            pts.append((np.mean([100 * r[axis] for r in rr]), np.mean([r[metric] for r in rr])))
    return sorted(pts)


def place(pts, x, y):
    """Ratio of a rival point (x delay, y rate) to the frontier; class: 'on', 'beyond' (more delay than any budget), 'below' (less delay),
    'dom' (less delay and no more rate than the frontier's lowest-delay point)."""
    xs, ys = zip(*pts)
    if x > xs[-1]:
        return float("nan"), "beyond"
    if x < xs[0]:
        return float("nan"), "dom" if y <= ys[0] else "below"
    return y / float(np.interp(x, xs, ys)), "on"


def compare(C, c, sp, F, metric, axis):
    ss = [s for s in C.starts(c, sp) if all(C.row(c, s, q) for q in F[:1])]
    # a start enters when the frontier budgets used are all there; frontier uses the budgets present at every one of these starts
    use = [q for q in F if all(C.row(c, s, q) for s in ss)]
    if len(ss) < 2 or len(use) < 2:
        return None
    rr = [C.row(c, s, sp) for s in ss]
    x, y = np.mean([100 * r[axis] for r in rr]), np.mean([r[metric] for r in rr])
    pts = fr(C, c, use, ss, metric, axis)
    rt, cl = place(pts, x, y)
    out = dict(n=len(ss), ratio=rt, cl=cl, x=x, y=y, forced=float(np.mean([r["forced_per_h"] for r in rr])), lo=np.nan, hi=np.nan, p=np.nan)
    if cl == "on":
        jk = []
        for q in range(len(ss)):
            sub = ss[:q] + ss[q + 1:]
            rs = [C.row(c, s, sp) for s in sub]
            v, cc = place(fr(C, c, use, sub, metric, axis), np.mean([100 * r[axis] for r in rs]), np.mean([r[metric] for r in rs]))
            if cc == "on":
                jk.append(v)
        if len(jk) == len(ss):
            ll = np.log(jk); n = len(ss)
            se = np.sqrt((n - 1) / n * np.sum((ll - ll.mean()) ** 2))
            tq = st.t.ppf(.975, n - 1)
            out.update(lo=np.exp(np.log(rt) - tq * se), hi=np.exp(np.log(rt) + tq * se),
                       p=(2 * st.t.sf(abs(np.log(rt)) / se, n - 1) if se > 0 else 1.0))
    return out


def frontier_macros(R, M, f, SP, CT):
    C = Cells(R)
    allc = [(c, CT[c]) for c in CELLS] + [(c, REALT[c]) for c in REALC]
    res = {}
    for fam, sets in FAM.items():
        for k, sp in enumerate(sets):
            for c, ct in allc:
                for fname, F in (("sel", LIFE), ("enf", HOLD)):
                    for axis, an in (("stretch_mean", "mean"), ("stretch_p99", "p99")):
                        o = compare(C, c, sp, F, "acq_per_h", axis)
                        if o is not None:
                            o.update(fam=fam, k=k, sp=sp)
                            res[(fam, k, c, fname, an)] = o
    LOSSFREE = 0.1
    testable = lambda c: c[0] in ("oneweb", "real_starlink", "real_oneweb")
    summary = {}
    for fname in ("sel", "enf"):
        for an in ("mean", "p99"):
            for fam in FAM:
                # lattices
                lat = [v for (fa, k, c, fn, a), v in res.items() if fa == fam and fn == fname and a == an and c in CELLS and v["forced"] < LOSSFREE]
                on = [v["ratio"] for v in lat if v["cl"] == "on"]
                best = []
                beat = 0
                for c in CELLS:
                    cand = [v for (fa, k, cc, fn, a), v in res.items() if fa == fam and fn == fname and a == an and cc == c and v["forced"] < LOSSFREE]
                    o2 = [v["ratio"] for v in cand if v["cl"] == "on"]
                    if any(v["cl"] == "dom" for v in cand):
                        best.append(0.0); beat += 1
                    elif o2:
                        best.append(min(o2))
                t = f"{fam}{'Sel' if fname == 'sel' else 'Enf'}{'Mean' if an == 'mean' else 'Pn'}"
                if on:
                    M[f"Med{t}"], M[f"Lo{t}"], M[f"Hi{t}"] = f(float(np.median(on)), 2), f(min(on), 2), f(max(on), 2)
                    M[f"N{t}"] = str(len(on)); M[f"Off{t}"] = str(len(lat) - len(on)); M[f"Dom{t}"] = str(sum(1 for v in lat if v["cl"] == "dom"))
                pos = [b for b in best if b > 0]
                if best:
                    M[f"Best{t}"] = f(float(np.median(best)), 2)
                    M[f"BestLo{t}"], M[f"BestHi{t}"] = f(min(pos), 2) if pos else "--", f(max(pos), 2) if pos else "--"
                    M[f"BestBelow{t}"] = str(sum(1 for b in best if b < 1.0))
                    M[f"BestDom{t}"] = str(beat)
                # best per cell, table entries
                for c, ct in allc:
                    cand = [v for (fa, k, cc, fn, a), v in res.items() if fa == fam and fn == fname and a == an and cc == c and v["forced"] < LOSSFREE]
                    if not cand:
                        M[f"B{t}{ct}"] = "?"
                        continue
                    o2 = [v for v in cand if v["cl"] == "on"]
                    if any(v["cl"] == "dom" for v in cand):
                        M[f"B{t}{ct}"] = "D"
                    elif o2:
                        v = min(o2, key=lambda z: z["ratio"])
                        M[f"B{t}{ct}"] = f(v["ratio"], 2)
                    else:
                        M[f"B{t}{ct}"] = "$>$" if all(v["cl"] == "beyond" for v in cand) else "$<$"
                rb = []
                for c in REALC:
                    cand = [v for (fa, k, cc, fn, a), v in res.items() if fa == fam and fn == fname and a == an and cc == c and v["forced"] < LOSSFREE and v["cl"] == "on"]
                    if cand:
                        rb.append(min(v["ratio"] for v in cand))
                if rb:
                    M[f"RealBestLo{t}"], M[f"RealBestHi{t}"] = f(min(rb), 2), f(max(rb), 2)
                # real
                rr = [v for (fa, k, c, fn, a), v in res.items() if fa == fam and fn == fname and a == an and c in REALC and v["forced"] < LOSSFREE and v["cl"] == "on"]
                if rr:
                    vs = [v["ratio"] for v in rr]
                    M[f"RealMed{t}"], M[f"RealLo{t}"], M[f"RealHi{t}"], M[f"RealN{t}"] = f(float(np.median(vs)), 2), f(min(vs), 2), f(max(vs), 2), str(len(vs))
    # the same comparison on route changes instead of acquisitions: ISL-only path updates count here
    for fam, sets in FAM.items():
        vs = []
        by = {c: [] for c in CELLS}
        for sp in sets:
            for c in CELLS:
                o = compare(C, c, sp, LIFE, "changes_per_h", "stretch_mean")
                if o and o["cl"] == "on" and o["forced"] < LOSSFREE:
                    vs.append(o["ratio"]); by[c].append(o["ratio"])
        if vs:
            M[f"ChgMed{fam}"], M[f"ChgLo{fam}"], M[f"ChgHi{fam}"] = f(float(np.median(vs)), 2), f(min(vs), 2), f(max(vs), 2)
        M[f"ChgBest{fam}"] = f(float(np.median([min(v) for v in by.values() if v])), 2)
    for g_ in (.3, 1., 3.):
        vs = []
        for c in CELLS:
            o = compare(C, c, dict(policy="isasr", gamma=g_), LIFE, "changes_per_h", "stretch_mean")
            if o and o["cl"] == "on":
                vs.append(o["ratio"])
        if vs:
            M["ChgIsasr" + {.3: "Low", 1.: "One", 3.: "High"}[g_]] = f(float(np.median(vs)), 2)
    # inference: independent units only, selection frontier, both axes; Holm over everything tested
    for an in ("mean", "p99"):
        keys = [k for k, v in res.items() if k[3] == "sel" and k[4] == an and testable(k[2]) and v["forced"] < LOSSFREE and v["cl"] == "on" and v["p"] == v["p"] and k[0] != "Trk"]
        adj = holm([res[k]["p"] for k in keys]) if keys else []
        tag = "Mean" if an == "mean" else "Pn"
        M[f"Tests{tag}"] = str(len(keys))
        M[f"SigFor{tag}"] = str(sum(1 for k, a in zip(keys, adj) if a < .05 and res[k]["ratio"] > 1))
        M[f"SigAgainst{tag}"] = str(sum(1 for k, a in zip(keys, adj) if a < .05 and res[k]["ratio"] < 1))
        M[f"AboveOne{tag}"] = str(sum(1 for k in keys if res[k]["ratio"] > 1))
        for fam in ("Elev", "Merged", "Isasr", "Alpr"):
            kk = [(k, a) for k, a in zip(keys, adj) if k[0] == fam]
            M[f"Tests{fam}{tag}"], M[f"SigFor{fam}{tag}"], M[f"SigAgainst{fam}{tag}"] = str(len(kk)), str(sum(1 for k, a in kk if a < .05 and res[k]["ratio"] > 1)), str(sum(1 for k, a in kk if a < .05 and res[k]["ratio"] < 1))
    # noisy track against noise-free track and against the search
    for nm, cfg in NS.items():
        for ee, et in ((.05, "E"), (.10, "T")):
            sp = dict(cfg, eps=ee)
            for c, ct in allc:
                o = compare(C, c, sp, LIFE, "acq_per_h", "stretch_mean")
                if o:
                    M[f"{nm}{et}R{ct}"] = f(o["ratio"], 2) if o["cl"] == "on" else {"beyond": "$>$", "below": "$<$", "dom": "D"}[o["cl"]]
    allr = []
    for nm, cfg in NS.items():
        vals = []
        for ee in (.05, .10):
            for c in CELLS + REALC:
                o = compare(C, c, dict(cfg, eps=ee), LIFE, "acq_per_h", "stretch_mean")
                if o and o["cl"] == "on":
                    vals.append(o["ratio"])
        if vals:
            M[f"{nm}Med"], M[f"{nm}Lo"], M[f"{nm}Hi"] = f(float(np.median(vals)), 2), f(min(vals), 2), f(max(vals), 2)
    # fresh starts (m1): Starlink and Kuiper, same comparison
    fresh = {}
    fC = Cells([dict(r, part="main") for r in R if r["part"] == "fresh"])
    for fam, sets in FAM.items():
        vs = []
        for sp in sets:
            for c in CELLS[:4]:
                o = compare(fC, c, sp, LIFE, "acq_per_h", "stretch_mean")
                if o and o["cl"] == "on" and o["forced"] < LOSSFREE:
                    vs.append(o["ratio"])
        if vs:
            M[f"FreshMed{fam}"], M[f"FreshLo{fam}"], M[f"FreshHi{fam}"] = f(float(np.median(vs)), 2), f(min(vs), 2), f(max(vs), 2)
    json.dump({"|".join(map(str, (k[0], k[1], k[2][0], k[2][1], k[3], k[4]))): {a: (b if isinstance(b, (str, int, float)) else None) for a, b in v.items() if a != "sp"} for k, v in res.items()},
              open(os.path.join(OUT, "front_tests.json"), "w"), indent=0)
    # equal change rate: merged windows against the frontier
    for win, tag in ((180.0, "AsrThree"), (120.0, "AsrLong")):
        sp = dict(policy="asr", window=win)
        fac, b = [], {}
        for c, ct in CT.items():
            ss = C.starts(c, sp)
            use = [q for q in LIFE if all(C.row(c, s, q) for s in ss)]
            pts = fr(C, c, use, ss, "acq_per_h", "stretch_mean")
            rr = [C.row(c, s, sp) for s in ss]
            sx, ny = np.mean([100 * r["stretch_mean"] for r in rr]), np.mean([r["acq_per_h"] for r in rr])
            cand = [x for x, y in pts if y <= ny]
            if cand:
                fac.append(sx / min(cand))
                M[f"EqRate{tag}{ct}A"], M[f"EqRate{tag}{ct}B"] = f(sx, 1), f(min(cand), 1)
        if fac:
            M[f"EqRate{tag}Lo"], M[f"EqRate{tag}Hi"], M[f"EqRate{tag}Med"] = f(min(fac), 1), f(max(fac), 1), f(float(np.median(fac)), 1)
    return res
