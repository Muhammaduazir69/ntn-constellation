"""Helpers shared by make_macros.py and make_figs.py: one loader for the campaign of the article, policy specs,
cell statistics and equal-stretch comparison. All rows come from out/p7_v8.jsonl, from one binary."""
import json, os
from collections import defaultdict
import numpy as np
from scipy import stats as st

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SHELLS = ("starlink", "kuiper", "oneweb")
PAIRS = ("lon-nyc", "syd-lax")
CELLS = [(s, p) for s in SHELLS for p in PAIRS]


def load(path=None):
    rows, shas = [], set()
    seen = {}
    for l in open(path or os.path.join(OUT, "p7_v8.jsonl")):
        r = json.loads(l)
        if "error" in r:
            continue
        r["cfgd"] = dict(r["cfg"])
        shas.add(r["bin_sha"])
        seen[(r["part"], r["shell"], r["pair"], r["start"], json.dumps(r["cfg"]))] = r
    assert len(shas) <= 1, shas
    return list(seen.values())


def match(r, spec, part=None):
    """A row matches when its configuration equals spec on spec's keys, and has no other knob than tAcq/noData/bbm."""
    d = r["cfgd"]
    if part is not None and r["part"] not in (part if isinstance(part, (tuple, list, set)) else (part,)):
        return False
    for k, v in spec.items():
        x = d.get(k)
        if x is None or (abs(x - v) > 1e-9 if isinstance(v, (int, float)) else x != v):
            return False
    for k in ("hold", "proactive", "trigger", "gamma", "guard", "predSigmaKm", "crossLatDeg", "floorGrid", "elevNoiseDeg", "trackSamples", "trackArc"):
        if k in d and k not in spec and not (k == "gamma" and d[k] == 1.0):
            return False
    return True


def by_start(rows, spec, cell, part="main", key="changes_per_h", scale=1.0):
    return {r["start"]: r[key] * scale for r in rows if r["shell"] == cell[0] and r["pair"] == cell[1] and match(r, spec, part)}


def tci(v):
    """Mean and 95 % t half-width."""
    v = np.asarray(v, float)
    if len(v) < 2:
        return float(v.mean()) if len(v) else float("nan"), float("nan")
    return float(v.mean()), float(st.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v)))


def strat(rows, spec, key, scale=1.0, part="main", cells=CELLS):
    """Mean stratified over cells, with a 95 % t half-width from the within-cell variances (Welch-Satterthwaite df)."""
    ms, vs, ns = [], [], []
    for c in cells:
        v = np.array(list(by_start(rows, spec, c, part, key, scale).values()), float)
        if len(v) == 0:
            continue
        ms.append(v.mean()); vs.append(v.var(ddof=1) if len(v) > 1 else 0.0); ns.append(len(v))
    if not ms:
        return None
    k = len(ms)
    m = float(np.mean(ms))
    var = sum(a / b for a, b in zip(vs, ns)) / k ** 2
    den = sum((a / b) ** 2 / (b - 1) for a, b in zip(vs, ns) if b > 1)
    df = (sum(a / b for a, b in zip(vs, ns)) ** 2 / den) if den > 0 else 1.0
    return m, float(st.t.ppf(0.975, max(df, 1.0)) * np.sqrt(var))


def frontier(rows, family, cell, start, parts=("main", "front")):
    """Stretch and change rate of one policy family at one start, sorted by stretch: list of (stretch %, changes/h)."""
    pts = []
    for spec in family:
        for r in rows:
            if r["shell"] == cell[0] and r["pair"] == cell[1] and r["start"] == start and match(r, spec, parts):
                pts.append((100 * r["stretch_mean"], r["changes_per_h"]))
                break
    return sorted(pts)


def interp(pts, x):
    """Change rate of a frontier at stretch x; None outside its range (no extrapolation)."""
    if not pts or x < pts[0][0] - 1e-9 or x > pts[-1][0] + 1e-9:
        return None
    xs, ys = zip(*pts)
    return float(np.interp(x, xs, ys))


def holm(ps):
    """Holm-adjusted p-values, in input order."""
    ps = np.asarray(ps, float)
    o = np.argsort(ps); adj = np.zeros(len(ps)); run = 0.0
    for i, k in enumerate(o):
        run = max(run, (len(ps) - i) * ps[k]); adj[k] = min(1.0, run)
    return adj
