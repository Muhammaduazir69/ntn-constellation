"""How well does the circular Walker-pair law predict ISL range and cap crossings on real Starlink orbits?

Pairs are drawn from the 53.16 deg shell near 463 km, restricted to altitudes within 5 km of its median: each satellite with its nearest satellite in an
adjacent plane (node difference 3 to 8 deg), the pair a +Grid inter-plane link would use. The law is
fed the osculating node, argument of latitude, inclination and radius at t0 (pair means for i and r),
and its range and cap-crossing epochs over one orbit are compared with SGP4.
"""
import json, numpy as np, sys
from sgp4.api import Satrec, jday
import isl_core as C

def load(fn):
    L = [l.rstrip() for l in open(fn) if l.strip()]
    out = []
    for k in range(0, len(L) - 2, 3):
        s = Satrec.twoline2rv(L[k + 1], L[k + 2])
        inc = np.degrees(s.inclo); n = s.no_kozai / 60
        h = ((C.MU / n ** 2) ** (1 / 3) - C.RE) / 1e3
        if 53.10 < inc < 53.25 and 455 < h < 480:   # the current 53.16 deg shell
            out.append((s, h))
    # one shell: altitudes within 5 km of the median of the 53.0 deg satellites
    med = np.median([h for _, h in out])
    return [s for s, h in out if abs(h - med) <= 5.0]

def state(s, jd, fr):
    e, r, v = s.sgp4(jd, fr)
    return (np.array(r) * 1e3, np.array(v) * 1e3) if e == 0 else (None, None)

def elems(r, v):
    h = np.cross(r, v); i = np.arccos(h[2] / np.linalg.norm(h))
    O = np.arctan2(h[0], -h[1]); nvec = np.array([np.cos(O), np.sin(O), 0])
    u = np.arctan2(r @ np.cross(h / np.linalg.norm(h), nvec), r @ nvec)
    return i, O, np.mod(u, 2 * np.pi), np.linalg.norm(r)

def main(cap_km=1700.0, npairs=300, seed=1):
    S = load("data/starlink_20261002.tle")
    jd, fr = jday(2026, 10, 2, 12, 0, 0)
    st = [state(s, jd, fr) for s in S]
    keep = [k for k, x in enumerate(st) if x[0] is not None]
    E = {k: elems(*st[k]) for k in keep}
    rng = np.random.default_rng(seed)
    pairs = []
    for a in rng.permutation(keep):
        ia, Oa, ua, ra = E[a]
        best = None
        for b in keep:
            dO = np.degrees(np.mod(E[b][1] - Oa + np.pi, 2 * np.pi) - np.pi)
            if 3 <= dO <= 8:
                d = np.linalg.norm(st[a][0] - st[b][0])
                if best is None or d < best[0]: best = (d, b)
        if best: pairs.append((a, best[1]))
        if len(pairs) >= npairs: break
    T = 2 * np.pi * np.sqrt(np.mean([E[k][3] for k in keep]) ** 3 / C.MU)
    ts = np.arange(0, T, 1.0)
    rerr, everr, nev_true, nev_pred = [], [], 0, 0
    for a, b in pairs:
        ia, Oa, ua, ra = E[a]; ib, Ob, ub, rb = E[b]
        inc, r = (ia + ib) / 2, (ra + rb) / 2
        w = 2 * np.pi / T
        M, N = C.MN(Ob - Oa, ub - ua, inc)
        dpred = r * np.sqrt(2 - 2 * (M - N * np.cos(2 * (ua + w * ts) + (ub - ua))))
        dtrue = np.array([np.linalg.norm(state(S[a], jd, fr + t / 86400)[0] - state(S[b], jd, fr + t / 86400)[0]) for t in ts])
        rerr.append(np.max(np.abs(dpred - dtrue)))
        up_t, up_p = dtrue <= cap_km * 1e3, dpred <= cap_km * 1e3
        et = ts[1:][np.diff(up_t.astype(int)) != 0]; ep = ts[1:][np.diff(up_p.astype(int)) != 0]
        nev_true += len(et); nev_pred += len(ep)
        for x in et:
            if len(ep): everr.append(np.min(np.abs(ep - x)))
    res = dict(shell_sats=len(keep), pairs=len(pairs), period_s=T, cap_km=cap_km,
               range_err_km_median=float(np.median(rerr) / 1e3), range_err_km_p95=float(np.percentile(rerr, 95) / 1e3),
               range_err_km_max=float(np.max(rerr) / 1e3), events_true=nev_true, events_pred=nev_pred,
               epoch_err_s_median=float(np.median(everr)) if everr else None,
               epoch_err_s_p95=float(np.percentile(everr, 95)) if everr else None)
    json.dump(res, open("out/real_shell.json", "w"), indent=1); print(res)

if __name__ == "__main__":
    main(*(float(x) for x in sys.argv[1:2]))
