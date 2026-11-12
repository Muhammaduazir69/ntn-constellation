#!/usr/bin/env python3
"""Lemma 2 with noisy elevations (analysis only, no simulator): 3 samples over 30 s, and a least-squares parabola through 30 samples
over 60 s, Gaussian elevation noise 0.01 and 0.05 deg per sample. Same exact-model passes as check_track.py. Writes out/track_noise.json."""
import json, numpy as np, isl_core as C
from check_track_core import SH, GW, gw, elev, psi, DROP, UP
rng = np.random.default_rng(12)
res = {}
for name in ("starlink", "oneweb"):
    P, S, F, inc, h, star = SH[name]
    sh = C.walker(P, S, F, inc, h, star)
    pairs = []
    for la, lo in GW.values():
        g, up = gw(la, lo)
        for t0 in rng.uniform(0, 10 * 86400, 12):
            for k in np.where(elev(sh, np.arange(P * S), t0, g, up) >= UP)[0]:
                ts = t0 + np.arange(0, 1500, 0.5); e = elev(sh, k, ts, g, up); b = np.where(e < DROP)[0]
                if len(b):
                    pairs.append((k, t0, g, up, ts[b[0]] - t0))
    pd2 = psi(DROP, sh["r"]) ** 2
    for sig in (0.01, 0.05):
        for tag, off in (("three", np.array([0.0, 15.0, 30.0])), ("lsq", np.linspace(0.0, 60.0, 30))):
            errs, bad = [], 0
            for k, t0, g, up, true in pairs:
                y = psi(elev(sh, k, t0 - off, g, up) + rng.normal(0, sig, len(off)), sh["r"]) ** 2
                c2, c1, c0 = np.polyfit(-off, y, 2)
                a = c0 - pd2
                if c2 > 0 and a < 0:
                    errs.append((-c1 + np.sqrt(max(0.0, c1 * c1 - 4 * a * c2))) / (2 * c2) - true)
                else:
                    bad += 1
            errs = np.abs(errs)
            res[f"{name}_{tag}_{sig}"] = dict(n=len(pairs), med=float(np.median(errs)), p95=float(np.percentile(errs, 95)), bad=bad / len(pairs))
            print(name, tag, sig, res[f"{name}_{tag}_{sig}"], flush=True)
json.dump(res, open("out/track_noise.json", "w"), indent=1)
