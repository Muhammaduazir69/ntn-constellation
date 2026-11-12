#!/usr/bin/env python3
"""Exact-model Monte Carlo of the short-arc residual estimate (Lemma 2): circular Walker orbits under a rotating
spherical Earth with the gateway on the WGS84 ellipsoid, as in the simulator up to J2. Random visible satellite and time;
the true residual to the drop elevation against the parabola through psi^2 at t, t-15 s, t-30 s. Writes out/track_check.json."""
import json, numpy as np, isl_core as C
from check_track_core import SH, GW, WE, DROP, UP, gw, elev, psi

res = {}
rng = np.random.default_rng(11)
for name, (P, S, F, inc, h, star) in SH.items():
    sh = C.walker(P, S, F, inc, h, star); errs = []
    for _, (la, lo) in GW.items():
        g, up = gw(la, lo)
        for t0 in rng.uniform(0, 10 * 86400, 60):
            e0 = elev(sh, np.arange(P * S), t0, g, up)
            for k in np.where(e0 >= UP)[0]:
                ts = t0 + np.arange(0, 1500, 0.5)
                e = elev(sh, k, ts, g, up)
                below = np.where(e < DROP)[0]
                if not len(below):
                    continue
                true = ts[below[0]] - t0
                y = [psi(elev(sh, k, t0 - x, g, up), sh["r"]) ** 2 for x in (0.0, 15.0, 30.0)]
                c = (y[2] - 2 * y[1] + y[0]) / (2 * 15.0 ** 2); b = (y[0] - y[1]) / 15.0 + c * 15.0
                a = y[0] - psi(DROP, sh["r"]) ** 2
                est = (-b + np.sqrt(max(0.0, b * b - 4 * a * c))) / (2 * c) if c > 0 and a < 0 else 0.0
                errs.append(est - true)
    errs = np.array(errs)
    res[name] = dict(n=len(errs), med_abs=float(np.median(np.abs(errs))), p95_abs=float(np.percentile(np.abs(errs), 95)), max_abs=float(np.max(np.abs(errs))), bias=float(errs.mean()))
    print(name, res[name], flush=True)
json.dump(res, open("out/track_check.json", "w"), indent=1)
