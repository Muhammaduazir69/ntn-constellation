"""Largest +Grid range over every phasing factor F = 0..P-1 of each shell (Lemma 1); writes out/grid_allF.json."""
import json, numpy as np, isl_core as C
from route_proto import plusgrid

SHELLS = {"starlink": (72, 22, 53, 550e3, False), "kuiper": (34, 34, 51.9, 630e3, False),
          "oneweb": (12, 49, 87.9, 1200e3, True)}
res = {}
for name, (P, S, i, h, star) in SHELLS.items():
    worst = 0.0
    for F in range(P):
        sh = C.walker(P, S, F, i, h, star); E = plusgrid(sh)
        dO = sh["Om"][E[:, 1]] - sh["Om"][E[:, 0]]; D = sh["u0"][E[:, 1]] - sh["u0"][E[:, 0]]
        M, N = C.MN(dO, D, sh["i"])
        worst = max(worst, (sh["r"] * np.sqrt(2 - 2 * (M - N))).max())
    res[name] = worst / 1e3
    print(name, "largest +Grid range over all F: %.0f km" % res[name])
json.dump(res, open("out/grid_allF.json", "w"))
