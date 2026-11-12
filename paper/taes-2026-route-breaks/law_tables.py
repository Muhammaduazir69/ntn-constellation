"""Closed-form ISL law on the three shells, checked against direct sampling; writes out/law.json."""
import json, numpy as np, isl_core as C
from route_proto import plusgrid

SHELLS = {"starlink": (72, 22, 39, 53.0, 550e3, False), "kuiper": (34, 34, 1, 51.9, 630e3, False),
          "oneweb": (12, 49, 1, 87.9, 1200e3, True)}
WE = 7.2921159e-5
CAP, HATM, MASK = 5000e3, 80e3, 25.0


def sampled(sh, rc, npts=6001):
    n = sh["P"] * sh["S"]; deg = ev = 0.0
    for a in C.refs(sh):
        prev = None
        for u in np.linspace(0, 2 * np.pi, npts):
            X = C.pos(sh, u); up = np.linalg.norm(X - X[a], axis=1) <= rc; up[a] = False
            if prev is not None: ev += np.sum(up != prev)
            deg += up.sum(); prev = up
    R = len(C.refs(sh))
    return deg / (npts * R), ev * n / R / 2


def pass_max(sh, lat_deg, mask=MASK):
    """Longest pass over a site: cap diameter over the slowest ground-track speed (Lemma 1 of the paper)."""
    r = sh["r"]; e = np.radians(mask)
    psi = np.arccos(C.RE * np.cos(e) / r) - e
    w = np.sqrt(C.MU / r ** 3); phi = np.radians(lat_deg)
    Om = np.sqrt(w ** 2 - 2 * w * WE * np.cos(sh["i"]) + (WE * np.cos(phi)) ** 2)
    return 2 * psi / Om


if __name__ == "__main__":
    out = {}
    for name, (P, S, F, i, h, star) in SHELLS.items():
        sh = C.walker(P, S, F, i, h, star); rc = C.rcap_eff(sh["r"], CAP, HATM)
        E = plusgrid(sh)
        dO = sh["Om"][E[:, 1]] - sh["Om"][E[:, 0]]; D = sh["u0"][E[:, 1]] - sh["u0"][E[:, 0]]
        M, N = C.MN(dO, D, sh["i"])
        dmax = sh["r"] * np.sqrt(2 - 2 * (M - N)); dmin = sh["r"] * np.sqrt(np.maximum(2 - 2 * (M + N), 0))
        _, _, evg = C.pair_law(dO, D, sh["i"], sh["r"], rc)
        L = C.shell_law(sh, rc)
        dS, eS = sampled(sh, rc)
        Ev = C.event_times(sh, rc); T = 2 * np.pi * np.sqrt(sh["r"] ** 3 / C.MU)
        ep = np.unique(np.round(np.concatenate([x[0] for x in Ev]) * T, 3)).size
        out[name] = dict(P=P, S=S, F=F, inc=i, alt_km=h / 1e3, n=P * S, period_s=T, rcap_km=rc / 1e3,
                         limb_km=2 * np.sqrt(sh["r"] ** 2 - (C.RE + HATM) ** 2) / 1e3,
                         grid_links=len(E), grid_churning=int(np.sum(evg > 0)),
                         grid_dmax_km=float(dmax.max() / 1e3), grid_dmin_km=float(dmin.min() / 1e3),
                         af_degree=L["mean_degree"], af_degree_sampled=dS, af_events=L["events_per_orbit"],
                         af_events_sampled=eS, af_epochs=ep, af_epoch_gap_ms=1e3 * T / ep,
                         pass_max_s={k: pass_max(sh, lat) for k, lat in
                                     {"lon": 51.51, "nyc": 40.71, "fra": 50.11, "sin": 1.35, "syd": -33.87, "lax": 34.05}.items()})
        print(name, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in out[name].items() if k != "pass_max_s"},
              {k: round(v) for k, v in out[name]["pass_max_s"].items()})
    json.dump(out, open("out/law.json", "w"), indent=1)
