#!/usr/bin/env python3
"""Data figures of the TAES manuscript, all from out/p7_v7.jsonl (one binary) and the closed-form checks."""
import json, os, sys
import numpy as np
import isl_core as C
from figstyle import *  # noqa: F401,F403
from p7util import load, match, by_start, frontier, CELLS, strat, OUT, HERE
import matplotlib.pyplot as plt
import run_p7 as _R7
S6 = set(_R7.S6)

setup()
FIG = os.path.join(HERE, "..", "figures")
SHELLS = {"starlink": ("Starlink", 72, 22, 39, 53.0, 550e3, False), "kuiper": ("Kuiper", 34, 34, 1, 51.9, 630e3, False),
          "oneweb": ("OneWeb", 12, 49, 1, 87.9, 1200e3, True)}
CN = {("starlink", "lon-nyc"): "Starlink, London-New York", ("starlink", "syd-lax"): "Starlink, Sydney-Los Angeles",
      ("kuiper", "lon-nyc"): "Kuiper, London-New York", ("kuiper", "syd-lax"): "Kuiper, Sydney-Los Angeles",
      ("oneweb", "lon-nyc"): "OneWeb, London-New York", ("oneweb", "syd-lax"): "OneWeb, Sydney-Los Angeles"}
PARTS = ("main", "front", "extra")


def save(fig, name):
    fig.savefig(os.path.join(FIG, name)); plt.close(fig); print("wrote", name)


def fig_law():
    """Range of three +Grid inter-plane links over one orbit: Lemma 1 lines, the simulator's propagator as markers."""
    D = [json.loads(l) for l in open(os.path.join(OUT, "law_check.jsonl"))]
    fig, axs = plt.subplots(1, 3, figsize=(PAGE, 2.0), sharey=True)
    for ax, (s, (lab, P, S, F, i, h, star)) in zip(axs, SHELLS.items()):
        sh = C.walker(P, S, F, i, h, star)
        u = np.linspace(0, 2 * np.pi, 401)
        for col, r in zip((BLUE, ORANGE, GREEN), [x for x in D if x["shell"] == s]):
            a, b = r["a"], r["b"]
            dO, Du = sh["Om"][b] - sh["Om"][a], sh["u0"][b] - sh["u0"][a]
            M_, N_ = C.MN(dO, Du, sh["i"])
            d = sh["r"] * np.sqrt(np.maximum(2 - 2 * (M_ - N_ * np.cos(2 * (sh["u0"][a] + u) + Du)), 0)) / 1e3
            ax.plot(np.degrees(u), d, color=col, lw=1.0)
            k = np.arange(0, 201, 12)
            ax.plot(360 * k / 200, np.array(r["range"])[k] / 1e3, ls="none", marker="o", mfc="white", mec=col, ms=3)
        rc = C.rcap_eff(sh["r"], 5000e3, 80e3) / 1e3
        ax.axhline(rc, **THRESH)
        ax.set_title(lab, fontsize=8); ax.set_xlabel("Orbit phase (deg)"); ax.set_xticks([0, 90, 180, 270, 360])
        ax.set_ylim(0, 5600)
    axs[0].set_ylabel("Inter-plane range (km)")
    axs[0].text(10, 5150, r"cap $R^{\ast}$", fontsize=7, color=RED)
    axs[2].plot([], [], color=INK2, lw=1.0, label="Lemma 1"); axs[2].plot([], [], ls="none", marker="o", mfc="white", mec=INK2, label="ns-3 propagator")
    axs[2].legend(loc="center right")
    fig.tight_layout(w_pad=0.6)
    save(fig, "fig2_law.pdf")


FAM = [("Lifetime search", [dict(policy="life", eps=e) for e in (0., .01, .02, .03, .05, .07, .10, .15, .20)], BLUE, "o", "-"),
       ("Tracked lifetime search", [dict(policy="track", eps=e) for e in (.02, .05, .07, .10, .15)], TEAL, "P", "--"),
       (r"Elevation rank, 24$^\circ$", [dict(policy="elevb", eps=e, trigger=24.0) for e in (.02, .05, .10)], INK, "v", ":"),
       ("Merged snapshot", [dict(policy="asr", window=w) for w in (30., 60., 120., 180., 300.)], GOLD, "D", "-."),
       ("ISASR", [dict(policy="isasr", gamma=g) for g in (.3, 1., 3.)], RED, "X", "none"),
       ("ALPR", [dict(policy="alpr")], PURPLE, "s", "none"),
       ("On-demand", [dict(policy="ondemand")], GREY, "^", "none")]


def load_corr():
    import make_macros_corr as MC  # noqa: F401 (builds the merged row set and the ISASR matcher)
    return MC.R, MC.match2


def fig_acq():
    """Correspondence: ground-link acquisitions per hour against mean stretch, per development cell, means over 6 starts."""
    R, m2 = load_corr()
    import frontier_stats as F
    from make_floor import FLOOR
    fam = [("Lifetime search", F.LIFE, BLUE, "o", "-"),
           (r"Elevation rank, 24$^\circ$", [dict(policy="elevb", eps=e, trigger=24.0) for e in (.02, .05, .10)] + [dict(policy="elev", trigger=24.0)], INK, "v", (0, (1, 1.5))),
           ("Merged snapshot", [dict(policy="asr", window=w) for w in (30., 60., 120., 180.)], GOLD, "D", "-."),
           ("ISASR, Alg. 3", [dict(policy="isasr", gamma=g, isasrIslAct=1) for g in (0.3, 1.0, 3.0, 1000.0)], RED, "X", "--"),
           ("ISASR, ISLs cost-free", [dict(policy="isasr", gamma=g) for g in (0.3, 1.0, 3.0)], PURPLE, "P", (0, (5, 1.2, 1, 1.2, 1, 1.2)))]
    fig, axs = plt.subplots(2, 3, figsize=(PAGE, 2.2), sharex=True, sharey=True)
    for ax, c in zip(axs.T.ravel(), CELLS):
        for lab, Fs, col, mk, ls in fam:
            pts = []
            for sp in Fs:
                rr = [r for r in R if r["shell"] == c[0] and r["pair"] == c[1] and r["part"] in ("main", "front", "trig") and r["start"] in S6 and m2(r, sp)]
                if rr:
                    pts.append((100 * np.mean([r["stretch_mean"] for r in rr]), np.mean([r["acq_per_h"] for r in rr])))
            pts.sort()
            if pts:
                x, y = zip(*pts)
                ax.plot(x, y, color=col, ls=ls, lw=0.9, marker=mk, ms=3.5, label=lab)
        ax.axhline(FLOOR[c[0]], color="0.55", lw=0.6, ls="-")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(CN[c], fontsize=7.5); ax.set_xlim(0.8, 40); ax.set_ylim(4.5, 140)
        ax.set_xticks([1, 2, 5, 10, 20]); ax.set_xticklabels(["1", "2", "5", "10", "20"]); ax.set_yticks([5, 10, 20, 50, 100]); ax.set_yticklabels(["5", "10", "20", "50", "100"]); ax.minorticks_off()
    for ax in axs[1]:
        ax.set_xlabel("Mean delay stretch (%)")
    fig.supylabel("Ground-link acquisitions per hour", fontsize=7.5, x=0.005)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.04), fontsize=6.5)
    fig.tight_layout(w_pad=0.4, h_pad=0.6); save(fig, "fig3_policies.pdf")


def fig_policies():
    """Route changes per hour against mean stretch, per cell: the frontier of each family, means over starts.
    Hollow: the route broke before it was replaced (forced changes), which is what loses packets; control-plane rows have no loss count."""
    R = load()
    from make_floor import FLOOR
    fig, axs = plt.subplots(2, 3, figsize=(PAGE, 3.9), sharex=True, sharey=True)
    hs = None
    for ax, c in zip(axs.T.ravel(), CELLS):
        for lab, F, col, mk, ls in FAM:
            pts = []
            for sp in F:
                a = by_start(R, sp, c, ("main", "front")); b = by_start(R, sp, c, ("main", "front"), "stretch_mean", 100); l = by_start(R, sp, c, ("main", "front"), "forced_per_h")
                if a:
                    pts.append((np.mean(list(b.values())), np.mean(list(a.values())), np.mean(list(l.values())) > 0))
            pts.sort()
            if not pts:
                continue
            x, y, lost = zip(*pts)
            if ls != "none":
                ax.plot(x, y, color=col, lw=0.9, ls=ls)
            for xi, yi, li in pts:
                ax.plot([xi], [yi], ls="none", marker=mk, ms=4, color=col, mfc="white" if li else col, mec=col)
            ax.plot([], [], ls=ls if ls != "none" else "none", color=col, marker=mk, ms=4, lw=0.9, label=lab)
        ax.axhline(FLOOR[c[0]], **THRESH)
        ax.text(0.97, FLOOR[c[0]] * 1.05, "floor", transform=ax.get_yaxis_transform(), fontsize=6, color=RED, ha="right")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(CN[c], fontsize=7.5)
        ax.set_xlim(0.8, 40); ax.set_ylim(6, 140)
        ax.set_xticks([1, 2, 5, 10, 20]); ax.set_xticklabels(["1", "2", "5", "10", "20"])
        ax.set_yticks([10, 20, 50, 100]); ax.set_yticklabels(["10", "20", "50", "100"]); ax.minorticks_off()
    for ax in axs[1]:
        ax.set_xlabel("Mean delay stretch (%)")
    for ax in axs[:, 0]:
        ax.set_ylabel("Route changes per hour")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=7, bbox_to_anchor=(0.5, 1.04), fontsize=6.5, handlelength=2.2, columnspacing=1.0)
    fig.tight_layout(w_pad=0.4, h_pad=0.6)
    save(fig, "fig3_policies.pdf")


def fig_stale():
    """Cell-mean loss of the reactive policies against Proposition 2 fed by a route trace from held-out start instants."""
    d = json.load(open(os.path.join(OUT, "stale_check.json")))
    p, m, lab = np.array(d["pred"]), np.array(d["meas"]), np.array(d["lab"])
    fig, ax = plt.subplots(figsize=(COL, 2.3))
    lim = [0.5, 4.0]
    ax.plot(lim, lim, color=INK2, lw=0.8, label="Proposition 1")
    ax.plot(p[lab == "Od"], m[lab == "Od"], ls="none", marker="^", mfc="none", mec=BLUE, ms=5, label="On-demand, ns-3")
    ax.plot(p[lab == "SnapSlow"], m[lab == "SnapSlow"], ls="none", marker="s", mfc="none", mec=ORANGE, ms=5, label="Snapshot 1 s, ns-3")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xticks([0.5, 1, 2, 4]); ax.set_xticklabels(["0.5", "1", "2", "4"]); ax.set_yticks([0.5, 1, 2, 4]); ax.set_yticklabels(["0.5", "1", "2", "4"]); ax.minorticks_off()
    ax.set_xlabel("Predicted loss from the route trace (%)"); ax.set_ylabel("Measured packet loss (%)")
    ax.legend(loc="upper left", fontsize=6.5)
    fig.tight_layout(); save(fig, "fig4_stale.pdf")


def fig_tacq():
    """Loss against acquisition time: two terminals and one terminal per gateway, ns-3 markers (95 % t), Prop. 2 lines."""
    from p7util import tci
    R = load()
    c = ("starlink", "lon-nyc")
    fig, axs = plt.subplots(1, 2, figsize=(COL * 1.0, 2.2), sharey=True)
    tas = np.array([0.0, 0.2, 1.0, 3.0])
    pols = [("On-demand", dict(policy="ondemand"), PURPLE, "^"), ("Snapshot 1 s", dict(policy="snap", period=1.0), GREY, "s"),
            ("Elevation rank", dict(policy="elevb", eps=.05, trigger=24.0), INK, "v"), ("Tracked search", dict(policy="track", eps=.05), TEAL, "P"),
            ("Lifetime search", dict(policy="life", eps=.05), BLUE, "o")]
    for ax, bbm in zip(axs, (0, 1)):
        for lab, sp, col, mk in pols:
            xs, ys, es = [], [], []
            for ta in tas:
                rr = [r for r in R if r["shell"] == c[0] and r["pair"] == c[1] and match(r, sp) and abs(r["cfgd"].get("tAcq", -1) - ta) < 1e-9
                      and (r["part"] == "bbm" if bbm else r["part"] in ("main", "tacq"))]
                if rr:
                    m, h = tci([100 * r["loss"] for r in rr]); xs.append(ta); ys.append(max(m, 1e-3)); es.append(0 if h != h else h)
            if xs:
                ax.errorbar(xs, ys, yerr=es, ls="none", marker=mk, ms=3.5, color=col, mfc="none", capsize=1.5, lw=0.7, label=lab)
                if (lab in ("On-demand", "Snapshot 1 s") and not bbm) or (bbm and lab != "Snapshot 1 s"):
                    # analysis line: Prop. 2 with the rates of the t_a = 1 s runs in this cell
                    q = [r for r in R if r["part"] == "main" and r["shell"] == c[0] and r["pair"] == c[1] and match(r, sp)]
                    if q:
                        if lab == "On-demand":
                            lam = np.mean([r["breaks_per_h"] for r in q]) / 3600; W = np.mean([r["route_ms_mean"] for r in q]) / 1e3
                            yl = 100 * lam * (0.05 + tas + W / 2)
                        elif lab == "Snapshot 1 s":
                            lam = np.mean([r["breaks_per_h"] for r in q]) / 3600; W = np.mean([r["route_ms_mean"] for r in q]) / 1e3
                            yl = 100 * lam * (0.5 + tas + W / 2)
                        else:
                            yl = 100 * np.mean([r["acq_per_h"] for r in q]) / 3600 * tas
                        ax.plot(tas[tas > 0] if bbm else tas, np.maximum((yl[tas > 0] if bbm else yl), 1e-3), color=col, lw=0.7)
        ax.set_yscale("log"); ax.set_xlabel(r"$t_{\mathrm{a}}$ (s)"); ax.set_title("One terminal" if bbm else "Two terminals", fontsize=7.5); ax.set_ylim(5e-4, 30)
        ax.set_xticks([0, 1, 2, 3])
    axs[0].set_ylabel("Packet loss (%)")
    axs[0].text(0.5, 0.12, "early switchers: zero loss", transform=axs[0].transAxes, fontsize=6, color=INK2, ha="center")
    h, l = axs[1].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, bbox_to_anchor=(0.54, 1.12), fontsize=6, handlelength=1.2, columnspacing=0.8)
    fig.tight_layout(w_pad=0.4); save(fig, "fig6_tacq.pdf")


def fig_stress():
    """Loss against the along-track error of the predictor per shell; means with range over runs; tracked search as a reference."""
    R = load()
    fig, axs = plt.subplots(1, 2, figsize=(COL, 2.2), sharey=True)
    for ax, sh in zip(axs, ("starlink", "oneweb")):
        for gd, col in ((0.5, BLUE), (2.0, GREEN)):
            xs, ys, lo, hi = [], [], [], []
            for sg in (5.0, 20.0):
                rr = [r for r in R if r["part"] == "pred" and r["shell"] == sh and r["cfgd"].get("policy") == "life" and abs(r["cfgd"].get("guard", -1) - gd) < 1e-9 and abs(r["cfgd"].get("predSigmaKm", -1) - sg) < 1e-9]
                if rr:
                    v = 100 * np.array([r["loss"] for r in rr]); xs.append(sg); ys.append(max(v.mean(), 1e-3)); lo.append(ys[-1] - max(v.min(), 1e-3)); hi.append(max(v.max(), 1e-3) - ys[-1])
            if xs:
                ax.errorbar(xs, ys, yerr=[lo, hi], marker="o", ms=3.5, color=col, capsize=2, lw=0.9, label=f"Lifetime, guard {gd:g} s")
        for pol, col, mk, lab in (("asr", GOLD, "d", "Merged 120 s"),):
            ys = []
            xs = []
            for sg in (5.0, 20.0):
                rr = [r for r in R if r["part"] == "pred" and r["shell"] == sh and r["cfgd"].get("policy") == pol and abs(r["cfgd"].get("predSigmaKm", -1) - sg) < 1e-9]
                if rr:
                    xs.append(sg); ys.append(max(100 * np.mean([r["loss"] for r in rr]), 1e-3))
            if xs:
                ax.plot(xs, ys, marker=mk, ms=4, color=col, ls=(0, (4, 1.6)), lw=0.9, label=lab)
        ax.set_yscale("log"); ax.set_xlabel("Along-track error, 1$\\sigma$ (km)"); ax.set_title(SHELLS[sh][0], fontsize=7.5); ax.set_xticks([5, 20]); ax.set_ylim(5e-4, 3)
    axs[0].set_ylabel("Packet loss (%)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, bbox_to_anchor=(0.54, 1.09), fontsize=6, handlelength=1.6, columnspacing=0.8)
    fig.tight_layout(w_pad=0.4); save(fig, "fig5_stress.pdf")


def fig_pred():
    """Correspondence figure: reactive loss in ns-3 against Prop. 2 fed by a route trace; development cells hollow, held-out filled."""
    PP = json.load(open(os.path.join(OUT, "pred_check.json")))
    fig, ax = plt.subplots(figsize=(COL, 1.72))
    lim = [0.3, 4.0]
    ax.plot(lim, lim, color=INK2, lw=0.8, label="Proposition 1")
    for pol, col, mk, lab in (("Od", BLUE, "^", "On-demand"), ("Snap", ORANGE, "s", "Snapshot 1 s")):
        d = [p for p in PP if p["pol"] == pol and p["kind"] == "design"]
        h = [p for p in PP if p["pol"] == pol and p["kind"] == "held"]
        ax.plot([p["pred"] for p in d], [p["meas"] for p in d], ls="none", marker=mk, mfc="none", mec=col, mew=1.0, ms=5.5, label=f"{lab}, development")
        if h:
            ax.errorbar([p["pred"] for p in h], [p["meas"] for p in h],
                        yerr=[[p["meas"] - p["measlo"] for p in h], [p["meashi"] - p["meas"] for p in h]],
                        ls="none", marker=mk, color=col, ms=5.5, capsize=1.5, lw=0.7, label=f"{lab}, held out")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    t = [0.5, 1, 2, 4]
    ax.set_xticks(t); ax.set_xticklabels([str(x) for x in t]); ax.set_yticks(t); ax.set_yticklabels([str(x) for x in t]); ax.minorticks_off()
    ax.set_xlabel("Predicted loss from the route trace (%)"); ax.set_ylabel("Measured packet loss (%)")
    ax.legend(loc="upper left", fontsize=6.8, handlelength=1.2, borderaxespad=0.3)
    fig.tight_layout(); save(fig, "fig4_stale.pdf")


if __name__ == "__main__":
    for w in (sys.argv[1:] or ["law", "policies", "stale", "stress", "tacq"]):
        globals()[f"fig_{w}"]()
