#!/usr/bin/env python3
"""Our PYTHIA 8.313 vs ALICE JHEP07(2023)136 Fig. 7 (pp 13 TeV): (Lambda + anti-Lambda)/2K0S inclusive and in jets.

  python3 macros/plot_alice_check.py   -> compare/plots_alice/LamK0s_ALICE_pp13.{pdf,png}, values.txt
Inputs: alice/pp13_minbias_<tune>/*.root (inclusive) and alice/pp13_hardqcd_<tune>/*.root (in jets), bin/aliceCheck.
In jets (JE) = jet cone - perpendicular cone, per species, then the ratio.

ALICE points: READ OFF Fig. 7 BY EYE (approximate, +-0.02); replace with the HEPData tables (ins2180474) when available.
"""
import glob, os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(D, "compare", "plots_alice")
os.makedirs(OUT, exist_ok=True)
TUNES = [("monash", "Monash (no CR-BLC)", "#8a8a8a", "x"), ("crblc_mode0", "CR-BLC mode 0", "#e87ba4", "D"),
         ("crblc_mode2", "CR-BLC mode 2", "#eb6834", "o"), ("crblc_mode3", "CR-BLC mode 3", "#1baf7a", "^")]
# ALICE pp 13 TeV, Fig. 7 top left, read off by eye: (pT, value)
ALICE_INCL = [(0.7, 0.18), (0.9, 0.22), (1.1, 0.27), (1.3, 0.32), (1.5, 0.37), (1.7, 0.41), (1.9, 0.44), (2.1, 0.47),
              (2.3, 0.50), (2.5, 0.51), (2.8, 0.52), (3.2, 0.52), (3.6, 0.51), (4.0, 0.48), (4.5, 0.44), (5.5, 0.34),
              (7.0, 0.28), (10.0, 0.27)]
ALICE_JE = [(1.0, 0.115), (2.0, 0.21), (2.6, 0.24), (3.3, 0.26), (4.3, 0.26), (6.5, 0.29), (10.0, 0.27)]
# ALICE's PYTHIA CR-BLC band (spread of modes 0/2/3), read off by eye: (pT, low, high)
ALICE_PY_INCL = [(0.7, 0.17, 0.19), (1.0, 0.24, 0.26), (1.5, 0.34, 0.37), (2.0, 0.42, 0.46), (2.5, 0.46, 0.50),
                 (3.0, 0.46, 0.50), (4.0, 0.42, 0.46), (5.0, 0.37, 0.40), (6.0, 0.32, 0.35), (8.0, 0.29, 0.31),
                 (10.0, 0.27, 0.29), (12.0, 0.25, 0.27)]
ALICE_PY_JE = [(0.7, 0.13, 0.17), (1.0, 0.17, 0.21), (1.5, 0.21, 0.25), (2.0, 0.24, 0.27), (3.0, 0.27, 0.29),
               (4.0, 0.28, 0.30), (6.0, 0.28, 0.31), (8.0, 0.28, 0.31), (12.0, 0.28, 0.31)]


def load(proc, tune):
    files = glob.glob(os.path.join(D, "alice", f"{proc}_{tune}", "*.root"))
    if not files:
        return None, 0
    acc = {}
    for p in files:
        f = uproot.open(p)
        for k in ("hEv", "hIncl_K0s", "hIncl_Lam", "hJC_K0s", "hJC_Lam", "hPC_K0s", "hPC_Lam"):
            h = f[k]
            v, s2 = h.values(), h.variances()
            if k in acc:
                acc[k][0] += v; acc[k][1] += s2
            else:
                acc[k] = [v.copy(), s2.copy(), h.axis().edges()]
    return acc, len(files)


def ratio(l, sl, k, sk):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = 0.5 * l / k
        return r, r * np.sqrt(sl / l ** 2 + sk / k ** 2)


fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), sharey=True, gridspec_kw={"wspace": 0.06})
lines = []
for ax, proc, title, data, band in [
        (axes[0], "pp13_minbias", "Inclusive (minimum bias)", ALICE_INCL, ALICE_PY_INCL),
        (axes[1], "pp13_hardqcd", "Particles in jets (JC − perp. cone)", ALICE_JE, ALICE_PY_JE)]:
    b = np.array(band)
    ax.fill_between(b[:, 0], b[:, 1], b[:, 2], color="0.75", alpha=0.6, lw=0, label="ALICE's PYTHIA CR-BLC band (read off)")
    d = np.array(data)
    ax.plot(d[:, 0], d[:, 1], "s", color="black", ms=6, label="ALICE data (read off Fig. 7)", zorder=10)
    for tune, lab, col, mk in TUNES:
        acc, nf = load(proc, tune)
        if acc is None:
            continue
        e = acc["hIncl_Lam"][2]
        if proc == "pp13_minbias":
            r, er = ratio(acc["hIncl_Lam"][0], acc["hIncl_Lam"][1], acc["hIncl_K0s"][0], acc["hIncl_K0s"][1])
        else:
            l = acc["hJC_Lam"][0] - acc["hPC_Lam"][0]; sl = acc["hJC_Lam"][1] + acc["hPC_Lam"][1]
            k = acc["hJC_K0s"][0] - acc["hPC_K0s"][0]; sk = acc["hJC_K0s"][1] + acc["hPC_K0s"][1]
            r, er = ratio(l, sl, k, sk)
        ok = np.isfinite(r) & (er < 0.3 * np.abs(r))
        c = 0.5 * (e[1:] + e[:-1])
        lo, hi = np.where(ok, r - er, np.nan), np.where(ok, r + er, np.nan)
        ax.fill_between(e, np.append(lo, lo[-1]), np.append(hi, hi[-1]), step="post", color=col, alpha=0.35, lw=0)
        ax.plot(c[ok], r[ok], ls="none", marker=mk, ms=4.5, color=col, mfc="white", mew=1.2, label=f"this work: {lab}")
        lines.append(f"{proc:13s} {tune:12s} files {nf:3d} | " + " ".join(f"{x:.3f}" if o else "  -  " for x, o in zip(r, ok)))
    ax.set_title(title, fontsize=11, loc="left")
    ax.set_xlabel(r"$p_T$ (GeV)")
    ax.set_xlim(0, 12); ax.set_ylim(0, 0.65)
    ax.grid(color="0.92", lw=0.6)
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)
axes[0].set_ylabel(r"$(\Lambda+\bar{\Lambda})/2K^0_S$, $|\eta|$ < 0.75")
axes[1].tick_params(labelleft=True)
h, l = axes[1].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=3, fontsize=8.8, frameon=False, bbox_to_anchor=(0.5, 1.08))
fig.text(0.5, -0.06, "pp 13 TeV, PYTHIA 8.313; jets: anti-kT R = 0.4 charged (pT > 0.15, |η| < 0.9), "
         "pT,jet > 10 GeV, |η_jet| < 0.35. ALICE values read off JHEP07(2023)136 Fig. 7 by eye (±0.02).",
         ha="center", fontsize=8, color="0.4")
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"LamK0s_ALICE_pp13.{ext}"), dpi=140, bbox_inches="tight")
e = acc["hIncl_Lam"][2] if acc else []
open(os.path.join(OUT, "values.txt"), "w").write("pT edges: " + " ".join(f"{x:g}" for x in e) + "\n" + "\n".join(lines) + "\n")
print("\n".join(lines))
print("->", os.path.join(OUT, "LamK0s_ALICE_pp13.png"))
