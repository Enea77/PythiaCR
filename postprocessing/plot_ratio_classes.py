#!/usr/bin/env python3
"""(Lambda + anti-Lambda) / 2 K0S vs pT in N_trk^offline classes: p-Pb data vs PYTHIA p-Pb (Angantyr) samples.

  python3 plot_ratio_classes.py inclusive [dijet|mb]   -> plots/ratio_classes/Inclusive_LamK0s_<set>.{pdf,png}
  python3 plot_ratio_classes.py incone                 -> plots/ratio_classes/InCone_LamK0s.{pdf,png}

inclusive: data = v0-jet-analysis/inclusive/plot_inclusive.C on Inclusive_v4, classes 0,20,40,80,185 + HM185, HM250
           (input_data_inclusive/inclusive_spectra_v4_0-20-40-80-185.root; corrected yields, stat. only).
Sets: pp (default; official CMS sample + our pp dijet samples, as compare/plot_overview.py) | ppb | mb (Angantyr).
           MC = all generated V0s with |eta| < 2 of all events (genCompare.C v2: hV0PtAll_*, hV0PtNtrk_*), class from
           N_trk^offline = 0.844 x gen N_ch (hiBin). "dijet": the p-Pb dijet samples (pTHat > 15, weighted, NOT minimum
           bias); "mb": the p-Pb minimum-bias samples (minbiasppb*, when produced).
incone   : data and MC from the jetshape outputs (Integral_vs_pt_innercone, stage sideband_corrected, dr < 0.4, as
           plot_tunes_vs_data.C), classes 0-60, 60-80, 80-185, HM185, HM250; MC = p-Pb dijet samples with pTHat > 30.
MC points with a relative error above MAXREL are not drawn.
"""
import os, sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
COMP = os.path.join(HERE, "..", "compare", "v2")
OUT = os.path.join(HERE, "plots", "ratio_classes")
os.makedirs(OUT, exist_ok=True)
NTRK_SCALE = 0.844
MAXREL = 0.5
MODE = sys.argv[1] if len(sys.argv) > 1 else "incone"
SET = sys.argv[2] if len(sys.argv) > 2 else "pp"   # pp (default) | ppb (Angantyr dijet) | mb (Angantyr minimum bias)
SIG120_FIX = 2.11e-4 / 2.1848571e-4   # official pTHat-120 window, as compare/plot_overview.py

# pp samples, as in compare/plot_overview.py (same colours and markers); official = CMS sample (pp, PYTHIA 8.212)
PP = [("official", "CMS official (CUETP8M1, PYTHIA 8.212)", "#555555", "-"),
      ("dijet_cuetp8m1", "CUETP8M1", "#2a78d6", "-"), ("dijet_crblc_mode0", "CR-BLC mode 0", "#e87ba4", "-"),
      ("dijet_crblc_mode2", "CR-BLC mode 2", "#eb6834", "-"), ("dijet_crblc_mode3", "CR-BLC mode 3", "#e34948", "-")]


# (sample stem, legend, colour, linestyle)
def mcset(prefix, suffix=""):
    return [(f"{prefix}_cuetp8m1{suffix}", "CUETP8M1", "#2a78d6", "-"),
            (f"{prefix}_crblc_mode2{suffix}", "CR-BLC mode 2", "#eb6834", "-"),
            (f"{prefix}_crblc_mode3{suffix}", "CR-BLC mode 3", "#e34948", "-"),
            (f"{prefix}central_cuetp8m1{suffix}", "CUETP8M1, central", "#2a78d6", "--"),
            (f"{prefix}central_crblc_mode2{suffix}", "CR-BLC mode 2, central", "#eb6834", "--")]


def ratio(l, k, el, ek):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = 0.5 * l / k
        return r, r * np.sqrt((el / l) ** 2 + (ek / k) ** 2)


def rebin(v, s2, e, new):
    idx = np.searchsorted(e, new)
    assert np.allclose(e[idx], new), "target edges must be a subset"
    return np.add.reduceat(v, idx[:-1])[:len(new) - 1], np.add.reduceat(s2, idx[:-1])[:len(new) - 1]


# ---------------------------------------------------------------- inputs
if MODE == "inclusive":
    fd = uproot.open(os.path.join(HERE, "input_data_inclusive", "inclusive_spectra_v4_0-20-40-80-185.root"))
    CLASSES = [("MB", "MB integrated", None), ("mult0to20", "0–20", (0, 20)), ("mult20to40", "20–40", (20, 40)),
               ("mult40to80", "40–80", (40, 80)), ("mult80to185", "80–185", (80, 185)),
               ("HM185", "185–250 (HM185)", (185, 250)), ("HM250", "≥ 250 (HM250)", (250, 1e9))]
    SAMPLES = {"pp": PP + [("dijetppb_cuetp8m1", "p–Pb CUETP8M1 (Angantyr)", "#2a78d6", "--")],
               "ppb": mcset("dijetppb"), "mb": mcset("minbiasppb")}[SET]

    def data(cls):
        h = fd[f"LamK0s_{cls}"]
        return h.values(), np.sqrt(h.variances()), h.axis().edges()

    def mc(stem, cls_range, edges):
        if stem == "official":   # 11 pTHat windows; hiBin of the CMS skims is already N_trk^offline
            import glob
            fl = [(uproot.open(p), SIG120_FIX if p.endswith("pthat120.root") else 1.)
                  for p in sorted(glob.glob(os.path.join(COMP, "official", "pthat*.root")))]
            scale = 1.0
        else:
            p = os.path.join(COMP, stem + ".root")
            if not os.path.exists(p):
                return None
            fl, scale = [(uproot.open(p), 1.)], NTRK_SCALE
        out = []
        for sp in ("Lam", "K0s"):
            v = s2 = None
            for f, sc in fl:
                if cls_range is None:
                    h = f[f"hV0PtAll_{sp}"]
                    hv, hs, e = h.values(), h.variances(), h.axis().edges()
                else:
                    h = f[f"hV0PtNtrk_{sp}"]
                    ye = h.axes[1].edges()
                    sel = (scale * ye[:-1] >= cls_range[0]) & (scale * ye[:-1] < cls_range[1])
                    hv, hs, e = h.values()[:, sel].sum(1), h.variances()[:, sel].sum(1), h.axes[0].edges()
                v = sc * hv if v is None else v + sc * hv
                s2 = sc * sc * hs if s2 is None else s2 + sc * sc * hs
            out.append(rebin(v, s2, e, edges))
        (l, sl), (k, sk) = out
        return ratio(l, k, np.sqrt(sl), np.sqrt(sk))

    TITLE = (r"Inclusive $(\Lambda+\bar{\Lambda})/2K^0_S$, $|\eta|$ < 2: p–Pb 8.16 TeV data vs " +
             {"pp": "PYTHIA 8 pp dijet (no nucleus) + p–Pb CUETP8M1 (Angantyr); all events, $\\hat{p}_T$ > 15 GeV, NOT min. bias",
              "ppb": "PYTHIA 8.313 Angantyr p–Pb dijet (all events, $\\hat{p}_T$ > 15 GeV, NOT minimum bias)",
              "mb": "PYTHIA 8.313 Angantyr p–Pb minimum bias"}[SET])
    NAME = f"Inclusive_LamK0s_{SET}"
    XLIM, YLIM = (1, 10), (0, 1.0)
else:
    DD = os.path.join(HERE, "input_data_nominal")
    CLASSES = [("", "MB integrated", "MB"), ("mult0to60_", "0–60", "MB"), ("mult60to80_", "60–80", "MB"),
               ("mult80to185_", "80–185", "MB"), ("mult185to250_", "185–250 (HM185)", "HM185"),
               ("mult250to1000_", "≥ 250 (HM250)", "HM250")]
    SAMPLES = (PP + [("dijetppb_cuetp8m1_pthat30", "p–Pb CUETP8M1 (Angantyr)", "#2a78d6", "--")] if SET == "pp"
               else mcset("dijetppb", "_pthat30"))
    H = "Integral_vs_pt_innercone/Integral_vs_pt_innercone_sideband_corrected"

    def yields(fl, fk):
        if not (os.path.exists(fl) and os.path.exists(fk)):
            return None
        hl, hk = uproot.open(fl)[H], uproot.open(fk)[H]
        r, e = ratio(hl.values(), hk.values(), np.sqrt(hl.variances()), np.sqrt(hk.variances()))
        return r, e, hl.axis().edges()

    def data(cls):
        tag, trig = cls
        return yields(f"{DD}/jetshapes_rebinned_ijc_{tag}reco_reco_Data_Lam_{trig}.root",
                      f"{DD}/jetshapes_rebinned_ijc_{tag}reco_reco_Data_K0s_{trig}.root")

    def mc(stem, tag, edges):
        # pp: MB integrated from input_<sample>, classes from the multiplicity-binned variant input_<sample>_mult
        d = os.path.join(HERE, "input_" + stem + ("_mult" if SET == "pp" and tag and stem.startswith("dijet_") else ""))
        r = yields(f"{d}/jetshapes_rebinned_{tag}gen_gen_MC_Lam.root", f"{d}/jetshapes_rebinned_{tag}gen_gen_MC_K0s.root")
        if r is None:
            return None
        assert np.allclose(r[2], edges)
        return r[0], r[1]

    TITLE = (r"In-cone $(\Lambda+\bar{\Lambda})/2K^0_S$, jets $p_T$ > 60 GeV, $\Delta r$ < 0.4, UE subtracted: p–Pb 8.16 TeV data vs " +
             ("PYTHIA 8 pp dijet (no nucleus) + p–Pb CUETP8M1 (Angantyr, $\\hat{p}_T$ > 30 GeV)" if SET == "pp" else r"PYTHIA 8.313 Angantyr p–Pb dijet ($\hat{p}_T$ > 30 GeV)"))
    NAME = "InCone_LamK0s" + ("_pp" if SET == "pp" else "_ppb")
    XLIM, YLIM = (1.5, 10), (0, 0.5)

# ---------------------------------------------------------------- figure
MARK = {"CMS official (CUETP8M1, PYTHIA 8.212)": "o", "p–Pb CUETP8M1 (Angantyr)": "s", "CR-BLC mode 0": "D", "CUETP8M1": "s", "CR-BLC mode 2": "P", "CR-BLC mode 3": "X", "CUETP8M1, central": "s", "CR-BLC mode 2, central": "P"}
# colours and markers as compare/plot_overview.py; central-enriched samples: same marker, filled, hatched band
MFC = lambda slab, colr: colr if ("central" in slab or "p–Pb" in slab) else "white"


def style(a):
    a.grid(color="0.92", lw=0.6)
    a.minorticks_on()
    a.tick_params(which="both", direction="in", top=True, right=True)
    a.tick_params(which="major", length=6)
    a.tick_params(which="minor", length=3)


def band(a, edges, v, e, colr, hatch, label=None):
    """MC as a band (value +- stat. error per bin); NaN bins are left empty."""
    lo, hi = np.append(v - e, (v - e)[-1]), np.append(v + e, (v + e)[-1])   # step='post': repeat the last bin
    a.fill_between(edges, lo, hi, step="post", color=colr, alpha=0.30 if not hatch else 0.0, lw=0)
    if hatch:   # central-enriched samples: hatched band with an outline
        a.fill_between(edges, lo, hi, step="post", facecolor="none", edgecolor=colr, hatch=hatch, lw=0.8)


ncl = len(CLASSES)
ncol = 4 if ncl > 6 else 3
nrow = (ncl + ncol - 1) // ncol
fig = plt.figure(figsize=(4.5 * ncol, 5.0 * nrow + 0.9))
outer = fig.add_gridspec(nrow, ncol, hspace=0.30, wspace=0.26)
table, handles, labels = [], [], []
for i, cls in enumerate(CLASSES):
    row, col = divmod(i, ncol)
    inner = outer[row, col].subgridspec(2, 1, height_ratios=[2.3, 1], hspace=0.0)
    ax = fig.add_subplot(inner[0]); rax = fig.add_subplot(inner[1], sharex=ax)
    if MODE == "inclusive":
        dv, de, edges = data(cls[0])
    else:
        dv, de, edges = data((cls[0], cls[2]))
    lab = cls[1]
    c = 0.5 * (edges[1:] + edges[:-1])
    row_txt = [f"{lab:18s} data " + " ".join(f"{x:.3f}" for x in dv)]
    for stem, slab, colr, ls in SAMPLES:
        r = mc(stem, cls[2] if MODE == "inclusive" else cls[0], edges)
        if r is None:
            continue
        v, e = r
        ok = np.isfinite(v) & (v > 0) & (e / np.where(v > 0, v, 1) < MAXREL)
        vv, ee = np.where(ok, v, np.nan), np.where(ok, e, np.nan)
        hatch = "////" if ls == "--" else None
        band(ax, edges, vv, ee, colr, hatch)
        ax.plot(c, vv, ls="none", marker=MARK[slab], ms=4.5, color=colr, mfc=MFC(slab, colr), mew=1.2)
        with np.errstate(divide="ignore", invalid="ignore"):
            band(rax, edges, vv / dv, ee / dv, colr, hatch)
            rax.plot(c, vv / dv, ls="none", marker=MARK[slab], ms=4, color=colr, mfc=MFC(slab, colr), mew=1.1)
        if i == 0:
            import matplotlib.patches as mpatches
            from matplotlib.lines import Line2D
            patch = (mpatches.Patch(facecolor="none", edgecolor=colr, hatch=hatch) if hatch else
                     mpatches.Patch(color=colr, alpha=0.30))
            handles.append((patch, Line2D([], [], ls="none", marker=MARK[slab], ms=5, color=colr, mfc=MFC(slab, colr), mew=1.2)))
            labels.append(slab)
        row_txt.append(f"{'':18s} {slab:24s} " + " ".join("  -  " if not o else f"{x:.3f}" for x, o in zip(v, ok)))
    h = ax.errorbar(c, dv, yerr=de, xerr=np.diff(edges) / 2, fmt="o", color="black", ms=5, lw=1.2, zorder=10)
    if i == 0:
        handles.insert(0, h); labels.insert(0, "Data p–Pb (stat.)")
    # data stat. uncertainty around 1 in the ratio pad
    rel = de / dv
    rax.fill_between(edges, np.append(1 - rel, 1 - rel[-1]), np.append(1 + rel, 1 + rel[-1]), step="post",
                     color="0.55", alpha=0.35, lw=0, zorder=0)
    table += row_txt
    ax.text(0.04, 0.95, lab, transform=ax.transAxes, va="top", fontsize=10.5, weight="bold")
    rax.axhline(1, color="black", lw=0.8)
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM); rax.set_ylim(0.2, 1.45)
    ax.tick_params(labelbottom=False)
    for a in (ax, rax):
        style(a)
    ax.set_ylabel(r"$(\Lambda+\bar{\Lambda})/2K^0_S$"); rax.set_ylabel("MC / data")
    rax.set_xlabel(r"$p_T^{V0}$ (GeV)")
import matplotlib.patches as mpatches
handles.append(mpatches.Patch(color="0.55", alpha=0.35)); labels.append("data stat. unc. (ratio pad)")
fig.legend(handles, labels, loc="upper center", ncol=4, fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.985))
fig.suptitle(TITLE, fontsize=10.5, x=0.01, ha="left", y=1.0)
fig.text(0.99, 0.003, f"N_trk^offline = {NTRK_SCALE} × gen N_ch" + (" (official CMS sample: its own N_trk^offline)" if SET == "pp" else "") + f"; MC bands = stat. unc.; MC bins with rel. error > {MAXREL:.0%} not shown",
         ha="right", fontsize=7.5, color="0.45")
fig.subplots_adjust(top=1 - 1.0 / (5.0 * nrow + 0.9))
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"{NAME}.{ext}"), dpi=130, bbox_inches="tight")
with open(os.path.join(OUT, f"{NAME}.txt"), "w") as f:
    f.write(f"pT edges: {' '.join(f'{x:g}' for x in edges)}\n" + "\n".join(table) + "\n")
print("->", os.path.join(OUT, NAME + ".png"))
