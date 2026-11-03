#!/usr/bin/env python3
"""Generator-level overview: official CMS PYTHIA (pgoing, "nominal") vs the standalone samples (pp, 8.16 TeV boosted).

  python3 compare/plot_overview.py          # genCompare outputs compare/<sample>.root  -> compare/plots_overview/
  python3 compare/plot_overview.py v2       # compare/v2/<sample>.root (adds unweighted pTHat and all-event
                                            # histograms)                               -> compare/plots_overview_v2/
Canvases: pthat (unweighted | weighted), jets (pT | eta), charged (Nch | eta | pT), V0 (K0S pT | Lambda pT | Lambda/2K0S),
each with the ratio to the official sample below. Samples whose file does not exist are skipped.

Nominal = sum of the 11 official pTHat-window files (compare[/v2]/official/pthat*.root) with the pTHat-120 window
rescaled by SIG120_FIX: the McM cross section of that sample is ~3.5% high (overlaps with the pTHat-80 and -170 samples
give sigma_rel(120) = 2.11e-4 instead of 2.1849e-4, LOGBOOK section 3).
Bins: consecutive bins are merged (same edges for all samples) until every sample has a relative stat. error <= thr.
Normalisation: weighted distributions per event with pTHat > 15 (official: sigma_rel weights, ours: bias weight / sum of
weights), jet-event distributions per event with a jet > 60 GeV; "jet-V0 pairs" (existing files) per jet.
"""
import glob, os, sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V2 = len(sys.argv) > 1 and sys.argv[1] == "v2"
HERE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(HERE, "v2") if V2 else HERE
OUT = os.path.join(HERE, "plots_overview_v2" if V2 else "plots_overview")
os.makedirs(OUT, exist_ok=True)
SIG120_FIX = 2.11e-4 / 2.1848571e-4

REF = "official"
REFLAB = "CMS official (CUETP8M1, PYTHIA 8.212)"
# CP5 family removed 2026-11-03 (different PDF/alpha_s: harder pTHat, not comparable to the official CUETP8M1)
# (file stem, legend, colour); the first is the reference (nominal); MARKERS: marker per sample
MARKERS = {"dijet_cuetp8m1": "s", "dijet_cp5": "^", "dijet_cp5cr1": "v", "dijet_cp5cr2": "<",
           "dijet_crblc_mode0": "D", "dijet_crblc_mode2": "P", "dijet_crblc_mode3": "X", "dijet_cuetp8m1_py8212": "v"}
SAMPLES = [
    (REF, REFLAB, "black"),
    ("dijet_cuetp8m1", "CUETP8M1", "#2a78d6"),
    ("dijet_crblc_mode0", "CR-BLC mode 0", "#e87ba4"),
    ("dijet_crblc_mode2", "CR-BLC mode 2", "#eb6834"),
    ("dijet_crblc_mode3", "CR-BLC mode 3", "#e34948"),
    ("dijet_cuetp8m1_py8212", "CUETP8M1, PYTHIA 8.212 (our copy)", "#6250d6"),   # jets canvas only
    # our p-Pb Angantyr dijet samples (same tune colour, dashed, filled markers); not in the charged canvases:
    # p-Pb files store gen particles only within dR < 0.8 of a jet
    ("dijetppb_cuetp8m1", "p–Pb CUETP8M1 (Angantyr)", "#2a78d6"),
    ("dijetppb_crblc_mode2", "p–Pb CR-BLC mode 2 (Angantyr)", "#eb6834"),
    ("dijetppb_crblc_mode3", "p–Pb CR-BLC mode 3 (Angantyr)", "#e34948"),
]
PPB = {"dijetppb_cuetp8m1", "dijetppb_crblc_mode2", "dijetppb_crblc_mode3"}
CANVASES = {"dijet_cuetp8m1_py8212": {"jets"}, **{s: {"pthat", "jets", "v0"} for s in PPB}}   # default: all canvases
MARKERS.update({"dijetppb_cuetp8m1": "s", "dijetppb_crblc_mode2": "P", "dijetppb_crblc_mode3": "X"})
files = {}   # stem -> list of (uproot file, scale)
wins = sorted(glob.glob(os.path.join(IN, "official", "pthat*.root")))
assert len(wins) == 11, f"expected 11 official window files in {IN}/official, found {len(wins)}"
files[REF] = [(uproot.open(p), SIG120_FIX if os.path.basename(p) == "pthat120.root" else 1.) for p in wins]
for stem, lab, col in SAMPLES[1:]:
    p = os.path.join(IN, stem + ".root")
    if os.path.exists(p):
        files[stem] = [(uproot.open(p), 1.)]
    else:
        print("skip (no file):", p)
STEMS_ALL = [s for s, _, _ in SAMPLES if s in files]
STEMS = []


def use(canvas):   # samples drawn on a canvas
    global STEMS
    STEMS = [s for s in STEMS_ALL if canvas in CANVASES.get(s, {canvas})]


def raw(stem, name):
    """(sum of weights, sum of w^2, edges) of histogram `name`, summed over the files of the sample."""
    v = s2 = edges = None
    for f, sc in files[stem]:
        if name not in f:
            return None
        h = f[name]
        if v is None:
            v, s2, edges = np.zeros(h.axis().edges().size - 1), np.zeros(h.axis().edges().size - 1), h.axis().edges()
        v += sc * h.values(flow=False)
        s2 += sc * sc * h.variances(flow=False)
    return v, s2, edges


def raw_unweighted(stem, name):   # hPthatRaw: plain counts, no window scale
    out = None
    for f, _ in files[stem]:
        if name not in f:
            return None
        h = f[name]
        v = h.values(flow=False)
        out = (v.copy(), v.copy(), h.axis().edges()) if out is None else (out[0] + v, out[1] + v, out[2])
    return out


def merged_edges(names, thr, getter=raw, xlim=None, maxw=None, maxabs=None):
    """Common edges: merge consecutive bins until every sample (and every histogram in `names`) has rel. error <= thr,
    or (maxw) the merged bin is wider than maxw x its lower edge (single large-weight events must not swallow a range)."""
    # bins from the pp samples only: the p-Pb samples carry a few very large weights (low-pTHat jet events, LOGBOOK 3d)
    # and would merge everything; they are drawn on the same bins with their larger errors
    data = [getter(s, n) for s in STEMS if s not in PPB for n in names]
    data = [d for d in data if d is not None]
    if not data:
        return None
    edges = data[0][2]
    out = [edges[0]]
    acc_v = np.zeros(len(data)); acc_s = np.zeros(len(data))
    for i in range(edges.size - 1):
        acc_v += [d[0][i] for d in data]; acc_s += [d[1][i] for d in data]
        inside = xlim is None or (edges[i + 1] > xlim[0] and edges[i] < xlim[1])
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.sqrt(acc_s) / np.abs(acc_v)
        wide = (maxw is not None and out[-1] > 0 and (edges[i + 1] - out[-1]) >= maxw * out[-1]) or \
               (maxabs is not None and (edges[i + 1] - out[-1]) >= maxabs)
        at_edge = xlim is not None and edges[i + 1] >= xlim[1]   # never merge across the plotted range (official V0s stop at 10 GeV)
        if not inside or wide or at_edge or (np.all(acc_v > 0) and np.all(rel <= thr)):
            out.append(edges[i + 1]); acc_v[:] = 0; acc_s[:] = 0
    if out[-1] != edges[-1]:      # leftover tail: merge into the last bin
        out[-1] = edges[-1]
    return np.array(out)


def rebin(d, new):
    v, s2, e = d
    idx = np.searchsorted(e, new)
    return (np.add.reduceat(v, idx[:-1]), np.add.reduceat(s2, idx[:-1]), new)


def norm_events(stem):   # weighted events with pTHat > 15
    return raw(stem, "hPthat")[0].sum()


def norm_jetevents(stem):
    return raw(stem, "hNjetEvt")[0].sum()


def norm_jets(stem):
    return raw(stem, "hJetEta")[0].sum()


def draw(ax, rax, series, xlabel, ylabel, logx=False, logy=False, xlim=None, rlim=(0.9, 1.1), title=None, empty=""):
    """series: list of (stem, values, errors, edges) with the reference first."""
    ref = None
    for stem, v, e, edges in series:
        lab, col = {s: (l, c) for s, l, c in SAMPLES}[stem]
        c = 0.5 * (edges[1:] + edges[:-1])
        with np.errstate(divide="ignore", invalid="ignore"):
            if stem == REF:
                ref = (v, e)
                ax.errorbar(c, v, yerr=e, xerr=np.diff(edges) / 2, fmt="o", ms=3.5, color="black", lw=1, label=lab, zorder=10)
                re = e / v
                ok = np.isfinite(re)
                rax.bar(edges[:-1][ok], 2 * re[ok], width=np.diff(edges)[ok], bottom=1 - re[ok], align="edge",
                        color="0.82", lw=0, label="CMS official, stat. unc.")
                continue
            mk = MARKERS.get(stem, "o")
            ls, mfc = ("--", col) if stem in PPB else ("-", "white")
            ax.stairs(v, edges, color=col, lw=1.4, ls=ls, baseline=None)
            ax.plot(c, v, ls="none", marker=mk, ms=4.5, color=col, mfc=mfc, mew=1.2)
            ax.plot([], [], ls=ls, marker=mk, ms=4.5, color=col, mfc=mfc, mew=1.2, label=lab)   # legend proxy
            if ref is not None:
                rr = v / ref[0]
                ax_e = rr * e / v   # MC stat. only; the official stat. is the grey band
                rax.stairs(rr, edges, color=col, lw=1.4, ls=ls, baseline=None)
                rax.errorbar(c, rr, yerr=ax_e, fmt=mk, ms=4, color=col, mfc=mfc, mew=1.1, lw=0.9)
    if not series:
        ax.text(0.5, 0.5, empty, transform=ax.transAxes, ha="center", va="center", color="0.4")
    rax.axhline(1, color="black", lw=0.8)
    for a in (ax, rax):
        if logx:
            a.set_xscale("log")
        if xlim:
            a.set_xlim(*xlim)
        a.grid(color="0.92", lw=0.6)
        a.minorticks_on()
        a.tick_params(which="both", direction="in", top=True, right=True)
        a.tick_params(which="major", length=6)
        a.tick_params(which="minor", length=3)
    if logy:
        ax.set_yscale("log")
    rax.set_ylim(*rlim)
    ax.set_ylabel(ylabel)
    rax.set_ylabel("ratio to\nofficial")
    rax.set_xlabel(xlabel)
    ax.tick_params(labelbottom=False)
    if title:
        ax.set_title(title, fontsize=10, loc="left")


def hist_panel(ax, rax, name, norm, thr, xlabel, ylabel, getter=raw, density=True, maxabs=None, ymin=None, **kw):
    edges = merged_edges([name], thr, getter, kw.get("xlim"), 0.25 if kw.get("logx") else None, maxabs)
    series = []
    if edges is not None:
        for s in STEMS:
            d = getter(s, name)
            if d is None:
                continue
            v, s2, e = rebin(d, edges)
            n = norm(s) if norm else 1.
            w = np.diff(e) if density else 1.
            series.append((s, v / n / w, np.sqrt(s2) / n / w, e))
    draw(ax, rax, series, xlabel, ylabel, empty=f"{name}: not in these files\n(v2 histograms on condor)", **kw)
    if ymin is not None:
        ax.set_ylim(bottom=ymin)


def canvas(ncol, width=5.4):
    fig = plt.figure(figsize=(width * ncol, 6.6))
    gs = fig.add_gridspec(2, ncol, height_ratios=[2.6, 1], hspace=0.04, wspace=0.28)
    return fig, [(fig.add_subplot(gs[0, i]), fig.add_subplot(gs[1, i])) for i in range(ncol)]


def finish(fig, name, header):
    fig.suptitle(header, fontsize=10.5, x=0.01, ha="left", y=0.995)
    seen = {}
    for a in fig.axes:   # first axes that has each label (an empty pad must not drop the legend)
        for hh, ll in zip(*a.get_legend_handles_labels()):
            seen.setdefault(ll, hh)
    order = [lab for _, lab, _ in SAMPLES if lab in seen] + [l for l in seen if l not in [s[1] for s in SAMPLES]]
    fig.legend([seen[l] for l in order], order, loc="upper center", ncol=5, fontsize=8.5, frameon=False,
               bbox_to_anchor=(0.5, 0.965))
    fig.text(0.99, 0.005, f"official pTHat-120 sample rescaled by {SIG120_FIX:.4f} (sigma from overlaps); "
             "bins merged to <= thr stat. error", ha="right", fontsize=7, color="0.45")
    fig.subplots_adjust(top=0.86)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("->", os.path.join(OUT, name + ".png"))


HDR = "PYTHIA 8 pp (no nucleus), 8.16 TeV with the p-going boost, generator level"

HDR2 = "PYTHIA 8, generator level: pp (no nucleus; 8.16 TeV, p-going boost) and p–Pb (Angantyr, dashed)"

# 1. pTHat: unweighted (area) | weighted (per event with pTHat > 15)
use("pthat")
fig, ((a1, r1), (a2, r2)) = canvas(2)
hist_panel(a1, r1, "hPthatRaw", lambda s: raw_unweighted(s, "hPthatRaw")[0].sum(), 0.02, r"$\hat{p}_T$ (GeV)",
           r"unweighted, area-normalised  $dN/d\hat{p}_T$", getter=raw_unweighted, logx=True, logy=True,
           xlim=(15, 1500), rlim=(0, 3), title="unweighted (generated events)")
hist_panel(a2, r2, "hPthat", norm_events, 0.01, r"$\hat{p}_T$ (GeV)", r"$(1/N_\mathrm{ev})\,dN/d\hat{p}_T$ (GeV$^{-1}$)",
           logx=True, logy=True, xlim=(15, 1500), rlim=(0.95, 1.05), title="weighted (physical spectrum)")
finish(fig, "pthat", HDR2 + r", all events ($\hat{p}_T$ > 15 GeV)")

# 2. jets: pT (all |eta| < 2 jets, per event) | eta (jets > 60) | pTHat of the events with a jet > 60
#    + our PYTHIA 8.212 build of the same cards: the 2-3% jet deficit of 8.313 is the PYTHIA version
use("jets")
fig, ((a1, r1), (a2, r2), (a3, r3)) = canvas(3, width=4.9)
hist_panel(a1, r1, "hJetPt", norm_events, 0.015, r"gen jet $p_T$ (GeV), anti-$k_T$ 0.4, $|\eta_\mathrm{WTA}|$ < 2",
           r"$(1/N_\mathrm{ev})\,dN_\mathrm{jet}/dp_T$ (GeV$^{-1}$)", logx=True, logy=True, xlim=(30, 1200), rlim=(0.9, 1.3))
hist_panel(a2, r2, "hJetEta", norm_events, 0.01, r"$\eta_\mathrm{WTA}$ (jets $p_T$ > 60 GeV)",
           r"$(1/N_\mathrm{ev})\,dN_\mathrm{jet}/d\eta$", rlim=(0.9, 1.25))
hist_panel(a3, r3, "hPthatJet", norm_events, 0.015, r"$\hat{p}_T$ (GeV) of events with a jet > 60 GeV",
           r"$(1/N_\mathrm{ev})\,dN/d\hat{p}_T$ (GeV$^{-1}$)", logx=True, logy=True, xlim=(15, 1500), rlim=(0.8, 1.1))
finish(fig, "jets", HDR2 + r", per event with $\hat{p}_T$ > 15 GeV")
use("charged")

# 3. charged particles: jet events (existing) and all events (v2)
for suffix, norm, sel in [("", norm_jetevents, "events with a jet > 60 GeV"), ("All", norm_events, r"all events ($\hat{p}_T$ > 15)")]:
    if suffix == "All" and raw(REF, "hNchAll") is None:
        continue
    fig, ((a1, r1), (a2, r2), (a3, r3)) = canvas(3, width=4.8)
    hist_panel(a1, r1, "hNch" + suffix, norm, 0.03, r"$N_\mathrm{ch}$ ($|\eta|$ < 2.4, $p_T$ > 0.4 GeV)",
               r"$(1/N_\mathrm{ev})\,dN_\mathrm{ev}/dN_\mathrm{ch}$", maxabs=10, ymin=1e-7, logy=True, xlim=(0, 200), rlim=(0.5, 2.0))
    y = 0.46
    for st, lab, col in SAMPLES:
        d = raw(st, "hNch" + suffix) if st in STEMS else None
        if d is None:
            continue
        v, _, e = d
        mean = (v * (e[:-1] + 0.5)).sum() / v.sum()   # N_ch is an integer: bin [lo, lo+2) holds lo, lo+1
        a1.plot([0.10], [y], transform=a1.transAxes, marker="o" if st == REF else MARKERS.get(st, "o"), ms=5,
                color=col, mfc=col if st == REF else "white", mew=1.2, ls="none")
        a1.text(0.135, y, f"{lab.split(' (')[0]}: " + r"$\langle N_\mathrm{ch}\rangle$" + f" = {mean:.1f}",
                transform=a1.transAxes, va="center", fontsize=7.8)
        y -= 0.055
    hist_panel(a2, r2, "hChEta" + suffix, norm, 0.01, r"$\eta$ (charged, $p_T$ > 0.4 GeV)",
               r"$(1/N_\mathrm{ev})\,dN_\mathrm{ch}/d\eta$", rlim=(0.9, 1.15))
    hist_panel(a3, r3, "hChPt" + suffix, norm, 0.01, r"$p_T$ (charged, $|\eta|$ < 2.4) (GeV)",
               r"$(1/N_\mathrm{ev})\,dN_\mathrm{ch}/dp_T$ (GeV$^{-1}$)", logy=True, rlim=(0.85, 1.15))
    finish(fig, "charged" + ("_all" if suffix else "_jetevents"), HDR + ", " + sel)


# 4. V0 pT spectra and Lambda/2K0S. Existing files: jet-V0 pairs (any dr) per jet; v2: per event.
def v0sets():
    if raw(REF, "hV0PtJetEvt_Lam") is not None:
        return [("hV0PtJetEvt_", norm_jetevents, "events with a jet > 60 GeV, per event", "v0_jetevents"),
                ("hV0PtAll_", norm_events, r"all events ($\hat{p}_T$ > 15), per event", "v0_all")]
    return [("hEvt_", norm_jets, "events with a jet > 60 GeV, all V0s of the event per jet (jet-V0 pairs)", "v0_jetevents")]


use("v0")
for pre, norm, sel, fname in v0sets():
    fig, ((a1, r1), (a2, r2), (a3, r3)) = canvas(3, width=4.8)
    hist_panel(a1, r1, pre + "K0s", norm, 0.01, r"$p_T$ (GeV)", r"$K^0_S$: $dN/dp_T$ (GeV$^{-1}$), $|\eta|$ < 2",
               logy=True, xlim=(1, 10), rlim=(0.8, 1.2))
    hist_panel(a2, r2, pre + "Lam", norm, 0.01, r"$p_T$ (GeV)", r"$\Lambda+\bar{\Lambda}$: $dN/dp_T$ (GeV$^{-1}$), $|\eta|$ < 2",
               logy=True, xlim=(1, 10), rlim=(0.8, 1.3))
    edges = merged_edges([pre + "Lam", pre + "K0s"], 0.01, xlim=(1, 10))
    series = []
    for s in STEMS:
        l, k = raw(s, pre + "Lam"), raw(s, pre + "K0s")
        if l is None or k is None:
            continue
        l, k = rebin(l, edges), rebin(k, edges)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = 0.5 * l[0] / k[0]
            series.append((s, r, r * np.sqrt(l[1] / l[0] ** 2 + k[1] / k[0] ** 2), edges))
    draw(a3, r3, series, r"$p_T$ (GeV)", r"$(\Lambda+\bar{\Lambda})/2K^0_S$", xlim=(1, 10), rlim=(0.7, 1.5))
    a3.set_ylim(0, 0.45)
    finish(fig, fname, HDR2 + ", " + sel)
