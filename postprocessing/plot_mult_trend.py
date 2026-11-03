#!/usr/bin/env python3
"""In-cone (L+AL)/2K0S (dr < 0.4, 3-10 GeV) vs N_trk^offline class: data vs p-Pb tunes.

Reads the values.txt files written by run_tunes_vs_data.sh.
  python3 plot_mult_trend.py [variant]      # variant "" or "_pthat30" (default "_pthat30")
Output: plots/tunes_vs_data_mult<variant>/MultTrend_InCone.{pdf,png} and MultTrend.txt
"""
import re, sys, os
import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V = sys.argv[1] if len(sys.argv) > 1 else "_pthat30"
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"plots/tunes_vs_data_mult{V}")
classes = [("mult0to60", "0–60"), ("mult60to80", "60–80"), ("mult80to185", "80–185"),
           ("mult185to250", "185–250\n(HM185)"), ("mult250to1000", "≥ 250\n(HM250)")]
# MC points with a relative error above this are not drawn (HM250: 40-100%)
MAXREL = 0.25

# label in values.txt -> (legend, color, marker, filled, x offset)
series = {
    "pPb CUETP8M1":         ("CUETP8M1",             "#2a78d6", "s", False, -0.30),
    "pPb CR-BLC 2":         ("CR-BLC mode 2",        "#eb6834", "P", False, -0.18),
    "pPb CR-BLC 3":         ("CR-BLC mode 3",        "#e34948", "X", False,  0.12),
    "pPb central CUETP8M1": ("CUETP8M1, central",    "#2a78d6", "s", True,   0.24),
    "pPb central CR-BLC 2": ("CR-BLC mode 2, central", "#eb6834", "P", True,  0.36),
}
rx = re.compile(r"^(.*?)\s+3-10 GeV:\s+([-\d.]+)\s+\+-\s+([\d.]+)")

vals = {}  # label -> list of (val, err) per class
for ic, (d, _) in enumerate(classes):
    with open(os.path.join(base, d, "values.txt")) as f:
        for line in f:
            m = rx.match(line)
            if m:
                lab = "Data" if m.group(1).startswith("Data") else m.group(1).strip()
                vals.setdefault(lab, [None] * len(classes))[ic] = (float(m.group(2)), float(m.group(3)))

fig, (ax, rx_) = plt.subplots(2, 1, figsize=(7, 6.5), sharex=True,
                              gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.05})
x = list(range(len(classes)))
dv = vals["Data"]
ax.errorbar(x, [v for v, _ in dv], yerr=[e for _, e in dv], fmt="o", color="black", ms=8,
            lw=1.5, capsize=3, label="Data p–Pb (stat.)", zorder=5)
rx_.axhline(1, color="black", lw=1)
rx_.fill_between([-0.5, len(classes) - 0.5], 0.95, 1.05, color="0.85", lw=0, label="±5%")

out = ["class  " + "  ".join(f"{c[0]:>14}" for c in classes)]
out.append("Data   " + "  ".join(f"{v:.4f}±{e:.4f}" for v, e in dv))
for lab, (leg, col, mk, filled, dx) in series.items():
    if lab not in vals:
        continue
    pts = [(i, v, e) for i, ve in enumerate(vals[lab]) if ve for v, e in [ve] if v > 0 and e / v < MAXREL]
    kw = dict(fmt=mk, color=col, ms=7, lw=1.5, capsize=2,
              mfc=col if filled else "white", mew=1.5)
    ax.errorbar([i + dx for i, _, _ in pts], [v for _, v, _ in pts], yerr=[e for _, _, e in pts],
                label=leg, **kw)
    rx_.errorbar([i + dx for i, _, _ in pts], [v / dv[i][0] for i, v, _ in pts],
                 yerr=[v / dv[i][0] * ((e / v) ** 2 + (dv[i][1] / dv[i][0]) ** 2) ** 0.5 for i, v, e in pts], **kw)
    out.append(f"{leg:22s} " + "  ".join(
        (f"{ve[0]:.4f}±{ve[1]:.4f}" if ve else "-") for ve in vals[lab]))

ax.set_ylabel(r"In-cone $(\Lambda+\bar{\Lambda})/2K^0_S$")
ax.set_ylim(0.15, 0.345)
ax.text(0.03, 0.95, "p–Pb 8.16 TeV, jets $p_T$ > 60 GeV, $\\Delta r$ < 0.4, 3 < $p_T^{V0}$ < 10 GeV\n"
        "PYTHIA 8.313 Angantyr, $\\hat{p}_T$ > 30 GeV" if V == "_pthat30" else
        "p–Pb 8.16 TeV, jets $p_T$ > 60 GeV, $\\Delta r$ < 0.4, 3 < $p_T^{V0}$ < 10 GeV\nPYTHIA 8.313 Angantyr",
        transform=ax.transAxes, va="top", fontsize=9, color="0.25")
ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.84), fontsize=8.5, frameon=False, ncol=2)
ax.grid(axis="y", color="0.92")
rx_.set_ylabel("MC / data")
rx_.set_ylim(0.45, 1.15)
rx_.grid(axis="y", color="0.92")
rx_.set_xticks(x, [c[1] for c in classes])
rx_.set_xlim(-0.5, len(classes) - 0.5)
rx_.set_xlabel(r"$N_\mathrm{trk}^\mathrm{offline}$ class")
rx_.text(len(classes) - 0.55, 0.48, f"MC points with rel. error > {MAXREL:.0%} not shown",
         ha="right", fontsize=7.5, color="0.4")
for a in (ax, rx_):
    a.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    a.tick_params(which="both", direction="in", top=True, right=True)
    a.tick_params(which="major", length=6)
    a.tick_params(which="minor", length=3)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(base, f"MultTrend_InCone.{ext}"), bbox_inches="tight", dpi=150)
with open(os.path.join(base, "MultTrend.txt"), "w") as f:
    f.write("\n".join(out) + "\n")
print("\n".join(out))
print("->", os.path.join(base, "MultTrend_InCone.png"))
