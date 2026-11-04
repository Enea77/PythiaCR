#!/usr/bin/env python3
"""jpct_edges.py -- generated-N_ch edges of the pp in-cone classes with the same fraction of jets as in the data.

  cd PythiaCR && python3 macros/jpct_edges.py dijet_cuetp8m1 dijet_crblc_mode0 dijet_crblc_mode2 dijet_crblc_mode3
For each sample: jets with 60 < pT < 8160 GeV and |eta_WTA| < 2 (forest WTA axis; the V0-free axis differs only slightly)
per generated N_ch (hiBin), weighted with the event weight, from the files of condor/corr_list_pp_mult.txt. The four edges
are chosen from the high-multiplicity end so that the cumulative jet fractions match the data.
Writes condor/multmap_<sample>.txt ("n1,n2,n3,n4", read by condor/corr_job.sh for <sample>_jpct) and
condor/multmap_<sample>_table.txt.
Data (MB trigger, v0-jet-analysis/correlations/files_v0shift/Data_K0s_MB.root, hist_jet_from_reco_reco_sig, jet pT > 60):
4.23M jets; N_trk^offline 0-60: 15.254%, 60-80: 17.462%, 80-185: 63.457%, 185-250: 3.615%, >= 250: 0.212%.
"""
import os, sys
import ROOT
ROOT.gROOT.SetBatch(True); ROOT.EnableImplicitMT()
HERE = os.path.dirname(os.path.abspath(__file__)); BASE = os.path.join(HERE, "..")
DATA = [("0-60", 0.15254), ("60-80", 0.17462), ("80-185", 0.63457), ("185-250", 0.03615), ("250+", 0.00212)]
ROOT.gInterpreter.Declare("""
#include "ROOT/RVec.hxx"
int nJets(const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, int n) {
  int k = 0; for (int i = 0; i < n && i < (int)pt.size(); ++i) if (pt[i] > 60 && pt[i] < 8160 && std::fabs(eta[i]) < 2) ++k; return k; }
""")
lists = {}
for line in open(os.path.join(BASE, "condor", "corr_list_pp_mult.txt")):
    s, files, _ = line.split()
    lists.setdefault(s.replace("_mult", ""), []).extend(files.split("+"))
for s in sys.argv[1:]:
    files = lists[s]
    ev, jt = ROOT.TChain("hiEvtAnalyzer/HiTree"), ROOT.TChain("ak4PFJetAnalyzer/t")
    for f in files: ev.Add(f); jt.Add(f)
    ev.AddFriend(jt)
    df = ROOT.RDataFrame(ev).Define("nj", "nJets(t.genpt, t.WTAgeneta, t.ngen)").Define("wj", "weight * nj")
    h = df.Histo1D(("h", ";gen N_{ch};jets", 500, 0, 500), "hiBin", "wj").GetValue()
    tot = h.Integral(1, 500)
    above = lambda b: h.Integral(b, 500) / tot          # fraction of jets with N_ch >= low edge of bin b
    edges, cum, b_hi = [], 0.0, 500
    for i in range(len(DATA) - 1, 0, -1):               # classes from the top: 250+, 185-250, 80-185, 60-80
        cum += DATA[i][1]
        b = b_hi
        while b > 1 and above(b) < cum: b -= 1
        if b < b_hi and abs(above(b + 1) - cum) < abs(above(b) - cum): b += 1
        edges.insert(0, int(h.GetBinLowEdge(b))); b_hi = b - 1
    e = [0] + edges + [500]
    with open(os.path.join(BASE, "condor", "multmap_%s.txt" % s), "w") as fo: fo.write(",".join(map(str, edges)))
    rows = ["# %s: %d files, %.4g weighted jets (pT > 60, |eta| < 2); gen N_ch edges %s" % (s, len(files), tot, ",".join(map(str, edges))),
            "# data class   data %    pp gen N_ch    pp %"]
    for i, (lbl, fr) in enumerate(DATA):
        b1, b2 = h.FindBin(e[i] + 0.5), h.FindBin(e[i + 1] - 0.5)
        rows.append("%-12s %8.3f   %4d-%-4d   %8.3f" % (lbl, 100 * fr, e[i], e[i + 1], 100 * h.Integral(b1, b2) / tot))
    open(os.path.join(BASE, "condor", "multmap_%s_table.txt" % s), "w").write("\n".join(rows) + "\n")
    print("\n".join(rows))
