// plot_incone.C -- in-cone (Delta R < 0.4) V0 yields per jet and (Lambda+anti-Lambda)/2K0S vs V0 pT in multiplicity
// classes, from the data outputs of jetshape.cpp (run_jetshapes_data.C), counts (not pT-weighted).
// Version of JetAnalyzer/Plot_InCone.cpp for the new chain: classes, colours and markers as inclusive/plot_inclusive.C,
// systematics from inclusive/incl_systematics.C (in-jet totals per pT range), right-hand panel with the 3-10 GeV bin
// (run_jetshapes_data.C with ptEdges = "3,10").
//
//   root -l -b -q 'plot_incone.C+("../correlations/files/","0,60,80,185")'
// (jets above 60 GeV sit at high N_trk: the inclusive classes 0,20,30,40,60 have too few jets below 60)
//
// stage: correction stage of jetshape.cpp, default "sideband_corrected" (Delta eta UE subtracted, no seagull; nominal)
// cone:  "innercone" (Delta R < 0.4, default) or "" (Delta R < 1)
// int3to10: "sum" (default) = 3-10 GeV bin as the sum of the per-pT-bin yields (3-5, 5-7, 7-10), consistent with the
//           per-bin mixing normalization and mass sidebands; "run" = dedicated jetshape.cpp run with one 3-10 GeV bin
// Delta r distributions from the jet shapes of the 5 V0 pT bins, summed (1.5-10 GeV): P(Delta r) = dN/dDelta r per jet
// (DeltaR_*) and rho(Delta r) normalized to unit integral for Delta r < 1 (JetShape_*); Lambda, K0S, ratio, double ratio.
// drMerge / drEdges: classes (default "HM250", low statistics) whose Delta r ratio and double ratio use coarser bins
//                    (default "0,0.1,0.2,0.4,1,1.5": 0.4-0.6 and 0.6-1 merged), double ratio to MB merged the same way
//                    (drMerge = "" for the jet-shape binning everywhere)
// Output plots/incone/: InCone_Ratio, InCone_DoubleRatio_MB, InCone_Yield_Lam, InCone_Yield_K0s (per jet per GeV, log),
//                       DeltaR_{Lam,K0s,Ratio,DoubleRatio_MB}, JetShape_{Lam,K0s,Ratio,DoubleRatio_MB} (.png/.pdf),
//                       incone_ratios.root, incone_values.txt
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"
#include "TArrayI.h"
#include "TGraphErrors.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include <vector>
#include "../include/PlotStyle.h"

namespace pic {

int TempColor(double f) {
    static TArrayI pal;
    if (pal.GetSize() == 0) { gStyle->SetPalette(kTemperatureMap); pal = TColor::GetPalette(); gStyle->SetPalette(kBird); }
    return pal[int(f * (pal.GetSize() - 1) + 0.5)];
}

struct Cls { TString name, label, sample; int lo, hi; int color, marker; double msize; };
double MSize(int m) { return m == 33 ? 2.1 : m == 29 ? 2.0 : m == 43 ? 2.2 : m == 47 ? 1.6 : m == 34 ? 1.6 : m == 21 ? 1.4 : 1.5; }   // similar visual size
// as inclusive/plot_inclusive.C: MB black; MB classes and HM on the kTemperatureMap gradient (near-white middle skipped);
// point-symmetric markers
std::vector<Cls> Classes(const std::vector<double>& e) {
    std::vector<Cls> c;
    c.push_back({"MB", "MB integrated", "MB", -1, -1, kBlack, 20, 1.4});
    const int markers[] = {21, 33, 34, 47, 29, 22, 23};
    const int n = e.size() - 1, nAll = n + 2;
    auto col = [&](int k) { const double t = nAll > 1 ? double(k) / (nAll - 1) : 0; return TempColor(t <= 0.5 ? 0.5 * t : 0.5 + 0.5 * t); };   // light middle of the map avoided
    for (int i = 0; i < n; ++i)
        c.push_back({Form("mult%gto%g", e[i], e[i + 1]), Form("MB: %g-%g", e[i], e[i + 1]), "MB", (int)e[i], (int)e[i + 1], col(i), markers[i % 7], MSize(markers[i % 7])});
    c.push_back({"HM185", "HM185: 185-250", "HM185", 185, 250, col(n), 47, 1.6});
    c.push_back({"HM250", "HM250: 250+", "HM250", 250, 1000, col(n + 1), 43, 2.2});
    return c;
}

// jetshape.cpp output name (same rule as in jetshape.cpp)
// in-jet efficiency corrected outputs (tag "ijc_", jetshape.cpp with injetEff) are used when present
TString gIJTag = "ijc_";
TString gPtw = "";                 // "ptw_": pT-weighted correlations (plot_incone(..., ptw = true))
TString gOut = "plots/incone/";    // output directory (plots/incone_ptw/ for the pT-weighted version)
TString JSFile(const TString& dir, const Cls& c, const char* v0, const TString& ptTag) {
    TString m = c.lo > -1 ? TString(Form("mult%dto%d_", c.lo, c.hi)) : TString("");
    const TString f = dir + "jetshapes_rebinned_" + ptTag + gIJTag + gPtw + m + "reco_reco_Data_" + v0 + "_" + c.sample + ".root";
    if (gIJTag != "" && gSystem->AccessPathName(f)) {
        std::cout << "WARNING: no in-jet efficiency corrected output " << f << ", using the uncorrected one" << std::endl;
        return dir + "jetshapes_rebinned_" + ptTag + gPtw + m + "reco_reco_Data_" + v0 + "_" + c.sample + ".root";
    }
    return f;
}
TH1D* GetYield(const TString& file, const TString& cone, const TString& stage, const TString& name) {
    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie()) { std::cout << "missing " << file << std::endl; return nullptr; }
    const TString d = cone == "" ? "Integral_vs_pt" : "Integral_vs_pt_" + cone;
    TH1* h = (TH1*)f->Get(d + "/" + d + "_" + stage);
    TH1D* o = h ? (TH1D*)h->Clone(name) : nullptr;
    if (o) o->SetDirectory(0); else std::cout << "missing " << d << "_" << stage << " in " << file << std::endl;
    f->Close();
    return o;
}

// systematics per pT range (lines "<quantity> <pT low> <pT high> <value>")
struct SysRange { double lo, hi, val; };
std::vector<SysRange> ReadSys(const TString& file, const TString& q) {
    std::vector<SysRange> s;
    std::ifstream in(file.Data()); std::string k; double lo, hi, v;
    while (in >> k) {
        if (k[0] == '#') { std::getline(in, k); continue; }
        in >> lo >> hi >> v;
        if (TString(k) == q) s.push_back({lo, hi, v});
    }
    return s;
}
// largest value over the ranges overlapping [lo, hi)
double SysIn(const std::vector<SysRange>& s, double lo, double hi) {
    double m = 0;
    for (auto& r : s) if (r.lo < hi - 1e-6 && r.hi > lo + 1e-6) m = std::max(m, r.val);
    return m;
}

// points with an x offset (fraction of the bin width) so that the classes do not sit on top of each other
TGraphErrors* Points(TH1D* h, double off, const Cls& c) {
    TGraphErrors* g = new TGraphErrors();
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        if (h->GetBinContent(b) == 0) continue;
        const int n = g->GetN();
        g->SetPoint(n, h->GetBinCenter(b) + off * h->GetBinWidth(b), h->GetBinContent(b));
        g->SetPointError(n, 0, h->GetBinError(b));
    }
    g->SetMarkerStyle(c.marker); g->SetMarkerSize(c.msize); g->SetMarkerColor(c.color); g->SetLineColor(c.color); g->SetLineWidth(2);
    return g;
}
TGraphErrors* Boxes(TH1D* h, double off, double halfw, const std::vector<SysRange>& s, const Cls& c) {
    TGraphErrors* g = new TGraphErrors();
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        if (h->GetBinContent(b) == 0) continue;
        const int n = g->GetN();
        g->SetPoint(n, h->GetBinCenter(b) + off * h->GetBinWidth(b), h->GetBinContent(b));
        g->SetPointError(n, halfw * h->GetBinWidth(b), h->GetBinContent(b) * SysIn(s, h->GetBinLowEdge(b), h->GetXaxis()->GetBinUpEdge(b)));
    }
    g->SetFillColorAlpha(c.color, 0.35); g->SetLineColor(c.color); g->SetLineWidth(0);
    return g;
}
void Extent(const std::vector<TGraphErrors*>& gs, double& lo, double& hi) {
    for (auto g : gs) for (int i = 0; i < g->GetN(); ++i) { lo = std::min(lo, g->GetY()[i] - g->GetEY()[i]); hi = std::max(hi, g->GetY()[i] + g->GetEY()[i]); }
}

// 3-10 GeV yield as the sum of the per-pT-bin yields (each bin processed with its own mixing normalization and mass
// sidebands). jetshape.cpp stores Integral_vs_pt per (width of the first pT bin): yield_i = stored_i * w_i / w_1.
TH1D* SumBins(TH1D* h, double lo, double hi, const TString& name) {
    if (!h) return nullptr;
    TH1D* o = new TH1D(name, "", 1, lo, hi); o->SetDirectory(0);
    double v = 0, e2 = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        if (h->GetBinCenter(b) < lo || h->GetBinCenter(b) > hi) continue;
        const double f = h->GetBinWidth(b) / h->GetBinWidth(1);
        v += h->GetBinContent(b) * f; e2 += std::pow(h->GetBinError(b) * f, 2);
    }
    o->SetBinContent(1, v); o->SetBinError(1, std::sqrt(e2));
    return o;
}
// dN/dDeltaR per jet summed over the V0 pT bins of jetshape.cpp (per-bin jet shapes are yields per jet, not unit-normalized)
TH1D* GetDR(const TString& file, const TString& stage, const TString& name, double& ptlo, double& pthi) {
    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie()) return nullptr;
    TH1D* s = nullptr;
    for (int i = 0;; ++i) {
        TH1* h = (TH1*)f->Get(Form("jetshape/JetShape_%s_c%d", stage.Data(), i));
        if (!h) break;
        if (!s) { s = (TH1D*)h->Clone(name); s->SetDirectory(0); } else s->Add(h);
    }
    if (TH1* p = (TH1*)f->Get("Integral_vs_pt/Integral_vs_pt_raw")) { ptlo = p->GetXaxis()->GetXmin(); pthi = p->GetXaxis()->GetXmax(); }
    f->Close();
    return s;
}

// density histogram (per unit Delta r) on coarser edges (each must be an edge of h): contents and errors summed with
// the bin widths, then divided by the new widths
TH1D* MergeDR(TH1D* h, const std::vector<double>& edges, const TString& name) {
    if (!h) return nullptr;
    TH1D* o = new TH1D(name, h->GetTitle(), edges.size() - 1, edges.data()); o->SetDirectory(0); o->Sumw2();
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const int k = o->FindBin(h->GetBinCenter(b));
        if (k < 1 || k > o->GetNbinsX()) continue;
        const double w = h->GetBinWidth(b), c = h->GetBinContent(b) * w, e = h->GetBinError(b) * w;
        o->SetBinContent(k, o->GetBinContent(k) + c);
        o->SetBinError(k, std::sqrt(o->GetBinError(k) * o->GetBinError(k) + e * e));
    }
    for (int k = 1; k <= o->GetNbinsX(); ++k) { const double w = o->GetBinWidth(k); o->SetBinContent(k, o->GetBinContent(k) / w); o->SetBinError(k, o->GetBinError(k) / w); }
    return o;
}

struct Res { Cls c; TH1D *lam = nullptr, *k0s = nullptr, *r = nullptr, *lamI = nullptr, *k0sI = nullptr, *rI = nullptr,
                       *dLam = nullptr, *dK0s = nullptr, *dR = nullptr; };

// one pad, x = Delta R (variable bins), points offset per class, constant relative systematic (0: none)
void DrawDR(const std::vector<Res>& res, const std::vector<TH1D*>& hs, double sysRel, const char* ytitle, double xhi, bool line1, bool logy,
            const TString& out, const TString& info, TH1D* mc = nullptr) {
    TCanvas* c = new TCanvas("c_" + out, "", 800, 800);
    const double L = 0.17, R = 0.04, T = 0.08, B = 0.13;
    gPad->SetLeftMargin(L); gPad->SetRightMargin(R); gPad->SetTopMargin(T); gPad->SetBottomMargin(B); gPad->SetTickx(1); gPad->SetTicky(1);
    if (logy) gPad->SetLogy();
    const int n = res.size();
    std::vector<TGraphErrors*> pts(n, nullptr), box(n, nullptr);
    for (int i = 0; i < n; ++i) {
        if (!hs[i]) continue;
        const double off = (i - 0.5 * (n - 1)) * 0.07;
        TGraphErrors* g = new TGraphErrors(); TGraphErrors* s = new TGraphErrors();
        for (int b = 1; b <= hs[i]->GetNbinsX(); ++b) {
            if (hs[i]->GetBinLowEdge(b) >= xhi - 1e-6 || hs[i]->GetBinContent(b) == 0) continue;
            if (logy && hs[i]->GetBinContent(b) <= 0) continue;
            // offset scale: finest bin at this Delta r over the classes; a merged (wider) bin is drawn left of the
            // leftmost class in that finest bin, so it does not overlap the others
            double w = hs[i]->GetBinWidth(b), xc = hs[i]->GetBinCenter(b);
            for (TH1D* o : hs) if (o) { const int k = o->FindBin(hs[i]->GetBinCenter(b)); if (o->GetBinWidth(k) < w - 1e-9) { w = o->GetBinWidth(k); xc = o->GetBinCenter(k); } }
            const bool merged = w < hs[i]->GetBinWidth(b) - 1e-9;
            const double x = xc + (merged ? -(0.5 * (n - 1) + 1) * 0.07 : off) * w, y = hs[i]->GetBinContent(b);
            g->SetPoint(g->GetN(), x, y); g->SetPointError(g->GetN() - 1, 0, hs[i]->GetBinError(b));
            s->SetPoint(s->GetN(), x, y); s->SetPointError(s->GetN() - 1, 0.03 * w,std::fabs(y) * sysRel);
        }
        g->SetMarkerStyle(res[i].c.marker); g->SetMarkerSize(res[i].c.msize); g->SetMarkerColor(res[i].c.color); g->SetLineColor(res[i].c.color); g->SetLineWidth(2);
        s->SetFillColorAlpha(res[i].c.color, 0.35); s->SetLineWidth(0);
        pts[i] = g; if (sysRel > 0) box[i] = s;
    }
    std::vector<TGraphErrors*> all;
    for (auto v : {pts, box}) for (auto g : v) if (g) all.push_back(g);
    TGraphErrors* gmc = nullptr;
    if (mc) {
        gmc = new TGraphErrors();
        for (int b = 1; b <= mc->GetNbinsX(); ++b) {
            if (mc->GetBinLowEdge(b) >= xhi - 1e-6 || mc->GetBinContent(b) <= 0) continue;
            const int k = gmc->GetN(); gmc->SetPoint(k, mc->GetBinCenter(b), mc->GetBinContent(b)); gmc->SetPointError(k, 0, mc->GetBinError(b));
        }
        gmc->SetFillColorAlpha(kViolet, 0.4); gmc->SetLineColor(kViolet); gmc->SetLineWidth(2); all.push_back(gmc);
    }
    double lo = 1e30, hi = -1e30; Extent(all, lo, hi);
    const int ncol = 2, nrow = (n + (mc ? 1 : 0) + ncol - 1) / ncol;
    const double legTop = 0.90, legBot = legTop - 0.042 * (nrow + 1);
    const double frac = (legBot - 0.02 - B) / (1 - T - B);
    double ymin, ymax;
    if (logy) {
        double lp = 1e30; for (auto g : all) for (int i = 0; i < g->GetN(); ++i) if (g->GetY()[i] > 0) lp = std::min(lp, g->GetY()[i] - g->GetEY()[i] > 0 ? g->GetY()[i] - g->GetEY()[i] : 0.5 * g->GetY()[i]);
        ymin = lp / 2; ymax = std::exp(std::log(ymin) + (std::log(hi) - std::log(ymin)) / frac);
    } else { ymin = lo - 0.08 * (hi - lo); ymax = ymin + (hi - ymin) / frac; }
    TH1D* fr = new TH1D("fr_" + out, Form(";#Deltar;%s", ytitle), 100, 0, xhi);
    fr->SetMinimum(ymin); fr->SetMaximum(ymax);
    fr->GetYaxis()->SetTitleOffset(1.6); fr->Draw("AXIS");
    if (line1) { TLine* l = new TLine(0, 1, xhi, 1); l->SetLineStyle(7); l->SetLineColor(kGray + 2); l->Draw(); }
    if (gmc) { gmc->Draw("3"); gmc->Draw("LX"); }
    for (auto g : box) if (g) g->Draw("2");
    for (auto g : pts) if (g) g->Draw("PZ");
    TLegend* leg = new TLegend(L + 0.03, legBot, 1 - R - 0.01, legTop); leg->SetNColumns(ncol); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
    leg->SetHeader(info);
    for (int i = 0; i < n; ++i) if (pts[i]) leg->AddEntry(pts[i], res[i].c.label, "p");
    if (gmc) leg->AddEntry(gmc, "PYTHIA 8", "lf");
    leg->Draw();
    PS::Text(L, 0.955, "#font[61]{CMS} #font[52]{Preliminary}", 0.045);
    PS::Text(1 - R, 0.955, "pPb 8.16 TeV", 0.038, 42, 32);
    gPad->RedrawAxis();
    PS::Save(c, gOut + out);
}

// two pads: pT dependence (left) and the 3-10 GeV bin (right), same y axis; legend below the CMS header, never on points
void DrawTwoPanel(const std::vector<Res>& res, const std::vector<TH1D*>& hs, const std::vector<TH1D*>& hsI, const std::vector<SysRange>* sys,
                  const char* ytitle, bool line1, const TString& out, const TString& info, bool logy = false,
                  TH1D* mc = nullptr, TH1D* mcI = nullptr) {
    TCanvas* c = new TCanvas("c_" + out, "", 1000, 800);
    const double split = 0.80;
    TPad* p1 = new TPad("p1_" + out, "", 0, 0, split, 1); TPad* p2 = new TPad("p2_" + out, "", split, 0, 1, 1);
    p1->SetLeftMargin(0.15 / split); p1->SetRightMargin(0.0); p1->SetTopMargin(0.08); p1->SetBottomMargin(0.13);
    p2->SetLeftMargin(0.0); p2->SetRightMargin(0.05 / (1 - split)); p2->SetTopMargin(0.08); p2->SetBottomMargin(0.13);
    for (TPad* p : {p1, p2}) { p->SetTickx(1); p->SetTicky(1); c->cd(); p->Draw(); }
    const int n = res.size();
    std::vector<TGraphErrors*> pts, box, ptsI, boxI;
    for (int i = 0; i < n; ++i) {
        if (!hs[i]) { pts.push_back(nullptr); box.push_back(nullptr); continue; }
        const double off = (i - 0.5 * (n - 1)) * 0.07;   // within +-0.25 of the bin width
        pts.push_back(Points(hs[i], off, res[i].c));
        box.push_back(sys ? Boxes(hs[i], off, 0.03, *sys, res[i].c) : nullptr);
    }
    for (int i = 0; i < n; ++i) {
        ptsI.push_back(nullptr); boxI.push_back(nullptr);
        if (!hsI[i] || hsI[i]->GetBinContent(hsI[i]->GetNbinsX()) == 0) continue;
        const int b = hsI[i]->GetNbinsX();
        const double x = (i + 0.5) / n, v = hsI[i]->GetBinContent(b);
        TGraphErrors* g = new TGraphErrors(1); g->SetPoint(0, x, v); g->SetPointError(0, 0, hsI[i]->GetBinError(b));
        g->SetMarkerStyle(res[i].c.marker); g->SetMarkerSize(res[i].c.msize); g->SetMarkerColor(res[i].c.color); g->SetLineColor(res[i].c.color); g->SetLineWidth(2);
        ptsI.back() = g;
        if (sys) {
            TGraphErrors* s = new TGraphErrors(1); s->SetPoint(0, x, v);
            s->SetPointError(0, 0.3 / n, v * SysIn(*sys, hsI[i]->GetBinLowEdge(b), hsI[i]->GetXaxis()->GetBinUpEdge(b)));
            s->SetFillColorAlpha(res[i].c.color, 0.35); s->SetLineWidth(0);
            boxI.back() = s;
        }
    }
    // y range: everything, then room for the legend (top of the left pad)
    double lo = 1e30, hi = -1e30;
    std::vector<TGraphErrors*> all;
    for (auto v : {pts, box, ptsI, boxI}) for (auto g : v) if (g) all.push_back(g);
    TGraphErrors *gmc = nullptr, *gmcI = nullptr;   // PYTHIA band (pT bins) and its 3-10 GeV value across the right pad
    if (mc) {
        gmc = new TGraphErrors();
        for (int b = 1; b <= mc->GetNbinsX(); ++b) if (mc->GetBinContent(b) > 0) { const int k = gmc->GetN(); gmc->SetPoint(k, mc->GetBinCenter(b), mc->GetBinContent(b)); gmc->SetPointError(k, 0, mc->GetBinError(b)); }
        gmc->SetFillColorAlpha(kViolet, 0.4); gmc->SetLineColor(kViolet); gmc->SetLineWidth(2); all.push_back(gmc);
    }
    if (mcI) {
        const int b = mcI->GetNbinsX();
        gmcI = new TGraphErrors(1); gmcI->SetPoint(0, 0.5, mcI->GetBinContent(b)); gmcI->SetPointError(0, 0.5, mcI->GetBinError(b));
        gmcI->SetFillColorAlpha(kViolet, 0.4); gmcI->SetLineColor(kViolet); gmcI->SetLineWidth(2); all.push_back(gmcI);
    }
    Extent(all, lo, hi);
    const int ncol = 3, nrow = (n + (mc ? 1 : 0) + ncol - 1) / ncol;
    const double legTop = 0.90, legBot = legTop - 0.045 * (nrow + 1);
    const double frac = (legBot - 0.02 - 0.13) / (1 - 0.08 - 0.13);   // free fraction of the frame below the legend
    double ymin, ymax;
    if (logy) {
        double lp = 1e30; for (auto g : all) for (int i = 0; i < g->GetN(); ++i) if (g->GetY()[i] > 0) lp = std::min(lp, g->GetY()[i] - g->GetEY()[i] > 0 ? g->GetY()[i] - g->GetEY()[i] : 0.5 * g->GetY()[i]);
        ymin = lp / 2; ymax = std::exp(std::log(ymin) + (std::log(hi) - std::log(ymin)) / frac);
    } else {
        ymin = lo - 0.08 * (hi - lo); ymax = ymin + (hi - ymin) / frac;
    }
    if (logy) { p1->SetLogy(); p2->SetLogy(); }
    p1->cd();
    TH1D* fr = new TH1D("fr_" + out, Form(";p_{T}^{V0} (GeV);%s", ytitle), 100, 1.0, 10.0);
    fr->SetMinimum(ymin); fr->SetMaximum(ymax);
    fr->GetXaxis()->SetTitleSize(0.05); fr->GetXaxis()->SetLabelSize(0.045);
    fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetTitleOffset(1.45);
    fr->GetXaxis()->ChangeLabel(-1, -1, 0);   // "10" sits on the pad border
    fr->Draw("AXIS");
    if (line1) { TLine* l = new TLine(1.0, 1.0, 10.0, 1.0); l->SetLineStyle(7); l->SetLineColor(kGray + 2); l->Draw(); }
    if (gmc) { gmc->Draw("3"); gmc->Draw("LX"); }   // band (stat.) connecting the bin centres
    for (auto g : box) if (g) g->Draw("2");
    for (auto g : pts) if (g) g->Draw("PZ");
    TLegend* leg = new TLegend(0.20, legBot, 0.99, legTop); leg->SetNColumns(ncol); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
    leg->SetHeader(info);
    for (int i = 0; i < n; ++i) if (pts[i]) leg->AddEntry(pts[i], res[i].c.label, "p");
    if (gmc) leg->AddEntry(gmc, "PYTHIA 8", "lf");
    leg->Draw();
    PS::Text(0.15 / split + 0.005, 0.955, "#font[61]{CMS} #font[52]{Preliminary}", 0.05);
    PS::Text(0.995, 0.955, "pPb 8.16 TeV", 0.042, 42, 32);
    p1->RedrawAxis();
    p2->cd();
    TH1D* f2 = new TH1D("f2_" + out, ";;", 1, 0, 1);
    f2->SetMinimum(ymin); f2->SetMaximum(ymax);
    f2->GetXaxis()->SetBinLabel(1, "3-10 GeV"); f2->GetXaxis()->SetLabelSize(0.05 * split / (1 - split)); f2->GetXaxis()->SetLabelOffset(0.005);
    f2->GetYaxis()->SetLabelSize(0); f2->GetYaxis()->SetTickLength(0.03 * split / (1 - split));
    f2->GetXaxis()->SetTickLength(0);
    f2->Draw("AXIS");
    if (line1) { TLine* l = new TLine(0, 1.0, 1, 1.0); l->SetLineStyle(7); l->SetLineColor(kGray + 2); l->Draw(); }
    if (gmcI) { gmcI->Draw("2"); TLine* l = new TLine(0, gmcI->GetY()[0], 1, gmcI->GetY()[0]); l->SetLineColor(kViolet); l->SetLineWidth(2); l->Draw(); }
    for (auto g : boxI) if (g) g->Draw("2");
    for (auto g : ptsI) if (g) g->Draw("PZ");
    p2->RedrawAxis();
    PS::Save(c, gOut + out);
}

} // namespace pic

void plot_incone(TString dir = "../correlations/files/", TString mbEdges = "0,60,80,185", TString stage = "sideband_corrected",
                 TString cone = "innercone", TString sysFile = "../inclusive/plots/systematics/sys_inclusive.txt", TString int3to10 = "sum",
                 TString mcDir = "input/", TString drMerge = "HM250", TString drEdges = "0,0.1,0.2,0.4,1,1.5", bool ptw = false)
{
    pic::gPtw = ptw ? "ptw_" : ""; pic::gOut = ptw ? "plots/incone_ptw/" : "plots/incone/";
    using namespace pic;
    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    PS::SetGlobal();
    PS::MkDirs(pic::gOut);
    if (!dir.EndsWith("/")) dir += "/";
    std::vector<double> e;
    { TObjArray* t = mbEdges.Tokenize(","); for (int i = 0; i < t->GetEntries(); ++i) e.push_back(((TObjString*)t->At(i))->GetString().Atof()); delete t; }
    const std::vector<Cls> cls = Classes(e);

    std::vector<Res> res;
    for (const auto& c : cls) {
        Res r; r.c = c;
        r.lam  = GetYield(JSFile(dir, c, "Lam", ""), cone, stage, "lam_" + c.name);
        r.k0s  = GetYield(JSFile(dir, c, "K0s", ""), cone, stage, "k0s_" + c.name);
        if (int3to10 == "sum") { r.lamI = SumBins(r.lam, 3, 10, "lamI_" + c.name); r.k0sI = SumBins(r.k0s, 3, 10, "k0sI_" + c.name); }
        else {
            r.lamI = GetYield(JSFile(dir, c, "Lam", "ptbins3to10_"), cone, stage, "lamI_" + c.name);
            r.k0sI = GetYield(JSFile(dir, c, "K0s", "ptbins3to10_"), cone, stage, "k0sI_" + c.name);
        }
        {   // comparison with the dedicated 3-10 GeV processing
            TH1D* a = GetYield(JSFile(dir, c, "Lam", "ptbins3to10_"), cone, stage, "cmpL_" + c.name), *b = GetYield(JSFile(dir, c, "K0s", "ptbins3to10_"), cone, stage, "cmpK_" + c.name);
            TH1D* sa = SumBins(r.lam, 3, 10, "sL_" + c.name), *sb = SumBins(r.k0s, 3, 10, "sK_" + c.name);
            if (a && b && sa && sb) std::cout << Form("3-10 GeV %-12s Lam sum/run %.3f  K0s sum/run %.3f  ratio sum %.4f run %.4f", c.name.Data(),
                                                     sa->GetBinContent(1) / a->GetBinContent(1), sb->GetBinContent(1) / b->GetBinContent(1),
                                                     0.5 * sa->GetBinContent(1) / sb->GetBinContent(1), 0.5 * a->GetBinContent(1) / b->GetBinContent(1)) << std::endl;
        }
        double plo = 0, phi = 0;
        r.dLam = GetDR(JSFile(dir, c, "Lam", ""), stage, "dLam_" + c.name, plo, phi);
        r.dK0s = GetDR(JSFile(dir, c, "K0s", ""), stage, "dK0s_" + c.name, plo, phi);
        if (r.dLam && r.dK0s) { r.dR = (TH1D*)r.dLam->Clone("dR_" + c.name); r.dR->Divide(r.dK0s); r.dR->Scale(0.5); }
        if (r.lam && r.k0s)   { r.r  = (TH1D*)r.lam->Clone("R_" + c.name);   r.r->Divide(r.k0s);   r.r->Scale(0.5); }
        if (r.lamI && r.k0sI) { r.rI = (TH1D*)r.lamI->Clone("RI_" + c.name); r.rI->Divide(r.k0sI); r.rI->Scale(0.5); }
        res.push_back(r);
    }
    if (!res[0].r) { std::cout << "no MB integrated result: nothing to normalize the double ratio to" << std::endl; }

    // PYTHIA 8: generated jets (V0-free WTA axis) x generated V0s, same processing (jetshape.cpp, run_jetshapes.C)
    // in-cone vs pT with the finest V0 pT binning of the correlation histograms (ptEdges "1.5,2,3,4,5,6,7,8,10"), if processed
    auto mcFile = [&](const char* v0) { TString f = mcDir + "jetshapes_rebinned_ptbins1.5to10n8_" + gPtw + "gen_gen_MC_" + v0 + ".root";
                                        return gSystem->AccessPathName(f) ? mcDir + "jetshapes_rebinned_" + gPtw + "gen_gen_MC_" + v0 + ".root" : f; };
    TH1D *mcL = GetYield(mcFile("Lam"), cone, stage, "mcL"), *mcK = GetYield(mcFile("K0s"), cone, stage, "mcK");
    TH1D *mcR = nullptr, *mcRI = nullptr, *mcLI = SumBins(mcL, 3, 10, "mcLI"), *mcKI = SumBins(mcK, 3, 10, "mcKI");
    if (mcL && mcK) { mcR = (TH1D*)mcL->Clone("R_PYTHIA8"); mcR->Divide(mcK); mcR->Scale(0.5); mcRI = (TH1D*)mcLI->Clone("RI_PYTHIA8"); mcRI->Divide(mcKI); mcRI->Scale(0.5); }
    double mplo = 0, mphi = 0;
    TH1D *mcdL = GetDR(mcDir + "jetshapes_rebinned_" + gPtw + "gen_gen_MC_Lam.root", stage, "mcdL", mplo, mphi), *mcdK = GetDR(mcDir + "jetshapes_rebinned_" + gPtw + "gen_gen_MC_K0s.root", stage, "mcdK", mplo, mphi);
    TH1D* mcdR = nullptr;
    if (mcdL && mcdK) { mcdR = (TH1D*)mcdL->Clone("dR_PYTHIA8"); mcdR->Divide(mcdK); mcdR->Scale(0.5); }
    if (mcR) std::cout << Form("PYTHIA 8 in-cone ratio 3-10 GeV: %.4f +- %.4f", mcRI->GetBinContent(1), mcRI->GetBinError(1)) << std::endl;

    // systematics: in-jet totals (inclusive V0 sources + UE subtraction), fallback the inclusive totals
    std::vector<SysRange> sR = ReadSys(sysFile, "R_jet"), sDR = ReadSys(sysFile, "DR_jet");
    if (sR.empty()) sR = ReadSys(sysFile, "R");
    if (sDR.empty()) sDR = ReadSys(sysFile, "DR");
    std::cout << "systematics (" << sysFile << "): ratio";
    for (auto& s : sR) std::cout << Form("  %g-%g: %.2f%%", s.lo, s.hi, 100 * s.val);
    std::cout << "; double ratio";
    for (auto& s : sDR) std::cout << Form("  %g-%g: %.2f%%", s.lo, s.hi, 100 * s.val);
    std::cout << std::endl;

    const TString coneTxt = cone == "innercone" ? "#Deltar < 0.4" : "#Deltar < 1";
    const TString info = "p_{T}^{jet} > 60 GeV, |#eta^{jet}| < 2, " + coneTxt + ", |#eta^{V0}| < 2";
    const TString coneY = TString(ptw ? "p_{T}-weighted " : "") + (cone == "innercone" ? "in-cone " : "in-jet ");   // y-axis prefix of the ratio and double ratio

    // ratio
    std::vector<TH1D*> hR, hRI;
    for (auto& r : res) { hR.push_back(r.r); hRI.push_back(r.rI); }
    DrawTwoPanel(res, hR, hRI, sR.empty() ? nullptr : &sR, coneY + "(#Lambda+#bar{#Lambda}) / 2K^{0}_{S}", false, "InCone_Ratio", info, false, mcR, mcRI);

    // double ratio to MB integrated (all classes except MB)
    if (res[0].r) {
        std::vector<Res> rs; std::vector<TH1D*> d, dI;
        for (size_t i = 1; i < res.size(); ++i) {
            rs.push_back(res[i]);
            TH1D* a = res[i].r ? (TH1D*)res[i].r->Clone("DR_" + res[i].c.name) : nullptr; if (a) a->Divide(res[0].r);
            TH1D* b = (res[i].rI && res[0].rI) ? (TH1D*)res[i].rI->Clone("DRI_" + res[i].c.name) : nullptr; if (b) b->Divide(res[0].rI);
            d.push_back(a); dI.push_back(b);
        }
        DrawTwoPanel(rs, d, dI, sDR.empty() ? nullptr : &sDR, coneY + "(#Lambda+#bar{#Lambda}) / 2K^{0}_{S}: mult range / MB", true, "InCone_DoubleRatio_MB", info);
    }

    // in-cone yields per jet and GeV (statistical uncertainties only). jetshape.cpp gives the yield per (width of the
    // first pT bin): converted to per GeV, the 3-10 GeV bin divided by its width
    auto perGeV = [](TH1D* h, const TString& n, bool single) {
        if (!h) return (TH1D*)nullptr;
        TH1D* o = (TH1D*)h->Clone(n);
        o->Scale(1.0 / (single ? o->GetBinWidth(o->GetNbinsX()) : o->GetBinWidth(1)));
        return o;
    };
    for (int s = 0; s < 2; ++s) {
        std::vector<TH1D*> h, hI;
        for (auto& r : res) {
            h.push_back(perGeV(s == 0 ? r.lam : r.k0s, Form("y%d_%s", s, r.c.name.Data()), false));
            hI.push_back(perGeV(s == 0 ? r.lamI : r.k0sI, Form("yI%d_%s", s, r.c.name.Data()), true));
        }
        TH1D* m = perGeV(s == 0 ? mcL : mcK, Form("my%d", s), false), *mI = perGeV(s == 0 ? mcLI : mcKI, Form("myI%d", s), true);
        DrawTwoPanel(res, h, hI, nullptr, s == 0 ? "in-cone #Lambda+#bar{#Lambda} per jet per GeV (a.u.)" : "in-cone K^{0}_{S} per jet per GeV (a.u.)", false,
                     s == 0 ? "InCone_Yield_Lam" : "InCone_Yield_K0s", info + ";  stat. only", true, m, mI);
    }

    // Delta r distributions, V0 pT bins summed (1.5-10 GeV), two normalizations, 4 plots each:
    //   P(Delta r)   = dN/dDelta r per jet                                   -> DeltaR_{Lam,K0s,Ratio,DoubleRatio_MB}
    //   rho(Delta r) = dN/dDelta r normalized to unit integral for Delta r < 1 -> JetShape_{Lam,K0s,Ratio,DoubleRatio_MB}
    {
        const TString infoDR = TString("p_{T}^{jet} > 60 GeV, |#eta^{jet}| < 2, 1.5 < p_{T}^{V0} < 10 GeV") + (ptw ? ", p_{T}-weighted" : "");
        const double xhi = 1.0;   // in-jet region (the UE-subtracted distributions are compatible with 0 beyond)
        const double sR15 = SysIn(sR, 1.5, 10), sDR15 = SysIn(sDR, 1.5, 10);   // largest of the pT ranges (conservative)
        auto unit = [&](TH1D* h, const TString& n) {
            if (!h) return (TH1D*)nullptr;
            TH1D* o = (TH1D*)h->Clone(n); o->SetDirectory(0);
            const double I = o->Integral(1, o->FindBin(xhi - 1e-6), "width");
            if (I > 0) o->Scale(1.0 / I);
            return o;
        };
        // classes in drMerge (comma separated; not MB integrated, the double ratio reference) get the ratio and double
        // ratio on the coarser drEdges, with MB integrated merged the same way for their double ratio
        std::vector<double> drE;
        { TObjArray* t = drEdges.Tokenize(","); for (int i = 0; i < t->GetEntries(); ++i) drE.push_back(((TObjString*)t->At(i))->GetString().Atof()); delete t; }
        auto merge = [&](const TString& name) {
            if (drE.size() < 2 || name == res[0].c.name) return false;
            TObjArray* t = drMerge.Tokenize(","); bool m = false;
            for (int i = 0; i < t->GetEntries(); ++i) if (((TObjString*)t->At(i))->GetString() == name) m = true;
            delete t; return m;
        };
        auto ratio = [](TH1D* l, TH1D* k, double f, const TString& n) {
            if (!l || !k) return (TH1D*)nullptr;
            TH1D* o = (TH1D*)l->Clone(n); o->SetDirectory(0); o->Divide(k); o->Scale(f); return o;
        };
        for (int norm = 0; norm < 2; ++norm) {   // 0: P (per jet), 1: rho (unit integral)
            const TString tag = norm == 0 ? "DeltaR_" : "JetShape_", sym = norm == 0 ? (ptw ? "P" : "#Rho") : "#rho";
            const double f = norm == 0 ? 0.5 : 1.0;   // (Lam+AntiLam)/2K0S for the yields; ratio of unit-normalized shapes
            std::vector<TH1D*> l, k, r, r0;   // r0: MB integrated ratio on the binning of r[i] (double ratio denominator)
            for (auto& x : res) {
                TH1D* a = norm == 0 ? x.dLam : unit(x.dLam, "uL_" + x.c.name);
                TH1D* b = norm == 0 ? x.dK0s : unit(x.dK0s, "uK_" + x.c.name);
                l.push_back(a); k.push_back(b);
                if (merge(x.c.name) && a && b && l[0] && k[0]) {   // ratio and double ratio on the coarser Delta r bins
                    r.push_back(ratio(MergeDR(a, drE, "mL_" + x.c.name), MergeDR(b, drE, "mK_" + x.c.name), f, tag + "R_" + x.c.name));
                    r0.push_back(ratio(MergeDR(l[0], drE, "mL0_" + x.c.name), MergeDR(k[0], drE, "mK0_" + x.c.name), f, tag + "R0_" + x.c.name));
                } else {
                    r.push_back(ratio(a, b, f, tag + "R_" + x.c.name)); r0.push_back(nullptr);
                }
            }
            TH1D* mL = norm == 0 ? mcdL : unit(mcdL, "uL_PYTHIA8");
            TH1D* mK = norm == 0 ? mcdK : unit(mcdK, "uK_PYTHIA8");
            TH1D* mR = ratio(mL, mK, f, tag + "R_PYTHIA8");
            const TString yL = sym + "_{#Lambda+#bar{#Lambda}}(#Deltar)", yK = sym + "_{K^{0}_{S}}(#Deltar)";
            const TString yR = norm == 0 ? sym + "_{#Lambda+#bar{#Lambda}} / 2" + sym + "_{K^{0}_{S}}" : sym + "_{#Lambda+#bar{#Lambda}} / " + sym + "_{K^{0}_{S}}";
            DrawDR(res, l, 0, yL, xhi, false, true, tag + "Lam", infoDR + ";  stat. only", mL);
            DrawDR(res, k, 0, yK, xhi, false, true, tag + "K0s", infoDR + ";  stat. only", mK);
            DrawDR(res, r, sR15, yR, xhi, false, false, tag + "Ratio", infoDR, mR);
            if (r[0]) {
                std::vector<Res> rs; std::vector<TH1D*> d;
                for (size_t i = 1; i < res.size(); ++i) {
                    rs.push_back(res[i]);
                    TH1D* a = r[i] ? (TH1D*)r[i]->Clone(tag + "DR_" + res[i].c.name) : nullptr; if (a) a->Divide(r0[i] ? r0[i] : r[0]);
                    d.push_back(a);
                }
                DrawDR(rs, d, sDR15, yR + ": mult range / MB", xhi, true, false, tag + "DoubleRatio_MB", infoDR);
            }
        }
    }

    // values
    std::ofstream txt((pic::gOut + "incone_values.txt").Data());
    txt << "# in-cone (Lambda+anti-Lambda)/2K0S, stage " << stage << ", cone " << cone << ": value +- stat (+- syst)\n";
    for (auto& r : res) {
        txt << Form("%-12s", r.c.name.Data());
        if (r.r) for (int b = 1; b <= r.r->GetNbinsX(); ++b)
            txt << Form("  [%g-%g] %.4f +- %.4f (%.4f)", r.r->GetBinLowEdge(b), r.r->GetXaxis()->GetBinUpEdge(b), r.r->GetBinContent(b), r.r->GetBinError(b),
                        r.r->GetBinContent(b) * SysIn(sR, r.r->GetBinLowEdge(b), r.r->GetXaxis()->GetBinUpEdge(b)));
        if (r.rI) txt << Form("  [3-10] %.4f +- %.4f (%.4f)", r.rI->GetBinContent(r.rI->GetNbinsX()), r.rI->GetBinError(r.rI->GetNbinsX()),
                              r.rI->GetBinContent(r.rI->GetNbinsX()) * SysIn(sR, 3, 10));
        txt << "\n";
    }
    txt.close();
    std::cout << std::ifstream((pic::gOut + "incone_values.txt").Data()).rdbuf();
    TFile* fo = TFile::Open(pic::gOut + "incone_ratios.root", "RECREATE");
    if (fo && !fo->IsZombie()) { for (auto& r : res) for (TH1D* h : {r.lam, r.k0s, r.r, r.lamI, r.k0sI, r.rI, r.dLam, r.dK0s, r.dR}) if (h) h->Write();
                                for (TH1D* h : {mcL, mcK, mcR, mcRI, mcdL, mcdK, mcdR}) if (h) h->Write(); fo->Close(); }
}
