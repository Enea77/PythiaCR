// plotCompare.C: overlay genCompare.C outputs (first file = reference, ratio panels vs it).
//   root -l -b -q 'macros/plotCompare.C("compare/official_pgoing.root,compare/dijet_cuetp8m1.root","CMS official,standalone CUETP8M1","compare/plots")'
// Normalisation: jet spectra to unit area above 60 GeV; event quantities per event with a jet > 60;
// V0 yields per jet (pT > 60, |eta_WTA| < 2). Prints the integrated in-cone Lambda/2K0S (3-10 GeV) and other numbers.
// Binning: neighbouring bins are merged (same edges for all samples) until every merged bin has a relative statistical
// error <= maxRel in every sample (2% by default, 5% for Xi); the raw weighted sums are merged before any normalisation,
// and Lambda/2K0S ratios are formed after merging Lambda and K0S separately.
// Ratio panels: y range from the lowest / highest ratio point (including its error), with 15% padding.
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TObjString.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

namespace {
const int kCol[8] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet - 3};

// variable bin edges for [xlo, xhi): merge original bins until each merged bin has relative error <= maxRel in all
// histograms (a tail that never reaches it is merged into the last good bin)
std::vector<double> VarEdges(const std::vector<TH1*>& hs, double maxRel, double xlo = -1e30, double xhi = 1e30) {
  const TAxis* ax = hs[0]->GetXaxis();
  const int b1 = std::max(1, ax->FindBin(xlo + 1e-9)), b2 = std::min(ax->GetNbins(), ax->FindBin(xhi - 1e-9));
  std::vector<double> e = {ax->GetBinLowEdge(b1)};
  std::vector<double> s(hs.size(), 0), e2(hs.size(), 0);
  for (int b = b1; b <= b2; ++b) {
    bool ok = true;
    for (size_t i = 0; i < hs.size(); ++i) {
      s[i] += hs[i]->GetBinContent(b); e2[i] += std::pow(hs[i]->GetBinError(b), 2);
      if (s[i] <= 0 || std::sqrt(e2[i]) / s[i] > maxRel) ok = false;
    }
    if (ok) { e.push_back(ax->GetBinUpEdge(b)); std::fill(s.begin(), s.end(), 0); std::fill(e2.begin(), e2.end(), 0); }
  }
  if (e.size() == 1) e.push_back(ax->GetBinUpEdge(b2));                // never reached: one bin
  else if (e.back() < ax->GetBinUpEdge(b2) - 1e-9) e.back() = ax->GetBinUpEdge(b2);   // tail into the last bin
  return e;
}
TH1D* Rebinned(TH1* h, const std::vector<double>& e, const char* name) {
  TH1D* r = (TH1D*)h->Rebin(e.size() - 1, name, e.data());
  r->SetDirectory(0);
  return r;
}

void overlay(std::vector<TH1D*> h, const std::vector<std::string>& lab, const std::string& out, bool logy, bool logx = false) {
  TCanvas c("c", "", 700, 700);
  TPad top("t", "", 0, 0.3, 1, 1), bot("b", "", 0, 0, 1, 0.3);
  top.SetBottomMargin(0.02); top.SetLeftMargin(0.13); bot.SetTopMargin(0.02); bot.SetBottomMargin(0.3); bot.SetLeftMargin(0.13);
  if (logx) { top.SetLogx(); bot.SetLogx(); }
  top.Draw(); bot.Draw();
  top.cd(); if (logy) top.SetLogy();
  TLegend leg(0.55, 0.72, 0.88, 0.88); leg.SetBorderSize(0);
  double mx = 0, mn = 1e300;
  for (auto x : h) for (int b = 1; b <= x->GetNbinsX(); ++b) if (x->GetBinContent(b) > 0) { mx = std::max(mx, x->GetBinContent(b)); mn = std::min(mn, x->GetBinContent(b)); }
  for (size_t i = 0; i < h.size(); ++i) {
    h[i]->SetLineColor(kCol[i % 8]); h[i]->SetMarkerColor(kCol[i % 8]); h[i]->SetMarkerStyle(20 + i); h[i]->SetMarkerSize(0.8);
    h[i]->SetMaximum(logy ? mx * 5 : mx * 1.3);
    h[i]->SetMinimum(logy ? mn * 0.5 : 0);
    h[i]->GetXaxis()->SetLabelSize(0);
    h[i]->Draw(i ? "e same" : "e");
    leg.AddEntry(h[i], lab[i].c_str(), "lp");
  }
  leg.Draw();
  bot.cd();
  std::vector<TH1D*> rs;
  double lo = 1e30, hi = -1e30;
  for (size_t i = 1; i < h.size(); ++i) {
    TH1D* r = (TH1D*)h[i]->Clone(Form("%s_r%zu", h[i]->GetName(), i));
    r->Divide(h[0]);
    for (int b = r->GetXaxis()->GetFirst(); b <= r->GetXaxis()->GetLast(); ++b) {
      if (r->GetBinContent(b) <= 0 || h[0]->GetBinContent(b) <= 0) continue;
      lo = std::min(lo, r->GetBinContent(b) - r->GetBinError(b));
      hi = std::max(hi, r->GetBinContent(b) + r->GetBinError(b));
    }
    rs.push_back(r);
  }
  if (lo > hi) { lo = 0.9; hi = 1.1; }
  lo = std::min(lo, 1.0); hi = std::max(hi, 1.0);
  const double pad = 0.15 * (hi - lo);
  for (size_t i = 0; i < rs.size(); ++i) {
    TH1D* r = rs[i];
    r->SetMinimum(lo - pad); r->SetMaximum(hi + pad);
    r->GetYaxis()->SetTitle(Form("ratio to %s", lab[0].c_str()));
    r->GetYaxis()->SetNdivisions(505);
    r->GetYaxis()->SetTitleSize(0.08); r->GetYaxis()->SetLabelSize(0.08); r->GetYaxis()->SetTitleOffset(0.6);
    r->GetXaxis()->SetLabelSize(0.1); r->GetXaxis()->SetTitleSize(0.1);
    r->Draw(i ? "e same" : "e");
  }
  if (!rs.empty()) {
    TLine* one = new TLine(rs[0]->GetXaxis()->GetBinLowEdge(rs[0]->GetXaxis()->GetFirst()), 1,
                           rs[0]->GetXaxis()->GetBinUpEdge(rs[0]->GetXaxis()->GetLast()), 1);
    one->SetLineStyle(2); one->Draw();
  }
  c.SaveAs((out + ".png").c_str());
  c.SaveAs((out + ".pdf").c_str());
}
}  // namespace

void plotCompare(const char* files, const char* labels, const char* outDir, double maxRel = 0.02) {
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);
  gSystem->mkdir(outDir, true);
  std::vector<TFile*> F; std::vector<std::string> L;
  for (auto o : *TString(files).Tokenize(",")) F.push_back(TFile::Open(((TObjString*)o)->GetString()));
  for (auto o : *TString(labels).Tokenize(",")) L.push_back(((TObjString*)o)->GetString().Data());
  const std::string d = outDir;
  const size_t n = F.size();

  auto raw = [&](const char* name) { std::vector<TH1*> v; for (auto f : F) v.push_back((TH1*)f->Get(name)); return v; };
  // merged, normalised (1/norm per unit x) versions of one histogram for all samples
  auto spectra = [&](const char* name, const std::vector<double>& norm, double rel, double xlo = -1e30, double xhi = 1e30) {
    const std::vector<TH1*> r = raw(name);
    const std::vector<double> e = VarEdges(r, rel, xlo, xhi);
    std::vector<TH1D*> o;
    for (size_t i = 0; i < n; ++i) { TH1D* h = Rebinned(r[i], e, Form("%s_%zu", name, i)); h->Scale(1. / norm[i], "width"); o.push_back(h); }
    return o;
  };
  // Lambda / 2K0S with common merged edges (decided on Lambda and K0S of all samples)
  auto lamK0s = [&](const std::vector<TH1*>& l, const std::vector<TH1*>& k, const char* tag, const char* ytitle) {
    std::vector<TH1*> both(l); both.insert(both.end(), k.begin(), k.end());
    const std::vector<double> e = VarEdges(both, maxRel);
    std::vector<TH1D*> o;
    for (size_t i = 0; i < n; ++i) {
      TH1D* r = Rebinned(l[i], e, Form("%s_%zu", tag, i));
      r->Divide(Rebinned(k[i], e, Form("%s_k_%zu", tag, i))); r->Scale(0.5);
      r->GetYaxis()->SetTitle(ytitle);
      o.push_back(r);
    }
    return o;
  };

  std::vector<double> nJet(n), nEvt(n), nJetPt(n), one(n, 1.);
  printf("%-28s %10s %8s %10s %10s %10s %10s\n", "sample", "jets>60", "<Nch>", "L/2K 3-10", "Sig0 frac", "K0S/jet", "Lam/jet");
  for (size_t i = 0; i < n; ++i) {
    TH1D* je = (TH1D*)F[i]->Get("hJetEta");
    nJet[i] = je->Integral(); nEvt[i] = ((TH1D*)F[i]->Get("hNjetEvt"))->Integral();
    TH1D* jp = (TH1D*)F[i]->Get("hJetPt");
    nJetPt[i] = jp->Integral(jp->FindBin(60.01), jp->GetNbinsX());
    TH1D* l = (TH1D*)F[i]->Get("hInCone_Lam"); TH1D* k = (TH1D*)F[i]->Get("hInCone_K0s");
    const int b1 = l->FindBin(3.01), b2 = l->FindBin(9.99);
    double el, ek;
    const double il = l->IntegralAndError(b1, b2, el), ik = k->IntegralAndError(b1, b2, ek);
    TH1D* s0 = (TH1D*)F[i]->Get("hLamSig0");
    printf("%-28s %10.0f %8.1f %6.4f+-%.4f %10.3f %10.4f %10.4f\n", L[i].c_str(), je->GetEntries(), ((TH1D*)F[i]->Get("hNch"))->GetMean(),
           il / (2 * ik), il / (2 * ik) * std::sqrt(std::pow(el / il, 2) + std::pow(ek / ik, 2)),
           s0->GetBinContent(2) / s0->Integral(), ik / nJet[i], il / nJet[i]);
  }

  // pthat: shapes (unit area) and sigma(pthat > threshold) / sigma(pthat > 15), compared with CorrCfg::SigmaRel
  const int nThr = 10;
  const double thr[nThr] = {30, 50, 80, 120, 170, 220, 280, 370, 460, 540};
  const double sigRel[nThr] = {7.4232802e-2, 9.5408304e-3, 1.2977183e-3, 2.1848571e-4, 4.0848861e-5,
                               1.1474575e-5, 3.3198812e-6, 7.3406939e-7, 2.0944659e-7, 7.9191586e-8};
  std::vector<double> nPt(n), nPtJ(n);
  printf("\nsigma(pthat > X) / sigma(pthat > 15):\n%-28s", "pthat >");
  for (double x : thr) printf(" %10.0f", x);
  printf("\n%-28s", "CorrCfg::SigmaRel (McM)");
  for (double s : sigRel) printf(" %10.3e", s);
  printf("\n");
  for (size_t i = 0; i < n; ++i) {
    TH1D* h = (TH1D*)F[i]->Get("hPthat");
    nPt[i] = h->Integral(); nPtJ[i] = ((TH1D*)F[i]->Get("hPthatJet"))->Integral();
    printf("%-28s", L[i].c_str());
    for (double x : thr) printf(" %10.3e", h->Integral(h->FindBin(x * 1.0001), h->GetNbinsX() + 1) / nPt[i]);   // thresholds are bin edges
    printf("\n");
  }
  overlay(spectra("hPthat", nPt, maxRel), L, d + "/pthat", true, true);
  overlay(spectra("hPthatJet", nPtJ, maxRel), L, d + "/pthat_jetEvents", true, true);
  overlay(spectra("hJetPt", nJetPt, maxRel, 60, 1500), L, d + "/jetPt", true, true);
  overlay(spectra("hJetEta", nJet, maxRel), L, d + "/jetEta", false);
  overlay(spectra("hNch", nEvt, maxRel), L, d + "/Nch", false);
  overlay(spectra("hChEta", nEvt, maxRel), L, d + "/chEta", false);
  overlay(spectra("hChPt", nEvt, maxRel), L, d + "/chPt", true);
  const char* v0n[3] = {"K0s", "Lam", "Xi"};
  for (int s = 0; s < 3; ++s) overlay(spectra(Form("hInCone_%s", v0n[s]), nJet, s == 2 ? 0.05 : maxRel), L, d + "/inCone_" + v0n[s], true);
  overlay(lamK0s(raw("hInCone_Lam"), raw("hInCone_K0s"), "lk", "#Lambda/2K^{0}_{S}, #Deltar < 0.4"), L, d + "/LamK0s_inCone", false);
  overlay(lamK0s(raw("hEvt_Lam"), raw("hEvt_K0s"), "lkEvt", "#Lambda/2K^{0}_{S}, jet events"), L, d + "/LamK0s_jetEvents", false);
  // Delta r profile, 1.5 < pT < 10 GeV
  std::vector<TH1*> pl, pk;
  for (size_t i = 0; i < n; ++i) {
    TH2D* l = (TH2D*)F[i]->Get("hDr_Lam"); TH2D* k = (TH2D*)F[i]->Get("hDr_K0s");
    const int b1 = l->GetYaxis()->FindBin(1.51), b2 = l->GetYaxis()->FindBin(9.99);
    pl.push_back(l->ProjectionX(Form("pl_%zu", i), b1, b2, "e")); pk.push_back(k->ProjectionX(Form("pk_%zu", i), b1, b2, "e"));
  }
  overlay(lamK0s(pl, pk, "lkDr", "#Lambda/2K^{0}_{S}, 1.5 < p_{T} < 10 GeV"), L, d + "/LamK0s_dr", false);
}
