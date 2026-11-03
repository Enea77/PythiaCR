// genCompare.C: generator-level histograms from the official CMS skims or from runPythia trees (same branch names),
// to check that the standalone generator reproduces the official PYTHIA sample.
//   root -l -b -q 'macros/genCompare.C+("<in.root>","<out.root>",<scale>,<pthatLo>,<pthatHi>)'
// Event weight = scale x weight branch, events with pthatLo < pthat <= pthatHi only.
//   official skim: scale = -SigmaRel(pthatGen) (divided by N_events(file) in the macro), window (pthatGen, next
//                  threshold], weight branch = 1 (see macros/runOfficial.sh)
//   runPythia    : scale = 1, no window (pTHat bias weight in the weight branch)
// Jets: gen jets with pT > 60, |eta_WTA| < 2 (forest WTA axis). V0s: |eta| < 2. No UE subtraction.
// hInCone_*, hEvt_*, hDr_* are filled per jet-V0 pair; divide by the jet count (hJetPt integral above 60) to get per-jet yields.
// Per-event histograms (2026-11-03): hPthatRaw (unweighted), hNchAll / hChEtaAll / hChPtAll (all events; p-Pb samples
// store gen particles only near jets, so these are not inclusive there), hHiBin (hiBin = N_trk^offline in the official
// skims, gen N_ch in runPythia), hV0PtAll_* (all events), hV0PtJetEvt_* (events with a jet > 60), hV0PtNtrk_* (pT vs hiBin,
// all events). Optional 6th argument: maximum number of events (tests).
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include <algorithm>
#include <cstdio>
#include <vector>
#include <string>

namespace {
const double kJetPtMin = 60., kJetEtaMax = 2., kV0EtaMax = 2., kR = 0.4;
double dR(double e1, double p1, double e2, double p2) {
  const double dp = std::remainder(p1 - p2, 2 * TMath::Pi());
  return std::sqrt((e1 - e2) * (e1 - e2) + dp * dp);
}
}  // namespace

void genCompare(const char* in, const char* out, double scale = 1., double pthatLo = -1., double pthatHi = 1e9,
                Long64_t nMax = -1) {
  TFile* f = TFile::Open(in);
  TTree* t = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
  t->AddFriend("ak4PFJetAnalyzer/t");
  t->AddFriend("HiGenParticleAna/hi");
  const char* v0n[3] = {"K0s", "Lam", "Xi"};
  for (auto n : v0n) t->AddFriend(Form("g%sTree/g%sTree", n, n));
  t->SetBranchStatus("*", 0);

  Float_t pthat, weight; Int_t ngen, hiBin = 0;
  Float_t genpt[9999], weta[9999], wphi[9999];
  std::vector<double> *gpEta = nullptr, *gpPt = nullptr;
  std::vector<int>* gpChg = nullptr;
  std::vector<double>* vpt[3] = {}; std::vector<double>* veta[3] = {}; std::vector<double>* vphi[3] = {};
  std::vector<double>* lamMom = nullptr;
  auto br = [&](const char* n, void* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
  br("pthat", &pthat); br("weight", &weight); br("hiBin", &hiBin);
  br("ngen", &ngen); br("genpt", genpt); br("WTAgeneta", weta); br("WTAgenphi", wphi);
  br("eta", &gpEta); br("pt", &gpPt); br("chg", &gpChg);
  for (int s = 0; s < 3; ++s) {
    br(Form("g%s_pt", v0n[s]), &vpt[s]); br(Form("g%s_eta", v0n[s]), &veta[s]); br(Form("g%s_phi", v0n[s]), &vphi[s]);
  }
  br("gLam_mom1", &lamMom);

  TFile fo(out, "RECREATE");
  std::vector<double> jb;
  for (double x = 30; x < 1500; x *= 1.08) jb.push_back(x);
  const double v0b[] = {1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6., 7., 8., 10.};
  // pthat bins: ~5% wide, with the official sample thresholds as exact edges
  const double thr[] = {15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540};
  std::vector<double> pb(std::begin(thr), std::end(thr));
  for (double x = 15 * 1.05; x < 2500; x *= 1.05) {
    bool near = false;
    for (double t : thr) if (std::fabs(x / t - 1) < 0.02) near = true;
    if (!near) pb.push_back(x);
  }
  std::sort(pb.begin(), pb.end());
  TH1D* hPthat = new TH1D("hPthat", ";#hat{p}_{T} (GeV);events", pb.size() - 1, pb.data());
  TH1D* hPthatJet = new TH1D("hPthatJet", ";#hat{p}_{T} (GeV), events with a jet > 60;events", pb.size() - 1, pb.data());
  TH1D* hJetPt = new TH1D("hJetPt", ";gen jet p_{T} (GeV);jets", jb.size() - 1, jb.data());
  TH1D* hJetEta = new TH1D("hJetEta", ";#eta_{WTA} (p_{T} > 60);jets", 40, -2, 2);
  TH1D* hNjetEvt = new TH1D("hNjetEvt", ";events with a jet > 60;", 1, 0, 1);
  TH1D* hNch = new TH1D("hNch", ";N_{ch} (|#eta| < 2.4, p_{T} > 0.4), events with a jet > 60;events", 150, 0, 300);
  TH1D* hChEta = new TH1D("hChEta", ";#eta (charged, p_{T} > 0.4), events with a jet > 60;particles", 48, -2.4, 2.4);
  TH1D* hChPt = new TH1D("hChPt", ";p_{T} (charged, |#eta| < 2.4), events with a jet > 60;particles", 100, 0, 20);
  TH1D *hInCone[3], *hEvt[3];
  TH2D* hDr[3];
  for (int s = 0; s < 3; ++s) {
    hInCone[s] = new TH1D(Form("hInCone_%s", v0n[s]), ";p_{T}^{V0} (GeV);jet-V0 pairs, #Deltar < 0.4", 11, v0b);
    hEvt[s] = new TH1D(Form("hEvt_%s", v0n[s]), ";p_{T}^{V0} (GeV);jet-V0 pairs, any #Deltar", 11, v0b);
    hDr[s] = new TH2D(Form("hDr_%s", v0n[s]), ";#Deltar;p_{T}^{V0} (GeV)", 20, 0, 1, 11, v0b);
  }
  TH1D* hLamSig0 = new TH1D("hLamSig0", ";in-cone #Lambda, 3 < p_{T} < 10: 0 other, 1 from #Sigma^{0};", 2, 0, 2);
  // per-event histograms
  TH1D* hPthatRaw = new TH1D("hPthatRaw", ";#hat{p}_{T} (GeV);events (unweighted)", pb.size() - 1, pb.data());
  TH1D* hNchAll = new TH1D("hNchAll", ";N_{ch} (|#eta| < 2.4, p_{T} > 0.4), all events;events", 150, 0, 300);
  TH1D* hChEtaAll = new TH1D("hChEtaAll", ";#eta (charged, p_{T} > 0.4), all events;particles", 48, -2.4, 2.4);
  TH1D* hChPtAll = new TH1D("hChPtAll", ";p_{T} (charged, |#eta| < 2.4), all events;particles", 100, 0, 20);
  TH1D* hHiBin = new TH1D("hHiBin", ";hiBin (N_{trk});events", 500, 0, 500);
  // inclusive analysis pT edges (v0-jet-analysis/inclusive/plot_inclusive.C) extended to 0.4 and 20 GeV
  const double ib[] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 3.0, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0};
  const int nib = sizeof(ib) / sizeof(ib[0]) - 1;
  TH1D *hV0All[3], *hV0Jet[3];
  TH2D* hV0Ntrk[3];
  for (int s = 0; s < 3; ++s) {
    hV0All[s] = new TH1D(Form("hV0PtAll_%s", v0n[s]), ";p_{T}^{V0} (GeV);V0s (|#eta| < 2), all events", nib, ib);
    hV0Jet[s] = new TH1D(Form("hV0PtJetEvt_%s", v0n[s]), ";p_{T}^{V0} (GeV);V0s (|#eta| < 2), events with a jet > 60", nib, ib);
    hV0Ntrk[s] = new TH2D(Form("hV0PtNtrk_%s", v0n[s]), ";p_{T}^{V0} (GeV);hiBin", nib, ib, 500, 0, 500);
    hV0All[s]->Sumw2(); hV0Jet[s]->Sumw2(); hV0Ntrk[s]->Sumw2();
  }
  for (auto h : {hPthat, hPthatJet, hJetPt, hJetEta, hNjetEvt, hNch, hChEta, hChPt, hPthatRaw, hNchAll, hChEtaAll, hChPtAll, hHiBin})
    h->Sumw2();
  for (int s = 0; s < 3; ++s) { hInCone[s]->Sumw2(); hEvt[s]->Sumw2(); hDr[s]->Sumw2(); }

  const Long64_t n = (nMax > 0 && nMax < t->GetEntries()) ? nMax : t->GetEntries();
  if (scale < 0) scale = -scale / n;   // official skims: pass -SigmaRel, divided here by the events of the file
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    if (pthat <= pthatLo || pthat > pthatHi) continue;
    const double w = scale * weight;
    hPthat->Fill(pthat, w);
    hPthatRaw->Fill(pthat);
    hHiBin->Fill(hiBin, w);
    std::vector<int> jets;
    for (int j = 0; j < ngen; ++j) {
      if (std::fabs(weta[j]) >= kJetEtaMax) continue;
      hJetPt->Fill(genpt[j], w);
      if (genpt[j] > kJetPtMin) { jets.push_back(j); hJetEta->Fill(weta[j], w); }
    }
    int nchAll = 0;
    for (size_t k = 0; k < gpEta->size(); ++k) {
      if ((*gpChg)[k] == 0 || (*gpPt)[k] <= 0.4 || std::fabs((*gpEta)[k]) >= 2.4) continue;
      ++nchAll; hChEtaAll->Fill((*gpEta)[k], w); hChPtAll->Fill((*gpPt)[k], w);
    }
    hNchAll->Fill(nchAll, w);
    for (int s = 0; s < 3; ++s)
      for (size_t k = 0; k < vpt[s]->size(); ++k) {
        if (std::fabs((*veta[s])[k]) >= kV0EtaMax) continue;
        hV0All[s]->Fill((*vpt[s])[k], w);
        hV0Ntrk[s]->Fill((*vpt[s])[k], hiBin, w);
        if (!jets.empty()) hV0Jet[s]->Fill((*vpt[s])[k], w);
      }
    if (jets.empty()) continue;
    hNjetEvt->Fill(0.5, w);
    hPthatJet->Fill(pthat, w);
    int nch = 0;
    for (size_t k = 0; k < gpEta->size(); ++k) {
      if ((*gpChg)[k] == 0 || (*gpPt)[k] <= 0.4 || std::fabs((*gpEta)[k]) >= 2.4) continue;
      ++nch; hChEta->Fill((*gpEta)[k], w); hChPt->Fill((*gpPt)[k], w);
    }
    hNch->Fill(nch, w);
    for (int s = 0; s < 3; ++s)
      for (size_t k = 0; k < vpt[s]->size(); ++k) {
        const double pt = (*vpt[s])[k], eta = (*veta[s])[k], phi = (*vphi[s])[k];
        if (std::fabs(eta) >= kV0EtaMax) continue;
        for (int j : jets) {
          const double d = dR(weta[j], wphi[j], eta, phi);
          hEvt[s]->Fill(pt, w);
          hDr[s]->Fill(d, pt, w);
          if (d < kR) {
            hInCone[s]->Fill(pt, w);
            if (s == 1 && pt > 3 && pt < 10) hLamSig0->Fill(std::abs((*lamMom)[k]) == 3212 ? 1.5 : 0.5, w);
          }
        }
      }
  }
  printf("%s: %lld events, jets > 60: %.4g (weighted), %.0f entries\n", in, n, hJetEta->Integral(), hJetEta->GetEntries());
  fo.Write();
  fo.Close();
}
