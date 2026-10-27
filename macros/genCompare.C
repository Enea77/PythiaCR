// genCompare.C: generator-level histograms from the official CMS skims or from runPythia trees (same branch names),
// to check that the standalone generator reproduces the official PYTHIA sample.
//   root -l -b -q 'macros/genCompare.C+("<in.root>","<out.root>",<scale>,<pthatLo>,<pthatHi>)'
// Event weight = scale x weight branch, events with pthatLo < pthat <= pthatHi only.
//   official skim: scale = -SigmaRel(pthatGen) (divided by N_events(file) in the macro), window (pthatGen, next
//                  threshold], weight branch = 1 (see macros/runOfficial.sh)
//   runPythia    : scale = 1, no window (pTHat bias weight in the weight branch)
// Jets: gen jets with pT > 60, |eta_WTA| < 2 (forest WTA axis). V0s: |eta| < 2. No UE subtraction.
// All V0 histograms are filled per jet-V0 pair; divide by the jet count (hJetPt integral above 60) to get per-jet yields.
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

void genCompare(const char* in, const char* out, double scale = 1., double pthatLo = -1., double pthatHi = 1e9) {
  TFile* f = TFile::Open(in);
  TTree* t = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
  t->AddFriend("ak4PFJetAnalyzer/t");
  t->AddFriend("HiGenParticleAna/hi");
  const char* v0n[3] = {"K0s", "Lam", "Xi"};
  for (auto n : v0n) t->AddFriend(Form("g%sTree/g%sTree", n, n));
  t->SetBranchStatus("*", 0);

  Float_t pthat, weight; Int_t ngen;
  Float_t genpt[9999], weta[9999], wphi[9999];
  std::vector<double> *gpEta = nullptr, *gpPt = nullptr;
  std::vector<int>* gpChg = nullptr;
  std::vector<double>* vpt[3] = {}; std::vector<double>* veta[3] = {}; std::vector<double>* vphi[3] = {};
  std::vector<double>* lamMom = nullptr;
  auto br = [&](const char* n, void* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
  br("pthat", &pthat); br("weight", &weight);
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
  for (auto h : {hPthat, hPthatJet, hJetPt, hJetEta, hNjetEvt, hNch, hChEta, hChPt}) h->Sumw2();
  for (int s = 0; s < 3; ++s) { hInCone[s]->Sumw2(); hEvt[s]->Sumw2(); hDr[s]->Sumw2(); }

  const Long64_t n = t->GetEntries();
  if (scale < 0) scale = -scale / n;   // official skims: pass -SigmaRel, divided here by the events of the file
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    if (pthat <= pthatLo || pthat > pthatHi) continue;
    const double w = scale * weight;
    hPthat->Fill(pthat, w);
    std::vector<int> jets;
    for (int j = 0; j < ngen; ++j) {
      if (std::fabs(weta[j]) >= kJetEtaMax) continue;
      hJetPt->Fill(genpt[j], w);
      if (genpt[j] > kJetPtMin) { jets.push_back(j); hJetEta->Fill(weta[j], w); }
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
