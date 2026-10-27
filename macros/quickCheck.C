// quickCheck.C: sanity checks of the runPythia trees (example of a reader; all fills use the event weight).
//   root -l -b -q 'macros/quickCheck.C("output/test_dijet.root")'
// Rough in-cone ratio: V0s within dR < 0.4 of the WTA axis of gen jets with pT > 60, |eta_WTA| < 2,
// 3 < pT(V0) < 10 GeV, |eta(V0)| < 2, no UE subtraction (the analysis uses the correlation method).
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <cstdio>
#include <vector>

void quickCheck(const char* fname) {
  TFile* f = TFile::Open(fname);
  TTree* t = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
  t->AddFriend("ak4PFJetAnalyzer/t");
  t->AddFriend("HiGenParticleAna/hi");
  t->AddFriend("gK0sTree/gK0sTree");
  t->AddFriend("gLamTree/gLamTree");

  Float_t pthat, weight; Int_t hiBin, ngen;
  Float_t genpt[500], weta[500], wphi[500];
  std::vector<double> *eta = nullptr, *k0pt = nullptr, *k0eta = nullptr, *k0phi = nullptr;
  std::vector<double> *lpt = nullptr, *leta = nullptr, *lphi = nullptr, *lmom = nullptr;
  std::vector<int>* pdg = nullptr;
  t->SetBranchAddress("pthat", &pthat); t->SetBranchAddress("weight", &weight); t->SetBranchAddress("hiBin", &hiBin);
  t->SetBranchAddress("ngen", &ngen); t->SetBranchAddress("genpt", genpt);
  t->SetBranchAddress("WTAgeneta", weta); t->SetBranchAddress("WTAgenphi", wphi);
  t->SetBranchAddress("eta", &eta); t->SetBranchAddress("pdg", &pdg);
  t->SetBranchAddress("gK0s_pt", &k0pt); t->SetBranchAddress("gK0s_eta", &k0eta); t->SetBranchAddress("gK0s_phi", &k0phi);
  t->SetBranchAddress("gLam_pt", &lpt); t->SetBranchAddress("gLam_eta", &leta); t->SetBranchAddress("gLam_phi", &lphi);
  t->SetBranchAddress("gLam_mom1", &lmom);

  auto dr = [](double e1, double p1, double e2, double p2) {
    const double dp = std::remainder(p1 - p2, 2 * TMath::Pi());
    return std::sqrt((e1 - e2) * (e1 - e2) + dp * dp);
  };
  double sw = 0, swEta = 0, swJet = 0, swJetLow = 0, nLam = 0, nLamSig0 = 0, nSig0 = 0, inK = 0, inL = 0;
  double wmin = 1e30, wmax = 0, sumMult = 0;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    sw += weight; wmin = std::min<double>(wmin, weight); wmax = std::max<double>(wmax, weight);
    sumMult += weight * hiBin;
    for (double e : *eta) swEta += weight * e / eta->size();
    for (int p : *pdg) if (std::abs(p) == 3212) ++nSig0;
    for (size_t k = 0; k < lpt->size(); ++k) { ++nLam; if (std::abs((*lmom)[k]) == 3212) ++nLamSig0; }
    for (int j = 0; j < ngen; ++j) {
      if (genpt[j] < 60 || std::fabs(weta[j]) > 2) continue;
      swJet += weight;
      if (pthat < 40) swJetLow += weight;
      for (size_t k = 0; k < k0pt->size(); ++k)
        if ((*k0pt)[k] > 3 && (*k0pt)[k] < 10 && std::fabs((*k0eta)[k]) < 2 && dr(weta[j], wphi[j], (*k0eta)[k], (*k0phi)[k]) < 0.4) inK += weight;
      for (size_t k = 0; k < lpt->size(); ++k)
        if ((*lpt)[k] > 3 && (*lpt)[k] < 10 && std::fabs((*leta)[k]) < 2 && dr(weta[j], wphi[j], (*leta)[k], (*lphi)[k]) < 0.4) inL += weight;
    }
  }
  printf("events %lld, weights %.3g - %.3g\n", t->GetEntries(), wmin, wmax);
  printf("<eta> of gen particles (weighted): %.3f (boost toward -z gives < 0)\n", swEta / sw);
  printf("<N_trk gen>: %.1f\n", sumMult / sw);
  printf("Sigma0 in the record: %.0f (should be 0); Lambda from Sigma0: %.1f%%\n", nSig0, 100. * nLamSig0 / nLam);
  printf("jets pT > 60, |eta_WTA| < 2: %.4g weighted; fraction from pthat < 40: %.2f%%\n", swJet, 100. * swJetLow / swJet);
  printf("in-cone (no UE subtraction) Lambda/2K0S, 3-10 GeV: %.3f\n", inL / (2. * inK));
}
