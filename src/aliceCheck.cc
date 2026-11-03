// aliceCheck.cc: reproduce the PYTHIA side of ALICE JHEP07(2023)136 (pp 13 TeV, CR-BLC modes), Figs. 5-7:
// (Lambda + anti-Lambda)/2K0S inclusive and in charged-particle jets (JE = jet cone - perpendicular cone).
//   bin/aliceCheck out.root nEvents seed card1.cmnd [card2.cmnd ...]
// Definitions (paper sections 3.1-3.3):
//   V0s: K0S, Lambda + anti-Lambda with |eta| < 0.75; Xi/Omega stable via common.cmnd (tau0Max = 10 mm), so the Lambda
//        are primary (Sigma0 included), as ALICE's feed-down-corrected Lambda.
//   jets: anti-kT R = 0.4 on charged final-state particles with pT > 0.15 GeV, |eta| < 0.9; pT,jet^ch > 10 GeV,
//         |eta_jet| < 0.35; no UE subtraction of the jet pT (pp).
//   JC: V0 within R < 0.4 of a jet axis; PC: cones of R = 0.4 at the jet eta and phi +- pi/2 (filled with weight 1/2,
//       i.e. the average of the two cones per jet). JE = JC - PC. Event counts: all events / events with >= 1 jet.
// All histograms are filled with the event weight (needed for pTHat-biased hard-QCD samples).
#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "TFile.h"
#include "TH1D.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace Pythia8;

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cerr << "usage: aliceCheck out.root nEvents seed card1.cmnd [card2.cmnd ...]" << std::endl;
    return 1;
  }
  const std::string out = argv[1];
  const long nEv = std::atol(argv[2]);
  const int seed = std::atoi(argv[3]);
  Pythia pythia;
  for (int i = 4; i < argc; ++i)
    if (!pythia.readFile(argv[i])) { std::cerr << "cannot read " << argv[i] << std::endl; return 1; }
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  pythia.readString("Next:numberCount = 20000");
  if (!pythia.init()) return 1;

  const double kV0Eta = 0.75, kTrkEta = 0.9, kTrkPt = 0.15, kR = 0.4, kJetPt = 10., kJetEta = 0.35;
  const double pb[] = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0};
  const int npb = sizeof(pb) / sizeof(pb[0]) - 1;
  TFile fo(out.c_str(), "RECREATE");
  TH1D* hEv = new TH1D("hEv", "1: all events, 2: events with >= 1 jet (sum of weights)", 2, 0, 2);
  TH1D* hJetPt = new TH1D("hJetPt", ";p_{T,jet}^{ch} (GeV);jets", 100, 0, 100);
  const char* sp[2] = {"K0s", "Lam"};
  TH1D *hIncl[2], *hJC[2], *hPC[2];
  for (int s = 0; s < 2; ++s) {
    hIncl[s] = new TH1D(Form("hIncl_%s", sp[s]), ";p_{T} (GeV);V0s, |#eta| < 0.75, all events", npb, pb);
    hJC[s] = new TH1D(Form("hJC_%s", sp[s]), ";p_{T} (GeV);V0s in jet cones (R < 0.4)", npb, pb);
    hPC[s] = new TH1D(Form("hPC_%s", sp[s]), ";p_{T} (GeV);V0s in perpendicular cones (per jet, avg of 2)", npb, pb);
    hIncl[s]->Sumw2(); hJC[s]->Sumw2(); hPC[s]->Sumw2();
  }
  hEv->Sumw2(); hJetPt->Sumw2();
  const fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, kR);

  struct V0 { int s; double pt, eta, phi; };
  for (long iEv = 0; iEv < nEv; ++iEv) {
    if (!pythia.next()) continue;
    const double w = pythia.info.weight();
    const Event& ev = pythia.event;
    std::vector<fastjet::PseudoJet> trk;
    std::vector<V0> v0s;
    for (int i = 0; i < ev.size(); ++i) {
      const Particle& p = ev[i];
      if (!p.isFinal()) continue;
      const int a = std::abs(p.id());
      if (a == 310 || a == 3122) {
        if (std::abs(p.eta()) < kV0Eta) v0s.push_back({a == 310 ? 0 : 1, p.pT(), p.eta(), p.phi()});
        continue;   // neutral: not a jet constituent
      }
      if (p.isCharged() && p.pT() > kTrkPt && std::abs(p.eta()) < kTrkEta) trk.emplace_back(p.px(), p.py(), p.pz(), p.e());
    }
    hEv->Fill(0.5, w);
    for (const auto& v : v0s) hIncl[v.s]->Fill(v.pt, w);
    fastjet::ClusterSequence cs(trk, jetDef);
    std::vector<fastjet::PseudoJet> jets;
    for (const auto& j : fastjet::sorted_by_pt(cs.inclusive_jets(kJetPt)))
      if (std::abs(j.eta()) < kJetEta) jets.push_back(j);
    if (jets.empty()) continue;
    hEv->Fill(1.5, w);
    for (const auto& j : jets) {
      hJetPt->Fill(j.pt(), w);
      for (const auto& v : v0s) {
        const double de = v.eta - j.eta();
        if (de * de + std::pow(std::remainder(v.phi - j.phi_std(), 2 * M_PI), 2) < kR * kR) hJC[v.s]->Fill(v.pt, w);
        for (double sgn : {-1., 1.}) {
          const double dp = std::remainder(v.phi - (j.phi_std() + sgn * M_PI / 2), 2 * M_PI);
          if (de * de + dp * dp < kR * kR) hPC[v.s]->Fill(v.pt, 0.5 * w);
        }
      }
    }
  }
  pythia.stat();
  fo.Write();
  fo.Close();
  std::cout << "aliceCheck: wrote " << out << std::endl;
  return 0;
}
