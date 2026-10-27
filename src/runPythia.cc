// Generator-level PYTHIA 8 trees in the layout of the CMS pPb 2016 V0+jet skims (gen branches only),
// so that the histograms are made afterwards by a separate reader (any binning, any selection).
//
// Usage: runPythia <out.root> <nEvents> <seed> <card1.cmnd> [card2 ...]
// Cards are read in order, so later cards override earlier ones.
//
// Trees (same directory / tree / branch names and types as the skims):
//   hiEvtAnalyzer/HiTree   pthat, weight (PYTHIA event weight, != 1 with the pTHat bias), vz (= 0),
//                          hiBin (gen N_trk: charged final-state particles, |eta| < 2.4, pT > 0.4 GeV), nMPI
//   ak4PFJetAnalyzer/t     gen jets: anti-kT R = 0.4 on all final-state particles except neutrinos (no eta cut),
//                          pT > kJetPtStore; E-scheme axis (geneta, genphi) and WTA axis (WTAgeneta, WTAgenphi)
//   HiGenParticleAna/hi    final-state particles with |eta| < 2.4, pT > 0.4 GeV (as in the skims); with
//                          PythiaCR:genParticleJetR > 0 only those within that Delta R of a stored jet (p-Pb)
//   g{K0s,Lam,Xi,Om}Tree   K0S, Lambda, Xi-, Omega- (and antiparticles) with |eta| < 2.4; mom1 = mother PDG id
//                          if the mother is a hadron (e.g. Sigma0 -> Lambda gamma), 0 if from a string
// Histogram "hWeights": bin 1 events, bin 2 sum of weights, bin 3 sigmaGen (mb) x events (average it after hadd).
// p-Pb (Angantyr): weight = impact-parameter weight x pTHat bias weight (see below).

#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace Pythia8;

namespace {
const double kEtaGen = 2.4;      // gen particles and gen V0s stored
const double kPtGen = 0.4;       // gen particles stored and N_trk
const double kJetR = 0.4;
const double kJetPtStore = 30.;  // gen jets stored (analysis cut 60 GeV is applied by the reader)
const int kMaxJets = 500;

bool isNeutrino(int id) {
  const int a = std::abs(id);
  return a == 12 || a == 14 || a == 16;
}

// PDG id of the mother if it is a hadron, 0 otherwise (string / cluster)
int hadronMother(const Event& ev, const Particle& p) {
  const int m = p.mother1();
  if (m <= 0) return 0;
  const int id = ev[m].id();
  return ev[m].isHadron() ? id : 0;
}

struct V0Branches {
  std::vector<double> pt, eta, phi, mass, id, mom1;
  void clear() { pt.clear(); eta.clear(); phi.clear(); mass.clear(); id.clear(); mom1.clear(); }
  void push(const Event& ev, const Particle& p) {
    pt.push_back(p.pT()); eta.push_back(p.eta()); phi.push_back(p.phi());
    mass.push_back(p.m()); id.push_back(p.id()); mom1.push_back(hadronMother(ev, p));
  }
  TTree* book(TFile& f, const std::string& n, bool withId) {
    f.mkdir((n + "Tree").c_str())->cd();
    TTree* t = new TTree((n + "Tree").c_str(), "");
    t->Branch((n + "_pt").c_str(), &pt);
    t->Branch((n + "_eta").c_str(), &eta);
    t->Branch((n + "_phi").c_str(), &phi);
    t->Branch((n + "_mass").c_str(), &mass);
    if (withId) t->Branch((n + "_id").c_str(), &id);
    t->Branch((n + "_mom1").c_str(), &mom1);
    return t;
  }
};
}  // namespace

int main(int argc, char* argv[]) {
  if (argc < 5) {
    std::cerr << "Usage: " << argv[0] << " <out.root> <nEvents> <seed> <card1.cmnd> [card2 ...]\n";
    return 1;
  }
  const std::string outName = argv[1];
  const int nEvents = std::atoi(argv[2]);
  const int seed = std::atoi(argv[3]);

  Pythia pythia;
  // PythiaCR:genParticleJetR > 0: store gen particles only within this Delta R of a stored gen jet (WTA axis); used
  // for p-Pb to keep the files small (the correlation code needs them only near jets; N_trk is in hiBin). 0 = all.
  pythia.settings.addParm("PythiaCR:genParticleJetR", 0., true, false, 0., 10.);
  for (int i = 4; i < argc; ++i)
    if (!pythia.readFile(argv[i])) return 1;
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  // With ion beams (Angantyr) info.weight() carries only the impact-parameter weight: the pTHat bias of the signal
  // sub-collision (PhaseSpace:bias2Selection) is not propagated, so its weight (pTHatRef / pTHat)^pow is applied here
  // (the same weight PYTHIA gives in pp, where it is already contained in info.weight()).
  // Read before init(): Angantyr resets Beams:idA/idB of the main object to nucleons during initialisation.
  const bool isIon = std::abs(pythia.settings.mode("Beams:idA")) > 1000000000 || std::abs(pythia.settings.mode("Beams:idB")) > 1000000000;
  const bool addBias = isIon && pythia.settings.flag("PhaseSpace:bias2Selection");
  const double biasRef = pythia.settings.parm("PhaseSpace:bias2SelectionRef"), biasPow = pythia.settings.parm("PhaseSpace:bias2SelectionPow");
  if (!pythia.init()) return 1;
  if (addBias) std::cout << "runPythia: ion beams, pTHat bias weight (" << biasRef << "/pTHat)^" << biasPow << " applied to the event weight" << std::endl;
  const double gpJetR = pythia.settings.parm("PythiaCR:genParticleJetR");
  if (gpJetR > 0) std::cout << "runPythia: gen particles stored only within Delta R < " << gpJetR << " of a gen jet" << std::endl;

  TFile fout(outName.c_str(), "RECREATE");

  // event
  Float_t pthat, weight, vz = 0;
  Int_t hiBin, nMPI;
  ULong64_t evt;
  fout.mkdir("hiEvtAnalyzer")->cd();
  TTree* tEvt = new TTree("HiTree", "");
  tEvt->Branch("evt", &evt, "evt/l");
  tEvt->Branch("vz", &vz, "vz/F");
  tEvt->Branch("hiBin", &hiBin, "hiBin/I");
  tEvt->Branch("pthat", &pthat, "pthat/F");
  tEvt->Branch("weight", &weight, "weight/F");
  tEvt->Branch("nMPI", &nMPI, "nMPI/I");

  // gen jets
  Int_t ngen;
  Float_t genpt[kMaxJets], geneta[kMaxJets], genphi[kMaxJets], wtaeta[kMaxJets], wtaphi[kMaxJets];
  fout.mkdir("ak4PFJetAnalyzer")->cd();
  TTree* tJet = new TTree("t", "");
  tJet->Branch("ngen", &ngen, "ngen/I");
  tJet->Branch("genpt", genpt, "genpt[ngen]/F");
  tJet->Branch("geneta", geneta, "geneta[ngen]/F");
  tJet->Branch("genphi", genphi, "genphi[ngen]/F");
  tJet->Branch("WTAgeneta", wtaeta, "WTAgeneta[ngen]/F");
  tJet->Branch("WTAgenphi", wtaphi, "WTAgenphi[ngen]/F");

  // gen particles
  std::vector<double> gpPt, gpEta, gpPhi;
  std::vector<int> gpPdg, gpChg;
  fout.mkdir("HiGenParticleAna")->cd();
  TTree* tGen = new TTree("hi", "");
  tGen->Branch("pt", &gpPt);
  tGen->Branch("eta", &gpEta);
  tGen->Branch("phi", &gpPhi);
  tGen->Branch("pdg", &gpPdg);
  tGen->Branch("chg", &gpChg);

  V0Branches k0s, lam, xi, om;
  TTree* tK0s = k0s.book(fout, "gK0s", false);
  TTree* tLam = lam.book(fout, "gLam", true);
  TTree* tXi = xi.book(fout, "gXi", true);
  TTree* tOm = om.book(fout, "gOm", true);

  fout.cd();
  TH1D* hWeights = new TH1D("hWeights", "1: events, 2: sum of weights, 3: sigmaGen (mb) x events", 3, 0., 3.);

  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, kJetR);
  fastjet::JetDefinition wtaDef(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R,
                                fastjet::WTA_pt_scheme);

  for (int iEv = 0; iEv < nEvents; ++iEv) {
    if (!pythia.next()) continue;
    const Event& ev = pythia.event;
    evt = iEv;
    pthat = pythia.info.pTHat();
    weight = pythia.info.weight();
    if (addBias) weight *= std::pow(biasRef / pthat, biasPow);
    nMPI = pythia.info.nMPI();

    std::vector<fastjet::PseudoJet> inputs;
    gpPt.clear(); gpEta.clear(); gpPhi.clear(); gpPdg.clear(); gpChg.clear();
    k0s.clear(); lam.clear(); xi.clear(); om.clear();
    hiBin = 0;
    for (int i = 0; i < ev.size(); ++i) {
      const Particle& p = ev[i];
      if (!p.isFinal()) continue;
      if (!isNeutrino(p.id())) inputs.emplace_back(p.px(), p.py(), p.pz(), p.e());
      const double aeta = std::abs(p.eta());
      if (aeta >= kEtaGen) continue;
      if (p.pT() > kPtGen) {
        gpPt.push_back(p.pT()); gpEta.push_back(p.eta()); gpPhi.push_back(p.phi());
        gpPdg.push_back(p.id()); gpChg.push_back(p.chargeType() / 3);
        if (p.isCharged()) ++hiBin;
      }
      switch (std::abs(p.id())) {
        case 310: k0s.push(ev, p); break;
        case 3122: lam.push(ev, p); break;
        case 3312: xi.push(ev, p); break;
        case 3334: om.push(ev, p); break;
      }
    }

    fastjet::ClusterSequence cs(inputs, jetDef);
    ngen = 0;
    for (const auto& j : fastjet::sorted_by_pt(cs.inclusive_jets(kJetPtStore))) {
      if (ngen >= kMaxJets) break;
      fastjet::ClusterSequence csWta(j.constituents(), wtaDef);
      const fastjet::PseudoJet wta = csWta.exclusive_jets(1).at(0);
      genpt[ngen] = j.pt(); geneta[ngen] = j.eta(); genphi[ngen] = j.phi_std();
      wtaeta[ngen] = wta.eta(); wtaphi[ngen] = wta.phi_std();
      ++ngen;
    }
    if (gpJetR > 0) {   // keep only the gen particles near a stored jet
      size_t k = 0;
      for (size_t i = 0; i < gpPt.size(); ++i) {
        bool near = false;
        for (int j = 0; j < ngen && !near; ++j) {
          const double dphi = std::remainder(gpPhi[i] - wtaphi[j], 2. * M_PI), deta = gpEta[i] - wtaeta[j];
          near = deta * deta + dphi * dphi < gpJetR * gpJetR;
        }
        if (!near) continue;
        gpPt[k] = gpPt[i]; gpEta[k] = gpEta[i]; gpPhi[k] = gpPhi[i]; gpPdg[k] = gpPdg[i]; gpChg[k] = gpChg[i]; ++k;
      }
      gpPt.resize(k); gpEta.resize(k); gpPhi.resize(k); gpPdg.resize(k); gpChg.resize(k);
    }

    tEvt->Fill(); tJet->Fill(); tGen->Fill();
    tK0s->Fill(); tLam->Fill(); tXi->Fill(); tOm->Fill();
    hWeights->Fill(0.5);
    hWeights->Fill(1.5, weight);
  }

  hWeights->SetBinContent(3, pythia.info.sigmaGen() * hWeights->GetBinContent(1));
  pythia.stat();
  fout.Write();
  fout.Close();
  return 0;
}
