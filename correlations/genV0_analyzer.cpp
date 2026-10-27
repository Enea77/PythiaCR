// genV0_analyzer.cpp -- gen jets x gen V0 correlations from the runPythia trees, with the gen path of
// v0-jet-analysis/correlations/jetV0_analyzer.cpp (commit 10af763, see ../COPIED_FROM_v0-jet-analysis.txt):
// same V0 and jet selection, V0-free WTA reclustering (Corr::ReclusterWTA), same correlation and jet sparses
// (MakeCorr / MakeJet below are copied from the analyzer), same filling (Corr::FillSame) and event mixing
// (Corr::MixPool, CorrCfg::NMix("MC"), kMixDVz, kMixDMult, MixDepth). Output layout as the analyzer for MC, gen part only:
//   QA_histograms/hist_gen_jet_weighted, Nevents, pthathist_weighted, multiplicity_weighted
//   correlation_gen_gen_histograms/hist_correlation_{signal,mixing}_jet_gen_V0_gen[_ptw], hist_jet_from_gen_gen_{sig,mix}
// so that postprocessing/jetshape.cpp("gen","gen",...) runs on it unchanged.
//
// Differences from the analyzer (generator-level input, no detector):
//   - event weight = PYTHIA event weight (pTHat bias; replaces cross-section x vz weight), no event filters, no vz cut
//   - vz = 0 in the trees: the mixing pool is then classified by multiplicity only
//   - multiplicity = ntrkScale x hiBin (generator-level N_trk: charged, |eta| < 2.4, pT > 0.4) instead of the reco N_trk^offline
//   - jet flavour = 2 (other/unmatched) for all jets (the analyzer takes it from the matched reco jet)
//   - both V0 species in one pass (separate mixing pools, as separate analyzer jobs per species)
//
//   root -l -b -q 'genV0_analyzer.cpp+("<in.root>[,<in2.root>...]","<out prefix>"[,maxEvents,multBinned,ntrkScale])'
//   p-Pb: multBinned = true, ntrkScale = 0.844   ->  <out prefix>_Lam.root, <out prefix>_K0s.root
// Condor: ../condor/corr.sub (10 files per job), merge with ../condor/merge_corr.sh
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMath.h"

#include "../include/AnalysisConfig.h"
#include "include/CorrConfig.h"
#include "include/CorrTools.h"

namespace gana {

// copied from jetV0_analyzer.cpp (namespace ana): one multiplicity bin in the correlation sparses (MC, as the
// analyzer) unless gMultBinned (p-Pb: the data binning CorrCfg::MultBins(), classes selected in jetshape.cpp)
bool gMultBinned = false;
THnSparseD* MakeCorr(const char* name) {   // CorrCfg::Ax::kCorr*
    const std::vector<double> pt = CorrCfg::CorrPtBins(), mu = gMultBinned ? CorrCfg::MultBins() : CorrCfg::MultBinsCorrMC();
    Int_t    nb[5] = {CorrCfg::kNDPhi, CorrCfg::kNDEta, (int)pt.size() - 1, 3, (int)mu.size() - 1};
    Double_t lo[5] = {-TMath::Pi() / 2, -CorrCfg::kDEtaMax, pt.front(), -0.5, mu.front()};
    Double_t hi[5] = {3 * TMath::Pi() / 2, CorrCfg::kDEtaMax, pt.back(), 2.5, mu.back()};
    auto* h = new THnSparseD(name, Form("%s;#Delta#phi;#Delta#eta;V0 p_{T};flavor;N_{trk}^{offline}", name), 5, nb, lo, hi);
    h->GetAxis(CorrCfg::Ax::kCorrV0Pt)->Set(pt.size() - 1, pt.data());
    h->GetAxis(CorrCfg::Ax::kCorrMult)->Set(mu.size() - 1, mu.data());
    h->Sumw2();
    return h;
}
THnSparseD* MakeJet(const char* name) {    // CorrCfg::Ax::kJet*: pT, eta, flavor, mult
    const std::vector<double> mu = CorrCfg::MultBins();
    Int_t    nb[4] = {600, 40, 3, (int)mu.size() - 1};
    Double_t lo[4] = {0., -2., -0.5, mu.front()};
    Double_t hi[4] = {3000., 2., 2.5, mu.back()};
    auto* h = new THnSparseD(name, Form("%s;p_{T};#eta;flavor;N_{trk}^{offline}", name), 4, nb, lo, hi);
    h->GetAxis(CorrCfg::Ax::kJetMult)->Set(mu.size() - 1, mu.data());
    h->Sumw2();
    return h;
}

// one V0 species: gen jets x gen V0, same-event and mixed-event (as CorrSet "gen" of the analyzer)
struct Species {
    std::string v0;
    Corr::CorrHists same, mix;
    THnSparseD* hGenJet = nullptr;
    TH1D *hNev = nullptr, *hMult = nullptr;
    TH2D* hPthat = nullptr;
    Corr::MixPool* pool = nullptr;
    std::vector<Corr::V0c> vGen;
    void Make(const std::string& name) {
        v0 = name;
        same.corr   = MakeCorr("hist_correlation_signal_jet_gen_V0_gen");
        mix.corr    = MakeCorr("hist_correlation_mixing_jet_gen_V0_gen");
        same.corrPt = MakeCorr("hist_correlation_signal_jet_gen_V0_gen_ptw");
        mix.corrPt  = MakeCorr("hist_correlation_mixing_jet_gen_V0_gen_ptw");
        same.jets   = MakeJet("hist_jet_from_gen_gen_sig");
        mix.jets    = MakeJet("hist_jet_from_gen_gen_mix");
        hGenJet     = MakeJet("hist_gen_jet_weighted");
        hNev   = new TH1D("Nevents", "events;;", 3, 0, 3);
        hNev->GetXaxis()->SetBinLabel(1, "all"); hNev->GetXaxis()->SetBinLabel(2, "sum of weights"); hNev->GetXaxis()->SetBinLabel(3, "with gen jet");
        hMult  = new TH1D("multiplicity_weighted", ";N_{trk} (gen);events", 500, 0, 500);
        hPthat = new TH2D("pthathist_weighted", ";#hat{p}_{T};N_{trk} (gen)", 200, 0, 1000, 50, 0, 500);
        const int nMix = CorrCfg::NMix("MC");
        pool = new Corr::MixPool(nMix, CorrCfg::kMixDVz, CorrCfg::kMixDMult, Cfg::kVzMax, CorrCfg::MultBins().back(), CorrCfg::MixDepth(nMix));
    }
    void Write(const TString& file) {
        TFile fo(file, "RECREATE");
        fo.mkdir("QA_histograms")->cd();
        hNev->Write(); hMult->Write(); hPthat->Write(); hGenJet->Write();
        fo.mkdir("correlation_gen_gen_histograms")->cd();
        for (auto* h : {same.corr, mix.corr, same.corrPt, mix.corrPt, same.jets, mix.jets}) h->Write();
        fo.Close();
        std::cout << "Saved " << file << std::endl;
    }
};

} // namespace gana

// multBinned: correlation sparses with the data multiplicity binning (p-Pb); ntrkScale: multiplicity = ntrkScale x
// gen N_trk (hiBin). 0.844 converts the generator-level count to the reconstructed N_trk^offline (mean ratio in the
// official PYTHIA skims, flat in N_trk, rms 6%), so that the data class edges (60, 80, 185, 250) apply directly.
void genV0_analyzer(TString inFile, TString outPrefix, Long64_t maxEvents = -1, bool multBinned = false, double ntrkScale = 1.0) {
    using namespace gana;
    gMultBinned = multBinned;
    std::cout << "multiplicity = " << ntrkScale << " x gen N_trk, correlation sparses " << (multBinned ? "in N_trk classes" : "multiplicity-integrated") << std::endl;
    TStopwatch sw; sw.Start();
    TH1::AddDirectory(kFALSE);

    // inFile: one file or a comma-separated list (one mixing pool over all of them, as one analyzer job)
    std::vector<TString> files;
    { TObjArray* tok = inFile.Tokenize(","); for (auto o : *tok) files.push_back(((TObjString*)o)->GetString()); delete tok; }
    TChain* t = new TChain("hiEvtAnalyzer/HiTree");
    std::vector<TChain*> fr;
    for (const char* n : {"ak4PFJetAnalyzer/t", "HiGenParticleAna/hi", "gK0sTree/gK0sTree", "gLamTree/gLamTree"}) {
        fr.push_back(new TChain(n));
        for (const auto& f : files) fr.back()->Add(f);
    }
    for (const auto& f : files) t->Add(f);
    for (auto* c : fr) {
        if (c->GetEntries() != t->GetEntries()) { std::cerr << "FATAL: " << c->GetName() << " entries differ" << std::endl; return; }
        t->AddFriend(c);
    }
    std::cout << files.size() << " file(s), " << t->GetEntries() << " events" << std::endl;
    t->SetBranchStatus("*", 0);
    Float_t vz = 0, pthat = 0, weight = 0; Int_t hiBin = 0, ngen = 0;
    static Float_t genpt[9999], geneta[9999], genphi[9999];
    std::vector<double> *gpPt = nullptr, *gpEta = nullptr, *gpPhi = nullptr;
    std::vector<int>* gpPdg = nullptr;
    std::vector<double> *kPt = nullptr, *kEta = nullptr, *kPhi = nullptr, *lPt = nullptr, *lEta = nullptr, *lPhi = nullptr;
    auto br = [&](const char* n, void* a) { t->SetBranchStatus(n, 1); if (t->SetBranchAddress(n, a) < 0) std::cerr << "FATAL: branch " << n << std::endl; };
    br("vz", &vz); br("pthat", &pthat); br("weight", &weight); br("hiBin", &hiBin);
    br("ngen", &ngen); br("genpt", genpt); br("WTAgeneta", geneta); br("WTAgenphi", genphi);   // forest WTA axis = seed, as in the analyzer
    br("pt", &gpPt); br("eta", &gpEta); br("phi", &gpPhi); br("pdg", &gpPdg);
    br("gK0s_pt", &kPt); br("gK0s_eta", &kEta); br("gK0s_phi", &kPhi);
    br("gLam_pt", &lPt); br("gLam_eta", &lEta); br("gLam_phi", &lPhi);

    Species sp[2];
    sp[0].Make("Lam"); sp[1].Make("K0s");
    TRandom3 rnd(12345);   // as the analyzer
    std::vector<Corr::Jet> jGen;

    const Long64_t n = maxEvents > 0 ? std::min(maxEvents, t->GetEntries()) : t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        const double w = weight;
        const double mult = ntrkScale * hiBin;
        // ---- gen V0 (as the analyzer: |eta| < kV0EtaMax, kV0PtMin <= pT < kV0PtMax; correlated from kCorrPtMin)
        for (int s = 0; s < 2; ++s) {
            auto& v = sp[s].vGen; v.clear();
            const std::vector<double>& pt = *(s == 0 ? lPt : kPt), &eta = *(s == 0 ? lEta : kEta), &phi = *(s == 0 ? lPhi : kPhi);
            for (size_t k = 0; k < pt.size(); ++k) {
                if (std::fabs(eta[k]) >= Cfg::kV0EtaMax || pt[k] < Cfg::kV0PtMin || pt[k] >= Cfg::kV0PtMax) continue;
                if (pt[k] >= CorrCfg::kCorrPtMin) v.push_back({pt[k], eta[k], phi[k], 1.});
            }
            sp[s].hNev->Fill(0.5); sp[s].hNev->Fill(1.5, w);
            sp[s].hMult->Fill(mult, w); sp[s].hPthat->Fill(pthat, mult, w);
        }
        // ---- gen jets: V0-free WTA axis (as the analyzer)
        jGen.clear();
        for (int j = 0; j < ngen; ++j) {
            const double pt = genpt[j];
            if (pt <= CorrCfg::kJetPtMin || pt >= CorrCfg::kJetPtMax) continue;
            double aEta, aPhi;
            if (!Corr::ReclusterWTA(geneta[j], genphi[j], CorrCfg::kJetR, *gpPt, *gpEta, *gpPhi, *gpPdg, aEta, aPhi)) continue;
            if (std::fabs(aEta) >= CorrCfg::kJetEtaMax) continue;
            jGen.push_back({pt, aEta, aPhi, 2});
        }
        // ---- same event, mixing (before this event's V0s enter the pool), pool update
        for (int s = 0; s < 2; ++s) {
            for (const auto& j : jGen) { double x[4] = {j.pt, j.eta, (double)j.flav, mult}; sp[s].hGenJet->Fill(x, w); }
            if (!jGen.empty()) sp[s].hNev->Fill(2.5, w);
            Corr::FillSame(jGen, sp[s].vGen, w, mult, sp[s].same);
            sp[s].pool->Mix(vz, mult, w, jGen, sp[s].mix, rnd);
            sp[s].pool->Add(vz, mult, sp[s].vGen);
        }
        if (i % 20000 == 0) std::cout << i << " / " << n << " events, " << sw.RealTime() << " s" << std::endl, sw.Continue();
    }
    for (int s = 0; s < 2; ++s) {
        sp[s].pool->Summary(sp[s].v0.c_str(), sp[s].mix);
        sp[s].Write(outPrefix + "_" + sp[s].v0.c_str() + ".root");
    }
    std::cout << "Done: " << n << " events, " << sw.RealTime() << " s" << std::endl;
}
