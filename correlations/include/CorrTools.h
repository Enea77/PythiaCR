// CorrTools.h
// Helpers of the jet-V0 correlation analyzer: kinematics, efficiency map, vz weights,
// V0-free WTA reclustering of gen jets, correlation filling and event mixing.
#ifndef CORR_TOOLS_H
#define CORR_TOOLS_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <deque>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TRandom3.h"
#include "../../include/AnalysisConfig.h"

namespace Corr {

// ------------------------------------------------------------------ kinematics
inline double DPhi(double a, double b) {                 // (-pi, pi]; no loops (a non-finite input returns NaN)
    double d = std::remainder(a - b, 2 * TMath::Pi());
    if (d <= -TMath::Pi()) d += 2 * TMath::Pi();
    return d;
}
inline double DPhi2PC(double jet, double v0) {           // [-pi/2, 3pi/2)
    double d = DPhi(jet, v0);
    if (d < -0.5 * TMath::Pi()) d += 2 * TMath::Pi();
    return d;
}
inline double DR(double eta1, double phi1, double eta2, double phi2) {
    const double de = eta1 - eta2, dp = DPhi(phi1, phi2);
    return std::sqrt(de * de + dp * dp);
}

struct Jet { double pt, eta, phi; int flav; };
struct V0c { double pt, eta, phi, w; };

// ------------------------------------------------------------------ efficiency map (EPOS, binned)
class EffMap {
public:
    bool Load(const std::string& file, const std::string& hist) {
        TFile* f = TFile::Open(file.c_str(), "READ");
        if (!f || f->IsZombie()) { std::cerr << "FATAL: cannot open efficiency file " << file << std::endl; return false; }
        TH2D* h = (TH2D*)f->Get(hist.c_str());
        if (!h) { std::cerr << "FATAL: " << hist << " not found in " << file << std::endl; return false; }
        h_ = (TH2D*)h->Clone("epos_eff_map"); h_->SetDirectory(0); f->Close();
        std::cout << "Efficiency map: " << file << " : " << hist << std::endl;
        return true;
    }
    // bin content of the (pT, eta) bin, no interpolation; returns -1 outside the map
    double Eff(double pt, double eta) const {
        const int ix = h_->GetXaxis()->FindFixBin(pt), iy = h_->GetYaxis()->FindFixBin(eta);
        if (ix < 1 || ix > h_->GetNbinsX() || iy < 1 || iy > h_->GetNbinsY()) return -1;
        return h_->GetBinContent(ix, iy);
    }
    // 1/eps, or 0 for eps <= Cfg::kEffMin (candidate dropped)
    double Weight(double pt, double eta) const {
        const double e = Eff(pt, eta);
        return e > Cfg::kEffMin ? 1. / e : 0.;
    }
    TH2D* Map() const { return h_; }
private:
    TH2D* h_ = nullptr;
};

// ------------------------------------------------------------------ vz weights (MC)
class VzWeights {
public:
    bool Load(const std::string& file, bool pgoing) {
        TFile* f = TFile::Open(file.c_str(), "READ");
        if (!f || f->IsZombie()) { std::cerr << "FATAL: cannot open vz weight file " << file << std::endl; return false; }
        const char* name = pgoing ? "Vz_weight_Pythia_pgoing" : "Vz_weight_Pythia_Pbgoing";
        TH1D* h = (TH1D*)f->Get(name);
        if (!h) { std::cerr << "FATAL: " << name << " not found in " << file << std::endl; return false; }
        h_ = (TH1D*)h->Clone("vz_weights"); h_->SetDirectory(0); f->Close();
        std::cout << "Vz weights: " << file << " : " << name << std::endl;
        return true;
    }
    double W(double vz) const {
        const int b = h_->FindFixBin(vz);
        if (b < 1 || b > h_->GetNbinsX()) return 1.;
        const double w = h_->GetBinContent(b);
        return w > 0 ? w : 1.;
    }
private:
    TH1D* h_ = nullptr;
};

// ------------------------------------------------------------------ V0-free WTA reclustering of gen jets
// Same result as ReclusterGenJets2.c (CATree.h, BuildCATree(nodes, p = -1, WTAScheme)): all gen particles
// within dR < R of the forest gen-jet axis, except K0S, Lambda and neutrinos, are clustered pairwise with distance
// dR^2 * min(pT_i^-2, pT_j^-2); at each merge the axis of the harder branch is kept (WTA).
// Implemented with nearest-neighbour bookkeeping, O(n^2) instead of the O(n^3) loop of CATree.h
// (checked identical to the brute-force loop). Returns false if no particle is found.
inline bool ReclusterWTA(double seedEta, double seedPhi, double R,
                             const std::vector<double>& pt, const std::vector<double>& eta, const std::vector<double>& phi,
                             const std::vector<int>& pdg, double& axisEta, double& axisPhi, int* nIn = nullptr) {
    static std::vector<double> P, E, F, D; static std::vector<int> NN; static std::vector<char> A;
    P.clear(); E.clear(); F.clear();
    for (size_t k = 0; k < pt.size(); ++k) {
        const int a = std::abs(pdg[k]);
        if (a == 310 || a == 3122) continue;                  // V0-free axis
        if (a == 12 || a == 14 || a == 16) continue;          // neutrinos are not part of gen jets
        if (DR(seedEta, seedPhi, eta[k], phi[k]) > R) continue;
        P.push_back(pt[k]); E.push_back(eta[k]); F.push_back(phi[k]);
    }
    const int n = P.size();
    if (nIn) *nIn = n;
    if (n == 0) return false;
    auto dist = [&](int i, int j) {
        const double dr = DR(E[i], F[i], E[j], F[j]);
        return dr * dr * std::min(1. / (P[i] * P[i]), 1. / (P[j] * P[j]));
    };
    D.assign(n, -1); NN.assign(n, -1); A.assign(n, 1);
    auto findNN = [&](int i) {
        D[i] = -1; NN[i] = -1;
        for (int j = 0; j < n; ++j) if (j != i && A[j]) { const double d = dist(i, j); if (D[i] < 0 || d < D[i]) { D[i] = d; NN[i] = j; } }
    };
    for (int i = 0; i < n; ++i) findNN(i);
    for (int left = n; left > 1; --left) {
        int bi = -1;
        for (int i = 0; i < n; ++i) if (A[i] && (bi < 0 || D[i] < D[bi])) bi = i;
        const int L = std::min(bi, NN[bi]), Rr = std::max(bi, NN[bi]);
        const double pl = P[L], pr = P[Rr];
        // merged branch keeps the axis of the harder one (moved = axis jumps to R's position) and gains pT, so
        // its distance to every particle at that position only decreases; only neighbours of the lost position
        // need a full search. Distances between other pairs do not change.
        const bool equal = (pl == pr), moved = (pl < pr);
        if (moved)      { E[L] = E[Rr]; F[L] = F[Rr]; }
        else if (equal) { E[L] += E[Rr]; F[L] += F[Rr]; }   // as CATree.h
        P[L] = pl + pr; A[Rr] = 0;
        const int kept = moved ? Rr : L, lost = moved ? L : Rr;
        for (int i = 0; i < n; ++i) {
            if (!A[i] || i == L) continue;
            if (equal ? (NN[i] == L || NN[i] == Rr) : NN[i] == lost) findNN(i);
            else if (NN[i] == kept) { D[i] = dist(i, L); NN[i] = L; }
            else { const double d = dist(i, L); if (d < D[i]) { D[i] = d; NN[i] = L; } }
        }
        findNN(L);
    }
    for (int i = 0; i < n; ++i) if (A[i]) { axisEta = E[i]; axisPhi = F[i]; }
    return true;
}

// ------------------------------------------------------------------ correlation filling
// Correlation sparses (axes CorrCfg::Ax::kCorr*): Delta phi, Delta eta, V0 pT, jet flavor, multiplicity.
// Each pair enters twice: corr with weight w_evt / eps (counts) and corrPt with weight w_evt / eps * pT_V0.
// Jet sparse (axes CorrCfg::Ax::kJet*): pT, eta, flavor, multiplicity.
struct CorrHists {
    THnSparseD *corr = nullptr, *corrPt = nullptr, *jets = nullptr;
    long nJetEvt = 0, nTrig = 0, nShort = 0;   // mixing bookkeeping (mixed-event histograms only)
};

inline void FillPairs(const std::vector<Jet>& jets, const std::vector<V0c>& v0s, double wEvt, double mult, CorrHists& h) {
    for (const auto& j : jets)
        for (const auto& v : v0s) {
            double xc[5] = {DPhi2PC(j.phi, v.phi), j.eta - v.eta, v.pt, (double)j.flav, mult};
            h.corr->Fill(xc, wEvt * v.w);
            h.corrPt->Fill(xc, wEvt * v.w * v.pt);
        }
}
inline void FillJets(const std::vector<Jet>& jets, double wEvt, double mult, CorrHists& h) {
    for (const auto& j : jets) { double xj[4] = {j.pt, j.eta, (double)j.flav, mult}; h.jets->Fill(xj, wEvt); }
}
inline void FillSame(const std::vector<Jet>& jets, const std::vector<V0c>& v0s, double wEvt, double mult, CorrHists& h) {
    FillJets(jets, wEvt, mult, h);
    FillPairs(jets, v0s, wEvt, mult, h);
}

// ------------------------------------------------------------------ event mixing (rolling pools)
// The V0s of the last `depth` events with at least one V0 are kept per (vz, multiplicity) class; classes
// are dVz wide in vz and dMult wide in multiplicity. An event with jets is mixed, before its own V0s are
// added, with up to nMix pooled events chosen at random among those with |Delta vz| <= dVz and
// |Delta mult| <= dMult (own class and neighbours). Pair weight: w_trigger / eps / n_found (single event
// weight, normalized by the number of partners actually found). Memory does not grow with the job size.
class MixPool {
public:
    MixPool(int nMix, double dVz, double dMult, double vzMax, double multMax, int depth)
        : nMix_(nMix), dVz_(dVz), dMult_(dMult), vzMax_(vzMax), depth_(depth) {
        nVz_ = (int)std::ceil(2 * vzMax / dVz); nMu_ = (int)std::ceil(multMax / dMult) + 1;
        pools_.resize(nVz_ * nMu_);
    }
    // Mix the jets of an event with up to nMix pooled events (the event's own V0s are NOT in the pool yet).
    // The same pool can be mixed with several jet collections (e.g. reco and gen jets) before Add().
    void Mix(float vz, float mult, double w, const std::vector<Jet>& jets, CorrHists& h, TRandom3& rnd) {
        if (jets.empty()) return;
        const int iv = IVz(vz), im = IMu(mult);
        cand_.clear();
        for (int a = std::max(0, iv - 1); a <= std::min(nVz_ - 1, iv + 1); ++a)
            for (int b = std::max(0, im - 1); b <= std::min(nMu_ - 1, im + 1); ++b)
                for (const auto& e : pools_[a * nMu_ + b])
                    if (std::fabs(e.vz - vz) <= dVz_ && std::fabs(e.mult - mult) <= dMult_) cand_.push_back(&e);
        const int n = std::min<int>(nMix_, cand_.size());
        for (int k = 0; k < n; ++k) std::swap(cand_[k], cand_[k + rnd.Integer(cand_.size() - k)]);   // partial shuffle
        ++h.nJetEvt;
        if (n == 0) return;
        ++h.nTrig; if (n < nMix_) ++h.nShort;
        FillJets(jets, w, mult, h);
        for (int k = 0; k < n; ++k) FillPairs(jets, cand_[k]->v0s, w / n, mult, h);
    }
    // Add the event's V0s to its (vz, multiplicity) class (after all Mix calls of this event)
    void Add(float vz, float mult, const std::vector<V0c>& v0s) {
        if (v0s.empty()) return;
        auto& p = pools_[IVz(vz) * nMu_ + IMu(mult)];
        p.push_back({vz, mult, v0s});
        if ((int)p.size() > depth_) p.pop_front();
    }
    void Process(float vz, float mult, double w, const std::vector<Jet>& jets, const std::vector<V0c>& v0s,
                 CorrHists& h, TRandom3& rnd) {
        Mix(vz, mult, w, jets, h, rnd);
        Add(vz, mult, v0s);
    }
    void Summary(const char* name, const CorrHists& h) const {
        std::cout << "  mixing " << name << ": " << h.nTrig << " of " << h.nJetEvt << " jet events mixed, " << h.nShort
                  << " with fewer than " << nMix_ << " partners" << std::endl;
    }
private:
    struct Entry { float vz, mult; std::vector<V0c> v0s; };
    int IVz(double vz) const { int i = (int)std::floor((vz + vzMax_) / dVz_); return std::min(std::max(i, 0), nVz_ - 1); }
    int IMu(double m) const  { int i = (int)std::floor(m / dMult_);            return std::min(std::max(i, 0), nMu_ - 1); }
    int nMix_; double dVz_, dMult_, vzMax_; int depth_, nVz_, nMu_;
    std::vector<std::deque<Entry>> pools_;
    std::vector<const Entry*> cand_;
};

} // namespace Corr

#endif
