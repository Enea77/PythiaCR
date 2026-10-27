// CorrConfig.h
// Settings of the jet-V0 correlation analyzer (data and MC). V0 selection, mass windows, vz cut and
// efficiency definition come from ../../include/AnalysisConfig.h, shared with the EPOS efficiency chain.
#ifndef CORR_CONFIG_H
#define CORR_CONFIG_H

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>

namespace CorrCfg {

// ---------------------------------------------------------------- events
// Event filters of the 2016 pPb forests (applied in data AND MC)
inline std::vector<std::string> EventFilters() {
    return {"pBeamScrapingFilter", "pPAprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose",
            "phfCoincFilter", "pVertexFilterCutdz1p0"};
}
// Multiplicity: hiBin of the skims, which is filled with N_trk^offline (get_Ntrkoff) in the skimmer.
// Edges allow MB classes 0-20, 20-40, 40-80, 80-185 or 0-20, 20-30, 30-40, 40-60, 60-185 (about equal fractions of the
// MB events, see inclusive/incl_summary.C; the classes are chosen in the post-processing),
// HM185 185-250, HM250 >= 250 (each HM sample is used only in its own class).
inline std::vector<double> MultBins() { return {0., 20., 30., 40., 60., 80., 185., 250., 1000.}; }
// PYTHIA correlation sparses: one multiplicity bin (the closure is multiplicity-integrated; this keeps
// the number of filled bins, the file size and the hadd time down). The axis is kept so the axis indices
// (CorrCfg::Ax) are the same in data and MC. Jet and V0 QA sparses keep MultBins() in MC.
inline std::vector<double> MultBinsCorrMC() { return {0., 10000.}; }

// ---------------------------------------------------------------- jets
constexpr double kJetPtMin  = 60.;     // corrected jet pT (reco) / gen-jet pT (gen)
constexpr double kJetPtMax  = 8160.;
constexpr double kJetEtaMax = 2.0;     // |eta_lab| of the WTA axis
constexpr double kTrkMaxOverRawMin = 0.01;   // jet ID: 0.01 < trackMax/rawpT < 0.98
constexpr double kTrkMaxOverRawMax = 0.98;
constexpr double kJetR = 0.4;
// Jets with a reco V0 candidate (signal window) in the cone have a lower response than the average jet the JEC is
// derived for (PYTHIA, ref pT 60-170 GeV: <reco/ref> 0.976 for Lambda and 0.982 for K0S jets, vs 0.993 for all jets; roughly
// flat in V0 pT). Their corrected pT is scaled by (1 + V0JetShift) before the jet pT cut, in data and MC, with the shift
// chosen to bring their mean response to that of all jets: 0.993/0.976 and 0.993/0.982 (postprocessing/jes_v0.C).
// (A first test with 2.5% / 1.5% overshot the inclusive response by 0.8% / 0.4%.)
inline double V0JetShift(const std::string& v0) { return v0 == "Lam" ? 0.017 : 0.011; }
// JEC files (in jec/), applied in sequence, per beam direction:
//   data: L2Relative of the embedded PYTHIA+EPOS sample + pp 2016 Run H L2L3Residual (JetMET recommendation for PARun2016C)
//   MC  : L2Relative of the standalone PYTHIA sample, the sample used for the closure (derived on that same sample)
inline std::vector<std::string> JECFiles(const std::string& direction, bool isMC) {
    const std::string d = direction == "pPb" ? "pgoing" : "Pbgoing";
    if (isMC) return {"Autumn16_HI_pPb_" + d + "_Unembedded_MC_L2Relative_AK4PF.txt"};
    return {"Autumn16_HI_pPb_" + d + "_Embedded_MC_L2Relative_AK4PF.txt", "Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PF.txt"};
}

// ---------------------------------------------------------------- correlations
// V0 pT bins of the correlation histograms. V0s with 1 < pT < 1.5 GeV are selected (QA) but not
// correlated: the Lambda efficiency there is too small to be reliable.
inline std::vector<double> CorrPtBins() { return {1.5, 2., 3., 4., 5., 6., 7., 8., 10.}; }
constexpr double kCorrPtMin = 1.5;

// Correlation map binning: 300 x 400 bins, i.e. about 0.021 x 0.020 in (Delta phi, Delta eta)
constexpr int    kNDPhi = 300;
constexpr int    kNDEta = 400;
constexpr double kDEtaMax = 4.0;

// ---------------------------------------------------------------- histogram axes (use these in the post-processing)
namespace Ax {
    // hist_correlation_*: Delta phi, Delta eta, V0 pT, jet flavor, multiplicity
    constexpr int kCorrDPhi = 0, kCorrDEta = 1, kCorrV0Pt = 2, kCorrFlav = 3, kCorrMult = 4;
    // hist_jet_from_*, hist_reco_jet_corr_weighted, hist_gen_jet_weighted: pT, eta, flavor, multiplicity
    constexpr int kJetPt = 0, kJetEta = 1, kJetFlav = 2, kJetMult = 3;
    // hist_reco_V0*: pT, eta, mass, multiplicity
    constexpr int kV0Pt = 0, kV0Eta = 1, kV0Mass = 2, kV0Mult = 3;
}
// jet flavor codes (refparton_flavorForB): 0 quark (|id| 1-6), 1 gluon, 2 other/unmatched (always 2 in data)
inline int Flavor(int partonId) {
    const int a = partonId < 0 ? -partonId : partonId;
    return (a >= 1 && a <= 6) ? 0 : (a == 21 ? 1 : 2);
}

// ---------------------------------------------------------------- JES tree (MC)
// One entry per reco jet (jet ID, |eta| < 2) matched to a gen jet with ref pT > kJesRefPtMin. No reco pT
// cut, so the response at ref pT ~ 60 GeV is not biased by the reco threshold.
constexpr double kJesRefPtMin = 30.;

// ---------------------------------------------------------------- event mixing
// Rolling pools: the V0s of the last MixDepth(N) events of each (vz, multiplicity) class are kept.
// A jet event is mixed with up to N pooled events with |Delta vz| <= kMixDVz and |Delta mult| <= kMixDMult
// (the class of the event and its neighbours are searched), before its own V0s enter the pool.
constexpr double kMixDVz   = 0.5;   // |Delta vz| (cm)
constexpr int    kMixDMult = 100;   // |Delta N_trk^offline|
inline int MixDepth(int nMix) { return nMix * 3 > 100 ? nMix * 3 : 100; }
inline int NMix(const std::string& sample) {   // number of events to mix per jet event
    if (sample == "HM250") return 50;
    if (sample == "HM185") return 25;
    if (sample == "MB")    return 10;
    return 20;                                  // MC
}

// ---------------------------------------------------------------- PYTHIA pthat weights
// Cross section above each generated pthat threshold, relative to pthat > 15 GeV
// (same numbers as weights.h / MC_unembedded.sub). The event weight of a job is
//   w = sigma_rel(pthat_gen of the sample) / N_events(sample),
// applied to the events with pthat in the job's window (pthatMin, pthatMax]. This also covers
// jobs that fill a window from a sample generated with a lower threshold (Pbp 120-170 from pthat80).
inline double SigmaRel(int pthatGen) {
    static const std::map<int, double> s = {
        {15, 1.0000000}, {30, 7.4232802e-2}, {50, 9.5408304e-3}, {80, 1.2977183e-3}, {120, 2.1848571e-4},
        {170, 4.0848861e-5}, {220, 1.1474575e-5}, {280, 3.3198812e-6}, {370, 7.3406939e-7},
        {460, 2.0944659e-7}, {540, 7.9191586e-8}};
    auto it = s.find(pthatGen);
    if (it == s.end()) { std::cerr << "FATAL: no cross section for pthat " << pthatGen << std::endl; return -1; }
    return it->second;
}

} // namespace CorrCfg

#endif
