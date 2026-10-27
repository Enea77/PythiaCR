// injet_eff.h -- in-jet V0 efficiency correction (relative to the inclusive EPOS efficiency), from PYTHIA.
//
// The EPOS efficiency eps(pT, eta) is measured for inclusive V0s. Near the jet axis the V0 reconstruction efficiency is
// lower (dense environment). In PYTHIA (correlations/jetV0_analyzer.cpp, efficiency_QA/h2_*_dr):
//   c(pT, dr) = eps_PYTHIA(pT, dr) / eps_PYTHIA(pT, all V0),   eps = (gen-matched reco sig - gen-matched reco bkg) / gen
// with dr the distance of the V0 to the nearest gen jet axis (pT > 60 GeV). The correction is the RATIO only: the absolute
// efficiency stays the EPOS one; c carries the relative loss near jets, assumed to be a detector effect that transfers.
// jetshape.cpp divides the acceptance- and mass-sideband-corrected correlation (signal/mixing - sideband/mixing), before the
// UE subtraction, by c(pT bin, dr = sqrt(deta^2 + dphi^2)) of each (deta, dphi) cell: the efficiency loss applies to all
// V0s near the jet axis, jet and UE alike, so it must be corrected before the UE (from large |deta|, c = 1) is subtracted.
// c = 1 for dr >= the last edge of kInJetDr.
#ifndef INJET_EFF_H
#define INJET_EFF_H
#include "TFile.h"
#include "TH2.h"
#include <vector>
#include <cmath>
#include <iostream>

namespace InJet {
const std::vector<double> kInJetDr = {0.0, 0.1, 0.2, 0.4};   // dr bins of the correction (coarse: statistics of jet V0s)
// systematic variation of the factors (set by jetshape.cpp from the suffix of its injetEff argument):
//   0 nominal, +1 / -1: c +- its statistical uncertainty (all bins together), 2: half correction, c -> (1 + c)/2
//   3 fit: c(pT) from a weighted straight-line fit over kFitPt per dr bin, evaluated at the centre of the requested pT bin
//   4 / 5 fitup / fitdown: fit +- its 1 sigma band (covariance scaled by chi2/ndf when > 1)
inline int gVar = 0;
const std::vector<double> kFitPt = {1.5, 2, 3, 5, 7, 10};   // pT bins of the binned c used in the fit

struct Maps { TH2 *gen = nullptr, *sig = nullptr, *bkg = nullptr; };
inline Maps Load(const TString& file) {
    Maps m;
    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie()) { std::cout << "InJet: cannot open " << file << std::endl; return m; }
    auto get = [&](const char* n) -> TH2* {
        TH2* h = (TH2*)f->Get(n);
        if (!h) h = (TH2*)f->Get(Form("efficiency_QA/%s", n));
        if (!h) { std::cout << "InJet: missing " << n << " in " << file << std::endl; return nullptr; }
        h = (TH2*)h->Clone(Form("injet_%s", n)); h->SetDirectory(0); return h;
    };
    m.gen = get("h2_gen_dr"); m.sig = get("h2_reco_matched_sig_dr"); m.bkg = get("h2_reco_matched_bkg_dr");
    f->Close();
    return m;
}
// efficiency (S - B) / G in pT [ptlo, pthi) (bin centres) and dr [drlo, drhi); err: statistical (treated as independent)
inline double Eff(const Maps& m, double ptlo, double pthi, double drlo, double drhi, double* err = nullptr) {
    double G = 0, S = 0, B = 0, eG = 0, eS = 0, eB = 0;
    for (int ix = 1; ix <= m.gen->GetNbinsX(); ++ix) {
        const double x = m.gen->GetXaxis()->GetBinCenter(ix);
        if (x < ptlo || x >= pthi) continue;
        for (int iy = 1; iy <= m.gen->GetNbinsY(); ++iy) {
            const double y = m.gen->GetYaxis()->GetBinCenter(iy);
            if (y < drlo || y >= drhi) continue;
            G += m.gen->GetBinContent(ix, iy); eG += std::pow(m.gen->GetBinError(ix, iy), 2);
            S += m.sig->GetBinContent(ix, iy); eS += std::pow(m.sig->GetBinError(ix, iy), 2);
            B += m.bkg->GetBinContent(ix, iy); eB += std::pow(m.bkg->GetBinError(ix, iy), 2);
        }
    }
    if (G <= 0) { if (err) *err = 0; return 0; }
    const double e = (S - B) / G;
    if (err) *err = std::sqrt((eS + eB) / std::pow(G, 2) + std::pow(e, 2) * eG / std::pow(G, 2));
    return e;
}
inline std::vector<double> Factors(const Maps& m, double ptlo, double pthi, std::vector<double>* errs = nullptr);
// linear fit c = a + b (pT - x0) of the binned factors (kFitPt) for each dr bin; value (and 1 sigma band) at pT
inline std::vector<double> FitFactors(const Maps& m, double pt, std::vector<double>* band = nullptr) {
    const int keep = gVar; gVar = 0;
    const size_t nd = kInJetDr.size() - 1;
    std::vector<std::vector<double>> c(nd), e(nd); std::vector<double> x;
    for (size_t j = 0; j + 1 < kFitPt.size(); ++j) {
        std::vector<double> ej; std::vector<double> cj = Factors(m, kFitPt[j], kFitPt[j + 1], &ej);
        x.push_back(0.5 * (kFitPt[j] + kFitPt[j + 1]));
        for (size_t k = 0; k < nd; ++k) { c[k].push_back(cj[k]); e[k].push_back(ej[k]); }
    }
    gVar = keep;
    std::vector<double> out; if (band) band->clear();
    for (size_t k = 0; k < nd; ++k) {
        double S = 0, Sx = 0, Sxx = 0, Sy = 0, Sxy = 0;
        for (size_t j = 0; j < x.size(); ++j) { const double w = e[k][j] > 0 ? 1. / (e[k][j] * e[k][j]) : 0;
            S += w; Sx += w * x[j]; Sxx += w * x[j] * x[j]; Sy += w * c[k][j]; Sxy += w * x[j] * c[k][j]; }
        const double D = S * Sxx - Sx * Sx, b = (S * Sxy - Sx * Sy) / D, a = (Sy - b * Sx) / S;
        double chi2 = 0; for (size_t j = 0; j < x.size(); ++j) if (e[k][j] > 0) chi2 += std::pow((c[k][j] - a - b * x[j]) / e[k][j], 2);
        const double sc = std::max(1., chi2 / (x.size() - 2.));
        const double vaa = sc * Sxx / D, vbb = sc * S / D, vab = -sc * Sx / D;   // covariance of (a, b) for c = a + b x
        const double v = a + b * pt, ev = std::sqrt(vaa + pt * pt * vbb + 2 * pt * vab);
        double cv = v; if (gVar == 4) cv += ev; else if (gVar == 5) cv -= ev;
        out.push_back(cv); if (band) band->push_back(ev);
    }
    return out;
}
// c(pT bin, dr bin k) = eps(dr bin) / eps(all); one value per kInJetDr bin (gVar 3-5: from the linear fit in pT)
inline std::vector<double> Factors(const Maps& m, double ptlo, double pthi, std::vector<double>* errs) {
    if (gVar >= 3) return FitFactors(m, 0.5 * (ptlo + pthi), errs);
    std::vector<double> c;
    if (errs) errs->clear();
    double ea = 0;
    const double all = Eff(m, ptlo, pthi, 0, 1e9, &ea);
    for (size_t k = 0; k + 1 < kInJetDr.size(); ++k) {
        double e = 0;
        const double v = Eff(m, ptlo, pthi, kInJetDr[k], kInJetDr[k + 1], &e);
        double ck = all > 0 && v > 0 ? v / all : 1.;
        const double ek = all > 0 && v > 0 ? std::sqrt(std::pow(e / all, 2) + std::pow(v * ea / (all * all), 2)) : 0.;
        if (gVar == 1) ck += ek; else if (gVar == -1) ck -= ek; else if (gVar == 2) ck = 0.5 * (1 + ck);
        c.push_back(ck);
        if (errs) errs->push_back(ek);
    }
    return c;
}
inline double FactorAt(const std::vector<double>& c, double dr) {
    for (size_t k = 0; k + 1 < kInJetDr.size(); ++k) if (dr >= kInJetDr[k] && dr < kInJetDr[k + 1]) return c[k];
    return 1.;
}
// divides each (dphi = x, deta = y) cell by c(dr); returns the number of modified cells
inline int Apply(TH2* h, const std::vector<double>& c) {
    int n = 0;
    for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
        const double dphi = h->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
            const double deta = h->GetYaxis()->GetBinCenter(iy);
            const double f = FactorAt(c, std::sqrt(dphi * dphi + deta * deta));
            if (f == 1. || f <= 0) continue;
            h->SetBinContent(ix, iy, h->GetBinContent(ix, iy) / f);
            h->SetBinError(ix, iy, h->GetBinError(ix, iy) / f);
            ++n;
        }
    }
    return n;
}
} // namespace InJet
#endif
