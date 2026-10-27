// AnalysisConfig.h
// Single source of truth for the V0 selection used by EVERY step of the analysis
// (EPOS efficiency skim, efficiency maps, closures, data analyzers).
// If a number here changes, everything downstream must be rerun.
#ifndef ANALYSIS_CONFIG_H
#define ANALYSIS_CONFIG_H

#include <string>
#include <vector>
#include <cmath>

namespace Cfg {

// ---------------------------------------------------------------- event
constexpr double kVzMax = 15.0;   // |vz| < 15 cm (AN)

// ---------------------------------------------------------------- V0 kinematics
constexpr double kV0EtaMax = 2.0;  // |eta| < 2
constexpr double kV0PtMin  = 1.0;  // GeV
constexpr double kV0PtMax  = 10.0; // GeV

// ---------------------------------------------------------------- mass windows
// Centres and widths from the Gaussian-fit study on EPOS matched candidates
// (Plot_Lam_Sigma_Analysis / Plot_K0s_Sigma_Analysis). Signal = +-3 sigma,
// background = both sidebands 5-8 sigma. Widths in GeV.
constexpr double kMassCenterLam = 1.11585;
constexpr double kMassCenterK0s = 0.4982;
constexpr double kSigmaLam      = 0.002316;  // nominal
constexpr double kSigmaLamLow   = 0.0015;    // systematic
constexpr double kSigmaLamHigh  = 0.00335;   // systematic
constexpr double kSigmaK0s      = 0.0058;    // nominal
constexpr double kSigmaK0sLow   = 0.0040;    // systematic
constexpr double kSigmaK0sHigh  = 0.0080;    // systematic
constexpr double kNSigSignal    = 3.0;
constexpr double kNSigSideLo    = 5.0;
constexpr double kNSigSideHi    = 8.0;

inline double MassCenter(const std::string& v0) { return v0 == "Lam" ? kMassCenterLam : kMassCenterK0s; }
inline double Sigma(const std::string& v0)      { return v0 == "Lam" ? kSigmaLam : kSigmaK0s; }

// ---------------------------------------------------------------- MC matching
constexpr double kDrMatchV0 = 0.1;  // reco V0 <-> gen V0 (Lambda or K0S)
constexpr double kDrMatchXi = 0.2;  // reco Lambda <-> gen Xi (feed-down Lambda, counted in numerator)
constexpr double kDrMatchOm = 0.3;  // reco Lambda <-> gen Omega- (feed-down Lambda, counted in numerator; EPOS skim axis 6)

// ---------------------------------------------------------------- efficiency weights
constexpr double kEffMin = 0.005;   // eff <= kEffMin -> candidate dropped (weight 0)

// ---------------------------------------------------------------- topological cuts
struct CutSet {
    std::string name;
    double dxyz;       // |dxy|,|dz| of both daughters >  dxyz (significance)
    double cospoint;   // cos(3D pointing angle)       >  cospoint
    double trkdca;     // DCA between daughters        <= trkdca
    double decaylen;   // 3D decay-length significance >  decaylen
    int    nhits;      // hits on each daughter        >= nhits
    double mismass;    // mis-ID veto half-width (Lam: |m_pipi - m_K0S|; K0S: |m_ppi - m_Lam|)
};

// Nominal mis-ID veto half-widths differ by species; the "mismass" variations
// scale it (Lam: 0.015 / 0.020 / 0.025, K0S: 0.0075 / 0.010 / 0.0125).
inline double NominalMisMass(const std::string& v0) { return v0 == "Lam" ? 0.020 : 0.010; }

inline std::vector<CutSet> CutVariations(const std::string& v0) {
    const double mm = NominalMisMass(v0);
    return {
        {"nominal",        1.0, 0.999,  0.5, 5.0, 4, mm},
        {"dxyz_tight",     1.2, 0.999,  0.5, 5.0, 4, mm},
        {"dxyz_loose",     0.8, 0.999,  0.5, 5.0, 4, mm},
        {"cospoint_tight", 1.0, 0.9995, 0.5, 5.0, 4, mm},
        {"cospoint_loose", 1.0, 0.998,  0.5, 5.0, 4, mm},
        {"dca_tight",      1.0, 0.999,  0.4, 5.0, 4, mm},
        {"dca_loose",      1.0, 0.999,  0.6, 5.0, 4, mm},
        {"decay_tight",    1.0, 0.999,  0.5, 6.0, 4, mm},
        {"decay_loose",    1.0, 0.999,  0.5, 4.0, 4, mm},
        {"nhits_tight",    1.0, 0.999,  0.5, 5.0, 5, mm},
        {"nhits_loose",    1.0, 0.999,  0.5, 5.0, 3, mm},
        {"mismass_tight",  1.0, 0.999,  0.5, 5.0, 4, mm * 1.25},
        {"mismass_loose",  1.0, 0.999,  0.5, 5.0, 4, mm * 0.75},
    };
}

} // namespace Cfg

#endif
