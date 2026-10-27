# PythiaCR logbook

Running record of the standalone PYTHIA 8 study for the CMS pPb 8.16 TeV in-jet (Λ+Λ̄)/2K0S analysis:
what was done, why, the numbers obtained, and what is still open. Newest status at the top of each section,
dated entries in the journal at the bottom. How-to commands are in `README.md`.

---

## 1. Goal

The in-jet ratio measured in data (jets pT > 60 GeV, Δr < 0.4, UE subtracted with the jet–V0 correlation method)
is 0.248 ± 0.002 (MB, 3 < pT < 10 GeV), while the official CMS PYTHIA 8 dijet sample gives 0.196 ± 0.001 at
generator level. Question: can colour reconnection (CR) in PYTHIA move the in-jet ratio towards the data, as the
CR-BLC modes did for ALICE (JHEP 07 (2023) 136, Fig. 7)?

Plan:
1. **Validate the generator**: reproduce the official sample (same tune, beams, decays, jets) at generator level.
   Only then are differences between tunes meaningful.  ← *done, passed (§3)*
2. Generate CR variants with the same setup: CR-BLC modes 0/2/3 (ALICE) and the CMS tunes CP5, CP5-CR1, CP5-CR2.
   ← *done (all 7 samples, 10M events each)*
3. Run the analysis correlation method (mixing, UE from 1.5 < |Δη| < 2.5) on the gen trees with a small gen-only
   analyzer that uses `v0-jet-analysis/correlations/include/CorrTools.h`, and compare with data.
   ← *done for all 7 samples (§3c)*
4. Later: multiplicity dependence needs pPb (PYTHIA Angantyr), pp cannot reach N_trk ≥ 185.

---

## 2. Setup and the reasons behind each choice

| Choice | Reason |
|---|---|
| Standalone PYTHIA 8.313 + FastJet 3.4.3 + ROOT 6.34 from the LCG_107a view on CVMFS (`setup.sh`) | No installation, no admin rights; same software on lxplus and the batch nodes. CMSSW not needed for gen level. |
| HTCondor on lxplus (`condor/submit.sub`) | Jobs survive disconnections; parallel production. |
| Output on **user EOS** `/eos/user/e/eprifti/PythiaCR/<process>_<tune>/` | Group area `phys_heavyions` is shared and 96.8% full (1.40 of 1.50 PB, 2026-10-21); user EOS has 2 TB quota, 0.72 TB used. Production ~3.5 GB per 1M events. (The first 4 test files, ~28 MB, are in the group area under `eprifti/PythiaCR/`, can be deleted.) |
| Generator writes **TTrees**, histograms are made afterwards | Binning, cuts and observables can change without regenerating. |
| Trees use the **same directory/tree/branch names and types as the CMS V0+jet skims** (gen part) | The same readers run on the official skims and on our samples; the analysis code can be reused. |

### Physics settings (`cards/`)
Cards are read in order: `common` → process (`dijet` or `minbias`) → [`crblc_common`] → tune.

- **Beams** (`common.cmnd`): pp, 6500 GeV × 2562.5 GeV (proton energy and Pb energy per nucleon), i.e.
  √s = 8.16 TeV with the lab frame boosted by y = +0.465 (6.5 TeV proton towards +z). This is the direction of the
  official `pgoing` skims (checked with the charged-particle η distribution, §3), so |η_lab| < 2 means the same thing.
  With symmetric |η| cuts both beam directions give identical results, so only one is generated.
  (The first CUETP8M1 production had the opposite direction, y = −0.465; it is kept as `dijet_cuetp8m1_ycmNeg`.)
- **Decays**: CMS defaults `ParticleDecays:limitTau0 = on`, `tau0Max = 10` mm. K0S, Λ, Ξ, Ω, Σ± are stable
  (as in the CMS gen record, where GEANT decays them), Σ0 → Λγ decays (that Λ is primary). This matches the
  analysis definition: final Λ = prompt, no Ξ/Ω feed-down. Check: ~36% of Λ come from Σ0.
- **Process** (`dijet.cmnd`): `HardQCD:all` with a pTHat^4.5 bias from pTHatMin = 15 GeV
  (`PhaseSpace:bias2Selection`). Events are generated roughly flat in pTHat and carry a weight
  (`weight` branch, = (15/pTHat)^4.5) that restores the physical spectrum: **always fill with the weight**.
  No pTHat bins or cross-section stitching needed.
- **Tunes**:
  - `cuetp8m1.cmnd`: CMS CUETP8M1 = Monash + MPI pT0Ref 2.4024, ecmPow 0.25208, expPow 1.6
    (from `Configuration/Generator/python/Pythia8CUEP8M1Settings_cfi.py`). **Confirmed** as the tune of the official
    sample (McM request `HIN-pPb816Spring16GS-00052`, found with DAS, 2026-10-21): fragment =
    Pythia8CommonSettings + Pythia8CUEP8M1Settings + `HardQCD:all = on`, `pTHatMin = 15`, comEnergy 8160 GeV.
    PDF: Monash default NNPDF2.3 QCD+QED LO in both.
  - `monash.cmnd`: Monash 2013 (PYTHIA default).
  - `crblc_common.cmnd` + `crblc_mode{0,2,3}.cmnd`: QCD-based CR "beyond leading colour" (ColourReconnection:mode = 1,
    junctions allowed), the CR-BLC modes used by ALICE.
- **Jets** (in `src/runPythia.cc`): anti-kT R = 0.4 on all final-state particles except neutrinos, no η cut
  (like CMS gen jets), stored if pT > 30 GeV, with the E-scheme axis and a WTA axis (C/A reclustering with the
  WTA pT scheme). The V0-free WTA axis of the analysis is rebuilt later from the gen particles with
  `Corr::ReclusterWTA`, exactly as for the official sample.

### Where the CR-BLC parameter values come from
1. PYTHIA 8.313 manual (`share/Pythia8/xmldoc/ColourReconnection.xml`): for mode 1 it says to use the defaults with
   `MultipartonInteractions:pT0Ref = 2.15` and `ColourReconnection:allowDoubleJunRem = off`, and points to
   Christiansen & Skands for the full list.
2. J.R. Christiansen, P.Z. Skands, "String formation beyond leading colour", JHEP 08 (2015) 003, arXiv:1505.01681:
   the table in the appendix lists Monash, Mode 0, Mode 2, Mode 3. ALICE cites this paper ([62]) for its CR-BLC modes.
   All values in our cards were checked against that table (2026-10-21):

| Parameter | Monash | Mode 0 | Mode 2 | Mode 3 |
|---|---|---|---|---|
| StringPT:sigma | 0.335 | 0.335 | 0.335 | 0.335 |
| StringZ:aLund / bLund | 0.68 / 0.98 | 0.36 / 0.56 | 0.36 / 0.56 | 0.36 / 0.56 |
| StringFlav:probQQtoQ | 0.081 | 0.078 | 0.078 | 0.078 |
| StringFlav:probStoUD | 0.217 | 0.2 | 0.2 | 0.2 |
| StringFlav:probQQ1toQQ0join | 0.5,0.7,0.9,1.0 | 0.0275 ×4 | 0.0275 ×4 | 0.0275 ×4 |
| MultipartonInteractions:pT0Ref | 2.28 | 2.12 | 2.15 | 2.05 |
| BeamRemnants:remnantMode / saturation | 0 / – | 1 / 5 | 1 / 5 | 1 / 5 |
| ColourReconnection:mode | 0 | 1 | 1 | 1 |
| ColourReconnection:allowDoubleJunRem | on | off | off | off |
| ColourReconnection:m0 | – | 2.9 | 0.3 | 0.3 |
| ColourReconnection:junctionCorrection | – | 1.43 | 1.20 | 1.15 |
| ColourReconnection:timeDilationMode / Par | – | 0 / – | 2 / 0.18 | 3 / 0.073 |

Mode 2 values equal the PYTHIA 8.3 defaults of these parameters, so they do not appear in the "changed settings" printout.
Note: the CR-BLC tune is built on Monash (not CUETP8M1); the CR modes are used as published.

### CMS tunes
The CR-BLC modes are **not** CMS tunes (they are the Christiansen–Skands tunes, built on Monash). The CMS tunes
generated, copied line by line from `CMSSW_13_0_5/src/Configuration/Generator/python/`:

| card | CMS tune | CR model | PDF |
|---|---|---|---|
| `cuetp8m1.cmnd` | CUETP8M1 (2016; tune of the official pPb 8.16 TeV sample) | MPI-based (mode 0, default) | NNPDF2.3 LO |
| `cp5.cmnd` | CP5 (2017+ default) | MPI-based (mode 0), range 5.176 | NNPDF3.1 NNLO |
| `cp5cr1.cmnd` | CP5-CR1 (`MCTunes2017/PythiaCP5CR1TuneSettings_cfi.py`) | QCD-based (mode 1), junctionCorrection 0.238, m0 1.721, timeDilationPar 8.58 | NNPDF3.1 NNLO |
| `cp5cr2.cmnd` | CP5-CR2 (`MCTunes2017/PythiaCP5CR2TuneSettings_cfi.py`) | gluon move (mode 2), m2Lambda 4.395, fracGluon 0.9896 | NNPDF3.1 NNLO |

CP5-CR1 uses the same CR model as CR-BLC but was tuned by CMS to UE data with a junction *suppression*
(0.238 vs 1.15–1.43 in CR-BLC), so it is expected to give fewer baryons than CR-BLC.
NNPDF3.1 comes from LHAPDF6; the LCG view already sets `LHAPDF_DATA_PATH` (sets + `lhapdf.conf`), do not override it.
With NNPDF3.1 NNLO PYTHIA prints a few "SimpleSpaceShower::pT2nearThreshold: stuck in loop / xPDF = 0" errors; these
are retried internally and are normal for NNLO PDFs.

### Output trees (one file per job)
| Tree | Branches |
|---|---|
| `hiEvtAnalyzer/HiTree` | `pthat`, `weight`, `vz` (= 0), `hiBin` (= gen N_trk: charged final-state, \|η\| < 2.4, pT > 0.4), `nMPI` |
| `ak4PFJetAnalyzer/t` | `ngen`, `genpt`, `geneta`, `genphi`, `WTAgeneta`, `WTAgenphi` |
| `HiGenParticleAna/hi` | `pt`, `eta`, `phi`, `pdg`, `chg` (final-state, \|η\| < 2.4, pT > 0.4, as in the skims) |
| `gK0sTree/gK0sTree`, `gLamTree/gLamTree`, `gXiTree/gXiTree`, `gOmTree/gOmTree` | `g*_pt`, `_eta`, `_phi`, `_mass`, `_id`, `_mom1` (\|η\| < 2.4; `mom1` = PDG id of a hadron mother, 0 if from a string) |
| `hWeights` (histogram) | events, sum of weights, σ_gen × events |

Differences from the official skims to keep in mind: `hiBin` is a generator-level N_trk (the skims store the
reconstructed N_trk^offline); `vz` is 0; gen V0 trees have no pT cut (the skims have pT > 1 GeV).

---

## 3. Validation against the official sample (step 1)

**Official sample**: `Dijet_pThat-X_pPb-Bst_8p16_Pythia8` Summer16 80X, unembedded (pp at 8.16 TeV boosted, no EPOS),
pTHat thresholds 15…540 GeV. Skims: `/eos/cms/store/group/phys_heavyions/eprifti/JetV0_skims/MC/unembedded/V0JetTrkSkims_pgoing_pthat*_wXi.root`.
Only the `pgoing` files are used (the `Pbgoing` pthat120 file is 23 kB; both directions are equivalent here).
Combination as in the analysis (`CorrCfg::SigmaRel`): each sample fills its window (threshold, next threshold] with
weight σ_rel / N_events.

**Method** (`macros/genCompare.C`, same code for both samples): gen jets pT > 60, |η_WTA| < 2 (forest-style WTA axis),
gen V0s |η| < 2, jet–V0 pairs per jet, **no UE subtraction** (this is a generator check, not the physics result).
Observables: jet pT and η; N_ch and charged η, pT in events with a jet (tests MPI/UE and the boost direction);
K0S, Λ, Ξ per jet in the cone vs pT; Λ/2K0S in the cone, in the whole jet event, and vs Δr.
Scripts: `macros/runOfficial.sh` → `compare/official_pgoing.root`; `macros/runGenerated.sh <sample>` →
`compare/<sample>.root`; `macros/plotCompare.C` → plots and a summary table.

**Status: validation PASSED (2026-10-21).** Standalone CUETP8M1, 10M events (cluster 9724509, 30 CPU h, 31 GB)
vs official pgoing sample (13.5M jets > 60 GeV), generator level, no UE subtraction:

| sample | jets > 60 (entries) | ⟨N_ch⟩ in jet events | in-cone Λ/2K0S, 3–10 GeV | Λ from Σ0 (in cone) | K0S / jet | Λ / jet |
|---|---|---|---|---|---|---|
| CMS official | 13.5 M | 61.5 | 0.1946 ± 0.0015 | 0.357 | 0.2245 | 0.0874 |
| standalone CUETP8M1 | 13.7 M | 61.7 | 0.1946 ± 0.0010 | 0.362 | 0.2255 | 0.0877 |

Shapes (plots in `compare/plots_validation/`): jet pT spectrum 60–1000 GeV, jet η, N_ch distribution,
charged pT, in-cone K0S/Λ/Ξ yields vs pT, Λ/2K0S vs pT (in cone and in jet events) and vs Δr all agree within
statistics (ratios flat at 1, Ξ within its larger errors). The in-cone Λ/2K0S of 0.1946 also agrees with the
analysis' gen-level PYTHIA value (0.196 with UE subtraction), as expected since the UE contribution is small.

One difference found and fixed: the charged-particle η distribution of the first production was the mirror image of
the official one (the official `pgoing` sample has more particles at +η, i.e. the 6.5 TeV proton goes to +z).
Irrelevant for |η| < 2 cuts (all integrated numbers agree) but the beams were swapped in `common.cmnd` for all later
samples, and the mirrored sample was renamed `dijet_cuetp8m1_ycmNeg` (EOS and `compare/`).

**Conclusion**: the generator setup reproduces the official sample; differences between tunes produced with it can
be attributed to the tune.

**pTHat check (2026-10-22)**, new CUETP8M1 production (y_cm = +0.465, cluster 9724871) vs official pgoing sample
(official windows weighted with `CorrCfg::SigmaRel`; ours with the bias weight; plots `compare/plots_validation/pthat*.png`,
numbers in `compare/plots_validation/summary.txt`):

| σ(pTHat > X)/σ(pTHat > 15) | 30 | 50 | 80 | 120 | 170 | 220 | 280 | 370 | 460 | 540 |
|---|---|---|---|---|---|---|---|---|---|---|
| CorrCfg::SigmaRel (McM) | 7.423e-2 | 9.541e-3 | 1.298e-3 | 2.185e-4 | 4.085e-5 | 1.147e-5 | 3.320e-6 | 7.341e-7 | 2.094e-7 | 7.919e-8 |
| official (combined) | 7.428e-2 | 9.554e-3 | 1.305e-3 | 2.172e-4 | 4.089e-5 | 1.149e-5 | 3.323e-6 | 7.349e-7 | 2.097e-7 | 7.926e-8 |
| standalone CUETP8M1 | 7.408e-2 | 9.534e-3 | 1.298e-3 | 2.112e-4 | 4.101e-5 | 1.153e-5 | 3.343e-6 | 7.383e-7 | 2.112e-7 | 7.996e-8 |

- The pTHat spectrum agrees to better than 1% everywhere except the window 120–170 GeV, where the official combined
  spectrum is ~3.5% higher, with steps exactly at the 120 and 170 thresholds. This points to the McM cross section
  used for the pTHat-120 sample (`SigmaRel(120)` ~3.5% high relative to its neighbours), not to our generator, whose
  spectrum is continuous. Effect on the analysis: that window holds 3.4% of the jets > 60 GeV → ~0.1% on the jet
  composition; negligible for ratios, but worth knowing for absolute jet spectra.
- Events with a jet > 60 GeV, share per pTHat window (official / ours): 30–50: 21.2 ± 0.2% / 20.3 ± 0.1%;
  50–80: 55.5 / 55.9%; 80–120: 18.7 / 19.2%; 120–170: 3.43 / 3.40%; 15–30: 0.3 / 0.4%.
  Ours has ~4% (relative) fewer 60 GeV jets from pTHat 30–50, i.e. slightly less upward migration from low pTHat
  to the 60 GeV threshold. Most likely PYTHIA 8.2 (official) vs 8.313 (ours) differences in shower/MPI; it changes
  the jet sample composition by ~1%, and the in-cone Λ/2K0S still agrees (0.1960 ± 0.0010 vs 0.1946 ± 0.0015).

**Residual differences standalone vs official, understood (2026-10-27)** — plots remade with merged bins
(≤ 2% stat. error per bin) and ratio panels zoomed to the points (`macros/plotCompare.C`):

1. *pTHat 120–170 GeV, official 3.5% high* → normalisation of the official pTHat-120 sample. Overlap test: the
   same pTHat window computed from different official samples (σ_rel(gen) × N_window / N_total):

   | window | from pTHat-50 | from pTHat-80 | from pTHat-120 | from pTHat-170 | standalone |
   |---|---|---|---|---|---|
   | 120–170 | 1.702e-4 ± 0.8% | 1.696e-4 ± 0.3% | **1.761e-4** ± 0.1% | – | 1.702e-4 |
   | 170–220 | 2.99e-5 ± 1.8% | 2.933e-5 ± 0.7% | **3.037e-5** ± 0.3% | 2.938e-5 ± 0.1% | |

   Only the pTHat-120 sample is high (+3.5%, in both windows) → `CorrCfg::SigmaRel(120)` = 2.1849e-4 (from McM) is
   ~3.5% too large; the overlaps give ≈ 2.11e-4. **Consequence for the analysis**: jets from that window (3.4% of the
   jets > 60 GeV; they dominate jet pT 120–200 GeV, where the official jet spectrum is ~2% high) are over-weighted by
   3.5% → ~0.1% on the jet composition, negligible for the ratios. Could be fixed by taking σ_rel(120) from the overlap.
2. *Official jet spectrum has one bin (280–302 GeV) 4% high with 4% error* → one event: pTHat-30 sample, evt 131282,
   pTHat = 38.5 GeV with two gen jets of ~290 GeV (unphysical for the hard process: a pTHat outlier). It carries
   the large pTHat 30–50 weight (≈ 6500 × a pTHat-280 event). Effect on the ratios ~10⁻⁵ of the jet weight. Standard
   cure: reject events with leading gen-jet pT ≫ pTHat (not applied here).
3. *Remaining*: 4% fewer jets > 60 GeV from pTHat 30–50 and a slow rise of the jet-pT ratio above 300 GeV (+1–3%).
   Candidate: PYTHIA version (official 8.212 in CMSSW_8_0_X, standalone 8.313). **Test running**: the same program
   and cards built against PYTHIA 8.212 (LCG_84, CentOS 7, el7 container: `Makefile.8212`, `setup8212.sh`,
   `condor/submit8212.sub`, `run_job8212.sh`), CUETP8M1, 10M events, cluster 9747148, seeds 5001–5100, output
   `dijet_cuetp8m1_py8212`. Second-order candidate if needed: CMS uses `Tune:preferLHAPDF = 2` (LHAPDF6 copy of
   NNPDF2.3 LO instead of PYTHIA's internal grid).
4. Observables of the measurement agree within statistics already: in-cone Λ/2K0S 3–10 GeV 0.1960 ± 0.0010 vs
   0.1946 ± 0.0015 (0.8σ); through the full correlation chain 0.1955 ± 0.0010 vs 0.1964 ± 0.0012 (0.6σ);
   Λ/2K0S vs pT and vs Δr within 1–2σ bin by bin.

---

## 3b. First tune comparison (step 2, generator level, NO UE subtraction)

2026-10-22, plots `compare/plots_tunes/`, table `compare/plots_tunes_summary.txt` (CR-BLC mode 0: 98/100 files,
mode 3: 99/100 at the time; CP5 family still running):

| sample | ⟨N_ch⟩ jet events | in-cone Λ/2K0S, 3–10 GeV | K0S / jet | Λ / jet |
|---|---|---|---|---|
| CMS official (CUETP8M1, PYTHIA 8.2) | 61.5 | 0.1946 ± 0.0015 | 0.2245 | 0.0874 |
| CUETP8M1 | 61.7 | 0.1960 ± 0.0010 | 0.2242 | 0.0879 |
| CR-BLC mode 0 | 66.7 | 0.1964 ± 0.0011 | 0.1980 | 0.0778 |
| CR-BLC mode 2 | 63.4 | 0.2212 ± 0.0012 | 0.2019 | 0.0893 |
| CR-BLC mode 3 | 66.2 | 0.2109 ± 0.0012 | 0.2041 | 0.0861 |
| data (UE-subtracted, MB) | – | 0.248 ± 0.002 | – | – |

Observations:
- CR-BLC modes 2 and 3 raise the in-cone ratio by 13% and 8% (3–10 GeV); mode 0 does not change it (it lowers K0S
  and Λ per jet by the same ~12%). The enhancement is largest at low V0 pT (+40% at 1–2.5 GeV) and fades towards
  10 GeV.
- The Λ/2K0S vs Δr profile in modes 2/3 rises from ~0.17 at the axis to 0.37–0.40 at Δr ≈ 0.8, qualitatively like
  the data (0.18 → 0.35–0.40), while CUETP8M1 rises only to ~0.28.
- **Caveat**: these numbers include the UE inside the cone. The CR modes enhance baryons most outside the jet core,
  so their UE-subtracted in-cone ratio will be lower than the values above. The comparison with data (0.248) needs
  the gen-only correlation analyzer (step 3).

---

## 3c. Analysis correlation method on the generated samples (step 3)

**Code** (2026-10-22):
- Copied unchanged from `v0-jet-analysis` commit `10af763` (record in `COPIED_FROM_v0-jet-analysis.txt`):
  `include/AnalysisConfig.h`, `include/PlotStyle.h`, `correlations/include/CorrConfig.h`, `CorrTools.h`,
  `postprocessing/jetshape.cpp`, `jet_shape_V0.h`, `call_libraries.h`, `injet_eff.h`, `plot_incone.C`,
  `postprocessing/files/InJet_Efficiency_{Lam,K0s}.root`.
- New `correlations/genV0_analyzer.cpp`: the gen × gen path of `jetV0_analyzer.cpp` on the runPythia trees, with the
  analyzer's own functions (`MakeCorr`/`MakeJet` copied, `Corr::ReclusterWTA`, `FillSame`, `MixPool` with
  `NMix("MC")` = 20, `kMixDVz`, `kMixDMult`, `MixDepth`), same V0 selection (|η| < 2, 1.5 ≤ pT < 10 correlated) and jets
  (pT > 60, V0-free WTA axis seeded by the forest-style WTA axis, |η_axis| < 2). Output = `QA_histograms/hist_gen_jet_weighted`
  + `correlation_gen_gen_histograms/...`, read by `jetshape.cpp("gen","gen",...)` unchanged.
  Differences: event weight = PYTHIA weight; no filters, no vz cut (vz = 0 → mixing classes by multiplicity only);
  multiplicity = gen N_trk; jet flavour = 2 for all jets; Λ and K0S in one pass with separate mixing pools.
- Condor: `condor/corr.sub` + `corr_job.sh` + `make_corr_list.sh` (10 files = 1M events per job; outputs on user EOS
  `PythiaCR/corr/<sample>/`); `postprocessing/run_tune.sh <sample>` merges (hadd) and runs `jetshape.cpp`;
  `postprocessing/run_refs.sh` runs `jetshape.cpp` on the references; `postprocessing/plot_tunes_vs_data.C` makes
  the data vs tunes plots with the `plot_incone.C` definitions (stage `sideband_corrected`, Δr < 0.4, 0.5 Λ/K0S,
  3–10 GeV = sum of the per-bin yields).

**References** (jetshape.cpp copy, V0 pT bins 1.5, 2, 3, 5, 7, 10; inputs symlinked, outputs in `postprocessing/input_*`):

| input | in-cone Λ/2K0S, 3–10 GeV | 1.5–2 | 2–3 | 3–5 | 5–7 | 7–10 |
|---|---|---|---|---|---|---|
| official PYTHIA gen × gen (`MC_unembedded_v0shift2`) | 0.1964 ± 0.0012 | 0.103 | 0.132 | 0.168 | 0.202 | 0.237 |
| data MB, Λ `Data_v0shift_effOm` + K0S `Data_v0shift`, in-jet eff. corr. | 0.2412 ± 0.0023 | 0.154 | 0.260 | 0.247 | 0.256 | 0.220 |
| data MB, `Data_v0shift` (Λ and K0S), in-jet eff. corr. | 0.2424 ± 0.0023 | 0.155 | 0.262 | 0.248 | 0.257 | 0.221 |
| data MB, `Data_v0shift`, **no** in-jet eff. corr. | 0.2591 ± 0.0024 | 0.199 | 0.259 | 0.271 | 0.257 | 0.244 |

- The copied chain reproduces the official PYTHIA value of the manuscript (0.196 ± 0.001) exactly.
- The manuscript data value (0.248 ± 0.002) is between the corrected and uncorrected numbers: it was made with an
  earlier version of the chain. **Nominal data result (decided by the user, 2026-10-27): the first row** —
  Λ from `Data_v0shift_effOm`, K0S from `Data_v0shift`, in-jet efficiency correction: 0.2412 ± 0.0023 (stat.).
  All tune comparisons use it.

**First results (2026-10-27), full analysis chain, nominal data** (`postprocessing/plots/tunes_vs_data_prelim/`):

| sample | in-cone Λ/2K0S, 3–10 GeV | MC/data | 1.5–2 | 2–3 | 3–5 | 5–7 | 7–10 |
|---|---|---|---|---|---|---|---|
| data MB (nominal, stat. only) | 0.2412 ± 0.0023 | – | 0.154 | 0.260 | 0.247 | 0.256 | 0.220 |
| official PYTHIA (CUETP8M1, PYTHIA 8.2) | 0.1964 ± 0.0012 | 0.814 | 0.103 | 0.132 | 0.168 | 0.202 | 0.237 |
| standalone CUETP8M1 | 0.1955 ± 0.0010 | 0.810 | 0.098 | 0.134 | 0.168 | 0.201 | 0.236 |
| CR-BLC mode 2 | 0.2200 ± 0.0012 | 0.912 | 0.139 | 0.170 | 0.195 | 0.227 | 0.252 |
| CR-BLC mode 3 | 0.2104 ± 0.0012 | 0.872 | 0.127 | 0.158 | 0.186 | 0.217 | 0.243 |
| CP5 | 0.1970 ± 0.0010 | 0.817 | 0.100 | 0.139 | 0.170 | 0.202 | 0.236 |

- Standalone CUETP8M1 reproduces the official sample through the full chain as well (0.1955 vs 0.1964).
- CR-BLC mode 2 closes about half of the gap to data (MC/data 0.81 → 0.91 for 3–10 GeV). The UE subtraction
  barely changes its in-cone value (0.221 without, 0.220 with), so the extra baryons are in the jet cone itself.
- Shape vs pT: data are flat at ~0.25 from 2 to 7 GeV and drop at 7–10; all PYTHIA versions rise monotonically with
  pT and exceed the data at 7–10 GeV. CR-BLC mode 2 shifts the curve up but does not change its slope.
- Δr profile (1.5–10 GeV): CR-BLC mode 2 follows the data much better than CUETP8M1 (MC/data 0.88–0.96 vs 0.65–0.83
  for Δr < 1), but is still below the data at 0.2 < Δr < 0.6.

**All tunes, full chain (2026-10-27)** — `postprocessing/plots/tunes_vs_data/` (`InCone_Ratio`, `DeltaR_Ratio`,
`InCone_Yield_{Lam,K0s}`, `values.txt`), nominal data, PYTHIA 8.313 unless stated:

| sample | in-cone Λ/2K0S, 3–10 GeV | MC/data | 1.5–2 | 2–3 | 3–5 | 5–7 | 7–10 |
|---|---|---|---|---|---|---|---|
| data MB (stat. only) | 0.2412 ± 0.0023 | – | 0.154 | 0.260 | 0.247 | 0.256 | 0.220 |
| official PYTHIA (CUETP8M1, 8.212) | 0.1964 ± 0.0012 | 0.814 | 0.103 | 0.132 | 0.168 | 0.202 | 0.237 |
| CUETP8M1 | 0.1955 ± 0.0010 | 0.810 | 0.098 | 0.134 | 0.168 | 0.201 | 0.236 |
| CP5 | 0.1970 ± 0.0010 | 0.817 | 0.100 | 0.139 | 0.170 | 0.202 | 0.236 |
| CP5-CR1 (QCD-based CR, CMS) | 0.1804 ± 0.0010 | 0.748 | 0.112 | 0.139 | 0.160 | 0.183 | 0.209 |
| CP5-CR2 (gluon move, CMS) | 0.1854 ± 0.0009 | 0.768 | 0.097 | 0.127 | 0.159 | 0.186 | 0.226 |
| CR-BLC mode 0 | 0.1962 ± 0.0011 | 0.813 | 0.136 | 0.155 | 0.176 | 0.196 | 0.229 |
| CR-BLC mode 2 | 0.2200 ± 0.0012 | 0.912 | 0.139 | 0.170 | 0.195 | 0.227 | 0.252 |
| CR-BLC mode 3 | 0.2104 ± 0.0012 | 0.872 | 0.127 | 0.158 | 0.186 | 0.217 | 0.243 |

Conclusions so far:
- No tune reaches the data (0.241; in-jet systematic on the ratio ~3–5%). The best is CR-BLC mode 2 (−9%), then
  mode 3 (−13%); the CMS CR tunes go the wrong way (CP5-CR1 −25%, CP5-CR2 −23%: CP5-CR1 suppresses junctions,
  junctionCorrection 0.238), CP5 ≈ CUETP8M1.
- Junction (baryon) enhancement in the CR is what moves the in-jet ratio: modes 2/3 (junctionCorrection 1.20/1.15,
  m0 = 0.3, time dilation) raise it; mode 0 (m0 = 2.9, no time dilation) only raises the lowest pT bin.
- pT shape: all tunes rise monotonically with pT and exceed the data at 7–10 GeV, while the data are flat at ~0.25
  between 2 and 7 GeV; no tune reproduces the 2–3 GeV value (MC/data 0.5–0.65).
- Δr: all tunes are 15–35% below the data for Δr < 1; CR-BLC mode 2 is closest (MC/data 0.8–0.95).

**Status 2026-10-27**: correlations done for CUETP8M1, CR-BLC 2, CR-BLC 3, CP5 (40 jobs, cluster 9726945, all OK);
CR-BLC 0, CP5-CR1, CP5-CR2 submitted (30 jobs, cluster 9744876). Merging + jetshape running for the first four.

Official sample, share of weighted jets > 60 GeV per pTHat window: 15–30: 0.2%, 30–50: 17.6%, 50–80: 52.4%,
80–120: 23.2%, 120+: 6.5%.

---

## 3d. p–Pb with Angantyr (step 4, started 2026-10-27)

All results above are **pp** at √s = 8.16 TeV (boosted), like the official CMS PYTHIA sample. The data are p–Pb, and
the multiplicity classes need p–Pb: pp jet events reach N_trk^offline ≈ 185 essentially never.

**N_trk mapping**: in the official PYTHIA skims N_trk^offline (hiBin) = 0.844 × gen N_ch (charged, |η| < 2.4,
pT > 0.4), flat from 10 to 150 (0.83–0.84), event-by-event rms 6%. Data class edges in gen N_ch: 60 → 71, 80 → 95,
185 → 220, 250 → 298. Used as `ntrkScale = 0.844` in `genV0_analyzer.cpp` (assumes the same tracking efficiency in p–Pb).

**pp jet events per data class** (CUETP8M1, 10M events): < 60: ~70%, 60–80: ~25%, 80–185: ~9%, 185–250: ~0 (32
effective events), ≥ 250: 0 → in pp only the three MB classes are possible.

**Setup** (`cards/ppb.cmnd`, `cards/sigfit/<tune>.cmnd`, `condor/run_job.sh` process `dijetppb`):
- Beams p (6.5 TeV, +z) + Pb (`1000822080`, 2562.5 GeV per nucleon), PYTHIA 8.313 Angantyr (`HeavyIon:mode = 1`).
  The dijet process (same `dijet.cmnd`) is generated in one NN sub-collision, the others are minimum bias /
  diffractive.
- **pTHat bias weight**: with ion beams `info.weight()` contains only the impact-parameter weight, not the signal-process
  bias (checked: weighted and unweighted pTHat spectra almost identical). `runPythia` now multiplies it in:
  w = w_b × (15/pTHat)^4.5 (the exact weight PYTHIA uses in pp). Check: weighted p–Pb pTHat spectrum = pp spectrum.
  The ion check reads `Beams:idA/idB` before `init()` (Angantyr resets them to nucleons).
- Angantyr NN cross-section fit (≈ 100 s at every start-up) done once and fixed in `cards/sigfit/`:
  `SigFitNGen = 0`, `SigFitDefAvNDb = 0.97`, `SigFitDefPar = 1.80,15.72,0.28` (identical for CUETP8M1 and CR-BLC
  mode 2: same Monash total cross sections). CP5-family tunes would need their own fit (different SigmaTotal).
- Gen particles stored only within ΔR < 0.8 of a stored jet (`PythiaCR:genParticleJetR`, p–Pb only): 2.1 kB/event
  instead of 4.5; correlations identical with and without (checked on 3000 events).
- **Colour reconnection acts only inside each NN sub-collision** (PYTHIA 8.3 default). CR between sub-collisions
  (Lönnblad & Shah, EPJC 83 (2023) 575) has no switch in 8.313. So the jet's own sub-collision hadronises as in pp;
  p–Pb changes the surrounding event (multiplicity, UE).
- Timing: ~35 ms/event (CR-BLC 2, lxplus), 2000-event test 89 s incl. 22 s start-up.

**p–Pb jet events per data class** (3000-event CUETP8M1 test, 1980 with a jet > 60 GeV; ⟨gen N_ch⟩ = 111 vs 62 in pp):
< 60: 13%, 60–80: 24%, 80–185: 63%, 185–250: 0.2–0.5%, ≥ 250: 0 → for 10M events ≈ 1.4M / 1.6M / 3.5M / ~50k / 0 jet
events. HM250 needs a central-enriched sample (smaller `HeavyIon:bWidth`, the impact-parameter weight keeps it unbiased).

**Correlations**: `genV0_analyzer.cpp(..., multBinned = true, ntrkScale = 0.844)` fills the correlation sparses in the
data N_trk binning (`CorrCfg::MultBins()`); `condor/corr_job.sh` sets this automatically for `dijetppb_*`;
`postprocessing/run_tune.sh` also runs `jetshape.cpp` for the classes 0–60, 60–80, 80–185, 185–250, 250+.
Data per class (nominal): `postprocessing/run_data_classes.sh` (as `run_jetshapes_data.C`), inputs symlinked in
`postprocessing/input_data_nominal/`.

**Production**: `condor_submit submit.sub PROC=dijetppb TUNES="cuetp8m1 crblc_mode2 crblc_mode3" NJOBS=200 NEV=50000 FLAVOUR=workday SEEDOFFSET=6000`
(cluster 9747155, seeds 6001–6600, 10M events per tune, ~21 GB per tune on user EOS).

**Central-enriched sample for HM185/HM250** (`cards/central.cmnd`: `HeavyIon:bWidth = 1.5` fm; process `dijetppbcentral`).
Impact-parameter width tests, CUETP8M1, 2000 events each (jet events per class < 60 / 60–80 / 80–185 / HM185 / HM250):

| bWidth | jet events | counts | max gen N_ch |
|---|---|---|---|
| default (3000 events) | 1980 | 427 / 477 / 1061 / 15 / 0 | 279 |
| 4.0 fm | 1352 | 217 / 288 / 829 / 18 / 0 | 262 |
| 2.5 fm | 1339 | 139 / 251 / 916 / 33 / 0 | 293 |
| 1.5 fm | 1366 | 108 / 231 / 981 / 44 / 2 | 316 |

1.5 fm: ~4× more HM185 jet events per event than the default and ~1e-3 HM250 events (≈ 10k HM250 jet events per 10M);
narrower widths gain little (the multiplicity tail comes from fluctuations, not geometry). The events keep the
impact-parameter weight, so the sample is unbiased within itself; it is used for the HM classes (and as a cross-check of
the MB classes), not merged with the default sample.
`condor_submit submit.sub PROC=dijetppbcentral TUNES="cuetp8m1 crblc_mode2" NJOBS=200 NEV=50000 FLAVOUR=workday SEEDOFFSET=7000`
(cluster 9747318, seeds 7001–7400, ~42 GB; user EOS had ~270 GB free before the p–Pb productions).

**Data per class** (nominal, `run_data_classes.sh`, in-cone Λ/2K0S, Δr < 0.4, stat. only):

| class (N_trk^offline) | 3–10 GeV | 1.5–2 | 2–3 | 3–5 | 5–7 | 7–10 |
|---|---|---|---|---|---|---|
| MB integrated | 0.2412 ± 0.0023 | 0.154 | 0.260 | 0.247 | 0.256 | 0.220 |
| MB 0–60 | 0.2138 ± 0.0056 | 0.159 | 0.227 | 0.204 | 0.237 | 0.205 |
| MB 60–80 | 0.2265 ± 0.0054 | 0.154 | 0.256 | 0.239 | 0.231 | 0.205 |
| MB 80–185 | 0.2494 ± 0.0029 | 0.151 | 0.267 | 0.255 | 0.264 | 0.227 |
| HM185 (185–250) | 0.2759 ± 0.0026 | 0.120 | 0.264 | 0.291 | 0.295 | 0.234 |
| HM250 (≥ 250) | 0.2893 ± 0.0081 | 0.107 | 0.293 | 0.297 | 0.314 | 0.254 |

---

## 4. Tests and numbers

| Date | Test | Result |
|---|---|---|
| 2026-10-21 | First local run (min. bias, 2k events) | 13 s; works. "diquark too close to junction" errors in CR mode 1 are normal (PYTHIA retries). |
| 2026-10-21 | Condor test, 4 jobs (min. bias, Monash + CR-BLC 0/2/3) | All OK. CPU 27 s (Monash) to 50–130 s (CR-BLC) per 5k events. Only ~0.4% of min.-bias events have a 10 GeV charged jet → jets need the dijet sample. |
| 2026-10-21 | Dijet CUETP8M1, 2k events | 8 ms and 3.5 kB per event; boost OK (⟨η⟩ < 0); Σ0 absent from the record, 36.5% of Λ from Σ0. |
| 2026-10-21 | pTHatMin 15 vs 30 (20k events each) | Jets > 60 GeV per event only 11% fewer with 15 (27.5k vs 30.9k jets / 20k events). No weighted jets from pTHat < 30 at this statistics; 30–40 gives ~5%. → pTHatMin = 15 (like the official samples). |
| 2026-10-21 | Condor test to EOS, dijet, 4 tunes × 2k | All OK; output and job logs copied to EOS. |

---

## 5. Open points / information needed

- ~~Tune of the official sample~~: confirmed CUETP8M1 (see §2). How it was found:
  `voms-proxy-init -voms cms -rfc -valid 192:00`;
  `dasgoclient -query="config dataset=/Dijet_pThat-15_pPb-Bst_8p16_Pythia8/pPb816Summer16DR-pPbBst_80X_mcRun2_pA_v4-v1/AODSIM"`
  → task `pdmvserv_task_HIN-pPb816Spring16GS-00052_...`; fragment from
  `https://cms-pdmv-prod.web.cern.ch/mcm/public/restapi/requests/get_fragment/HIN-pPb816Spring16GS-00052`.
- `ParticleDecays:allowPhotonRadiation = on` (CMS common settings in 80X) does not exist in PYTHIA 8.3 (the run stops
  with "input string not found"), so it cannot be set. It only concerned QED radiation in hadron decays; the
  validation (§3) shows no effect on the observables used here.
- The official sample was made with PYTHIA 8.2 (CMSSW 8_0_X); some defaults changed in 8.3. If small differences
  appear in the validation this is the first suspect.
- Gen N_trk counts all charged final-state particles (including Σ±, Ξ±, Ω±, which are stable in the record);
  the data N_trk counts primary reconstructed tracks. Only relevant for multiplicity classes (pPb step).

---

## 6. Journal

**2026-10-21**
- Created `/afs/cern.ch/work/e/eprifti/PythiaCR` (LCG_107a view, Makefile, cards, condor). Checked condor access, no admin rights needed.
- First version filled ALICE-like histograms directly; replaced by tree output in skim layout (decision: trees + readers).
- Verified CR-BLC parameters against the PYTHIA manual and arXiv:1505.01681 (table above).
- Switched beams to the boosted lab frame, CMS decay settings, dijet with pTHat bias from 15 GeV, CUETP8M1 card.
- Fixed `condor/submit.sub` so that `NAME=value` on the command line overrides the defaults (`$(NAME:default)` syntax).
- Output moved from the group EOS area (96.8% full) to user EOS.
- Submitted CUETP8M1 production: `condor_submit submit.sub TUNES=cuetp8m1 NJOBS=100 NEV=100000 FLAVOUR=longlunch SEEDOFFSET=1000` (cluster 9724509, seeds 1001–1100).
- Filled the official-sample comparison histograms (`compare/official_pgoing.root`).
- Tune of the official sample confirmed as CUETP8M1 from McM. `allowPhotonRadiation` from the CMS common settings
  tried and removed (not a PYTHIA 8.3 setting).
- CUETP8M1 production finished (100/100 jobs exit 0, 30 CPU h, 31 GB). Validation against the official sample passed (§3).
- Found the beam direction mirrored w.r.t. the official pgoing sample; swapped the beams in `common.cmnd`, renamed the
  first production to `dijet_cuetp8m1_ycmNeg`.
- Submitted the main production with identical cards for all tunes:
  `condor_submit submit.sub NJOBS=100 NEV=100000 FLAVOUR=workday SEEDOFFSET=2000`
  (cluster 9724871: cuetp8m1, crblc_mode0, crblc_mode2, crblc_mode3, 10M events each, ~125 GB on user EOS).

**2026-10-22**
- Production 9724871 finished: 397/400 files OK (all > 328 MB). Mean wall time per 100k-event job: CUETP8M1 0.3 h,
  CR-BLC mode 0 2.6 h (max 3.3 h), mode 2 0.5 h, mode 3 0.7 h. Three jobs (mode 0 seeds 2142, 2160; mode 3 seed 2384)
  were removed by condor after 8 h (workday limit) — far above the normal time, so stuck in some event rather than
  slow. Same seed would hang again → replaced with new seeds: clusters 9726939 (mode 0, seeds 3001–3002) and 9726940
  (mode 3, seed 3003).
- Added `hPthat` (all events) and `hPthatJet` (events with a jet > 60) to `genCompare.C`, with the official pthat
  thresholds as exact bin edges; `plotCompare.C` prints σ(pthat > X)/σ(pthat > 15) next to `CorrCfg::SigmaRel`.
- Histogramming moved to condor (`condor/compare.sub`, `compare_job.sh`, list `compare_list.txt`: 11 official files +
  397 generated files), cluster 9726941; lxplus was too loaded (load ~35).
- Added the CMS CP5, CP5-CR1, CP5-CR2 tunes (cards copied from CMSSW, table in §2); 300-event local tests OK.
  Submitted: `condor_submit submit.sub TUNES="cp5 cp5cr1 cp5cr2" NJOBS=100 NEV=100000 FLAVOUR=workday SEEDOFFSET=4000`
  (cluster 9726942). User EOS before this: 529 GB of 1 TB (logical).

- Histogramming jobs 408/408 OK. pTHat check and first tune comparison: §3 and §3b.

- Copied the analysis correlation/post-processing code and wrote the gen-only analyzer (§3c). References processed:
  the copied chain gives 0.1964 for the official PYTHIA (manuscript 0.196); data 0.241–0.259 depending on the version (§3c).
- Correlation jobs for CUETP8M1, CR-BLC 2/3, CP5 (cluster 9726945).

**2026-10-27**
- All productions complete (100 files each for CUETP8M1, CR-BLC 0/2/3, CP5, CP5-CR1, CP5-CR2; the 3 hung CR-BLC jobs replaced).
- Correlation jobs for CR-BLC 0, CP5-CR1, CP5-CR2 submitted (cluster 9744876); merging + jetshape of the first four samples.
- Nominal data fixed (§3c). Plots: variable bins and zoomed ratio panels (`plotCompare.C`, `plot_tunes_vs_data.C`).
- Project put under git and pushed to https://github.com/Enea77/PythiaCR (private; `.gitignore` keeps events, ROOT
  outputs, logs and the EOS symlinks out).
- Residual standalone/official differences traced (§3): official σ_rel(120) ~3.5% high, one pTHat outlier event in the
  official pTHat-30 sample; PYTHIA 8.212 build + production submitted to test the version (cluster 9747148, seeds 5001–5100).

**Seeds used so far**: 1–4 (tests), 1001–1100 (CUETP8M1 ycmNeg), 2001–2400 (cluster 9724871, CUETP8M1 + CR-BLC 0/2/3),
3001–3003 (replacements), 4001–4300 (cluster 9726942, CP5 / CP5-CR1 / CP5-CR2).
5001–5100 (cluster 9747148, CUETP8M1 with PYTHIA 8.212), 6001–6600 (cluster 9747155, p–Pb CUETP8M1 / CR-BLC 2 / CR-BLC 3), 7001–7400 (cluster 9747318, central p–Pb CUETP8M1 / CR-BLC 2). Next production must start at SEEDOFFSET ≥ 8000.
