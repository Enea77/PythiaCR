# PythiaCR

Detailed record of decisions and results: `LOGBOOK.md`.

Generator-level PYTHIA 8 (8.313, LCG_107a) pp at 8.16 TeV in the pPb lab frame, to compare colour-reconnection
models with the CMS pPb in-jet (Λ+Λ̄)/2K0S measurement. The generator writes trees in the layout of the
V0+jet skims (gen branches only); histograms are made afterwards by readers.

```
source setup.sh && make            # build bin/runPythia
bin/runPythia out.root 2000 1 cards/common.cmnd cards/dijet.cmnd cards/cuetp8m1.cmnd
root -l -b -q 'macros/quickCheck.C+("out.root")'
cd condor && condor_submit submit.sub    # see the header of submit.sub for production arguments
```

## Cards (read in order: common, process, [crblc_common], tune)
| card | content |
|---|---|
| common | beams 6.5 x 2.5625 TeV (y_cm = +0.465, as the official pgoing sample), CMS decays (limitTau0, tau0Max = 10 mm) |
| dijet / minbias | HardQCD with pTHat^4.5 bias from pTHatMin = 15 (fill with `weight`) / SoftQCD nonDiffractive |
| cuetp8m1 | CMS CUETP8M1 (tune of the official Summer16 sample, confirmed in McM) |
| cp5 / cp5cr1 / cp5cr2 | CMS CP5, CP5-CR1 (QCD-based CR), CP5-CR2 (gluon move), copied from CMSSW; NNPDF3.1 via LHAPDF |
| monash | Monash 2013 |
| crblc_common + crblc_mode{0,2,3} | QCD-based CR beyond leading colour, Christiansen & Skands, JHEP 08 (2015) 003, appendix table (as in ALICE JHEP07(2023)136) |

## Output trees (same names and types as the skims)
`hiEvtAnalyzer/HiTree` (pthat, weight, vz = 0, hiBin = gen N_trk, nMPI), `ak4PFJetAnalyzer/t` (ngen, genpt, geneta,
genphi, WTAgeneta, WTAgenphi; anti-kT 0.4, all final-state particles except neutrinos, pT > 30), `HiGenParticleAna/hi`
(|eta| < 2.4, pT > 0.4), `g{K0s,Lam,Xi,Om}Tree` (|eta| < 2.4, mom1 = hadron mother id or 0).
