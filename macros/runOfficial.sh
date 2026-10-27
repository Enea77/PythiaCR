#!/bin/bash
# Gen-level histograms of the official pgoing PYTHIA skims, each pthat sample filling its window
# (pthatGen, next threshold] with weight SigmaRel(pthatGen) / N_events (as CorrCfg::SigmaRel), then hadd.
#   cd /afs/cern.ch/work/e/eprifti/PythiaCR && macros/runOfficial.sh
set -e
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc13-opt/setup.sh
D=root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/eprifti/JetV0_skims/MC/unembedded
THR=(15 30 50 80 120 170 220 280 370 460 540)
SIG=(1.0000000 7.4232802e-2 9.5408304e-3 1.2977183e-3 2.1848571e-4 4.0848861e-5 1.1474575e-5 3.3198812e-6 7.3406939e-7 2.0944659e-7 7.9191586e-8)
mkdir -p compare/official
root -l -b -q -e '.L macros/genCompare.C+' > /dev/null
for i in "${!THR[@]}"; do
  lo=${THR[$i]}; hi=${THR[$((i+1))]:-100000}
  root -l -b -q "macros/genCompare.C+(\"$D/V0JetTrkSkims_pgoing_pthat${lo}_wXi.root\",\"compare/official/pthat${lo}.root\",-${SIG[$i]},$lo,$hi)" \
    > compare/official/pthat${lo}.log 2>&1 &
  (( (i + 1) % 6 == 0 )) && wait
done
wait
grep -h -E 'events|Error' compare/official/pthat*.log
hadd -f compare/official_pgoing.root compare/official/pthat*.root > /dev/null
echo "wrote compare/official_pgoing.root"
