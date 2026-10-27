#!/bin/bash
# Merge the genV0_analyzer outputs of one sample and run jetshape.cpp (gen x gen, default V0 pT bins) on them
# (p-Pb samples dijetppb* incl. dijetppbcentral*: also per N_trk class, outputs jetshapes_rebinned_mult<lo>to<hi>_gen_gen_MC_<V0>.root).
#   ./run_tune.sh dijet_crblc_mode2    -> input_dijet_crblc_mode2/jetshapes_rebinned_gen_gen_MC_{Lam,K0s}.root
# Merged correlations: /eos/user/e/eprifti/PythiaCR/corr/<sample>/merged/MC_<V0>.root (symlinked as input_<sample>/MC_<V0>.root)
set -e
S=$1
cd /afs/cern.ch/work/e/eprifti/PythiaCR/postprocessing
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc13-opt/setup.sh
C=/eos/user/e/eprifti/PythiaCR/corr/$S
mkdir -p $C/merged input_$S
for v in Lam K0s; do
  n=$(ls $C/${S}_corr*_$v.root | wc -l)
  echo "$S $v: $n parts"
  hadd -f -k $C/merged/MC_$v.root $C/${S}_corr*_$v.root > input_$S/hadd_$v.log 2>&1
  ln -sf $C/merged/MC_$v.root input_$S/MC_$v.root
  root -l -b -q "jetshape.cpp+(\"gen\",\"gen\",\"MC_$v\",\"input_$S/\",true,60,2,-1,-1,false,\"\")" > input_$S/js_MC_$v.log 2>&1
  if [[ $S == dijetppb* ]]; then   # p-Pb: also the data N_trk classes (edges of CorrCfg::MultBins)
    for c in 0:60 60:80 80:185 185:250 250:1000; do
      root -l -b -q "jetshape.cpp+(\"gen\",\"gen\",\"MC_$v\",\"input_$S/\",true,60,2,${c%%:*},${c##*:},false,\"\")" > input_$S/js_MC_${v}_mult${c%%:*}.log 2>&1
    done
  fi
done
ls input_$S/jetshapes_rebinned_*
