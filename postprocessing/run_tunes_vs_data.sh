#!/bin/bash
# Data vs PYTHIA tunes with plot_tunes_vs_data.C: MB integrated (pp + p-Pb samples) and per N_trk class (p-Pb samples).
# Samples whose jetshape outputs do not exist yet are skipped.
#   ./run_tunes_vs_data.sh           -> plots/tunes_vs_data_mult/<class>/{InCone_Ratio,DeltaR_Ratio,...,values.txt}
#   ./run_tunes_vs_data.sh _pthat30  -> p-Pb samples <sample>_pthat30 (pTHat > 30), plots/tunes_vs_data_mult_pthat30/
set -e
cd /afs/cern.ch/work/e/eprifti/PythiaCR/postprocessing
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc13-opt/setup.sh
D=input_data_nominal
V=${1:-}   # p-Pb variant suffix
O=plots/tunes_vs_data_mult$V
# "label:sample"; pp samples only enter the integrated MB plot
PP=("pp CUETP8M1:dijet_cuetp8m1" "pp CR-BLC 2:dijet_crblc_mode2")
PPB=("pPb CUETP8M1:dijetppb_cuetp8m1$V" "pPb CR-BLC 2:dijetppb_crblc_mode2$V" "pPb CR-BLC 3:dijetppb_crblc_mode3$V"
     "pPb central CUETP8M1:dijetppbcentral_cuetp8m1$V" "pPb central CR-BLC 2:dijetppbcentral_crblc_mode2$V")

mclist() {  # mclist <jetshape tag ("" or "mult<lo>to<hi>_")> <label:sample>...
  local tag=$1 out="" fl fk
  shift
  for ls in "$@"; do
    fl=input_${ls##*:}/jetshapes_rebinned_${tag}gen_gen_MC_Lam.root
    fk=input_${ls##*:}/jetshapes_rebinned_${tag}gen_gen_MC_K0s.root
    [[ -f $fl && -f $fk ]] && out="$out;${ls%%:*}:$fl:$fk"
  done
  echo "${out#;}"
}

run() {  # run <outdir> <data Lam> <data K0s> <mc list>
  [ -z "$4" ] && { echo "$1: no MC samples yet"; return; }
  root -l -b -q "plot_tunes_vs_data.C(\"$2,$3\",\"$4\",\"$1\")" > /dev/null 2>&1
  echo "== $1"; cat $1/values.txt
}

run $O/MB $D/jetshapes_rebinned_ijc_reco_reco_Data_Lam_MB.root \
    $D/jetshapes_rebinned_ijc_reco_reco_Data_K0s_MB.root "$(mclist "" "${PP[@]}" "${PPB[@]}")"
for c in 0to60:MB 60to80:MB 80to185:MB 185to250:HM185 250to1000:HM250; do
  m=${c%%:*} t=${c##*:}
  run $O/mult$m $D/jetshapes_rebinned_ijc_mult${m}_reco_reco_Data_Lam_$t.root \
      $D/jetshapes_rebinned_ijc_mult${m}_reco_reco_Data_K0s_$t.root "$(mclist mult${m}_ "${PPB[@]}")"
done
