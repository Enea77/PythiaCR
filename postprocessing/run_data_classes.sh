#!/bin/bash
# jetshape.cpp on the nominal data (Lambda: Data_v0shift_effOm, K0S: Data_v0shift; symlinks in input_data_nominal/), in-jet
# efficiency correction, classes as run_jetshapes_data.C: MB integrated, MB N_trk 0-60, 60-80, 80-185, HM185 (185-250 in
# the HM185 sample), HM250 (>= 250 in the HM250 sample). Outputs next to the symlinks.
cd /afs/cern.ch/work/e/eprifti/PythiaCR/postprocessing
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc13-opt/setup.sh
D=input_data_nominal/
for v in Lam K0s; do
  ij="files/InJet_Efficiency_$v.root"
  for c in MB:-1:-1 MB:0:60 MB:60:80 MB:80:185 HM185:185:250 HM250:250:1000; do
    IFS=: read s lo hi <<< "$c"
    root -l -b -q "jetshape.cpp+(\"reco\",\"reco\",\"Data_${v}_$s\",\"$D\",true,60,2,$lo,$hi,false,\"\",\"$ij\")" > $D/js_${v}_${s}_${lo}.log 2>&1
    echo "$v $s $lo-$hi rc=$?"
  done
done
ls $D/jetshapes_rebinned_*
