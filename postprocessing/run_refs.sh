#!/bin/bash
# jetshape.cpp (copied unchanged) on the reference inputs: official PYTHIA gen x gen and MB data (in-jet efficiency
# correction as run_jetshapes_data.C), V0 pT bins 1.5,2,3,5,7,10 (default). Outputs next to the (symlinked) inputs.
cd /afs/cern.ch/work/e/eprifti/PythiaCR/postprocessing
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc13-opt/setup.sh
root -l -b -q -e '.L jetshape.cpp+' > /dev/null 2>&1
for v in Lam K0s; do
  root -l -b -q "jetshape.cpp+(\"gen\",\"gen\",\"MC_$v\",\"input_official/\",true,60,2,-1,-1,false,\"\")" > input_official/js_MC_$v.log 2>&1
  root -l -b -q "jetshape.cpp+(\"reco\",\"reco\",\"Data_${v}_MB\",\"input_data_v0shift/\",true,60,2,-1,-1,false,\"\",\"files/InJet_Efficiency_$v.root\")" > input_data_v0shift/js_Data_${v}_MB.log 2>&1
done
root -l -b -q "jetshape.cpp+(\"reco\",\"reco\",\"Data_Lam_MB\",\"input_data_v0shift_effOm/\",true,60,2,-1,-1,false,\"\",\"files/InJet_Efficiency_Lam.root\")" > input_data_v0shift_effOm/js_Data_Lam_MB.log 2>&1
ls -la input_*/jetshapes_rebinned_*
