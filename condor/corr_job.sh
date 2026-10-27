#!/bin/bash
# Condor wrapper for correlations/genV0_analyzer.cpp: corr_job.sh <inputs joined by +> <output prefix>
set -e
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc13-opt/setup.sh
cd /afs/cern.ch/work/e/eprifti/PythiaCR/correlations
# genV0_analyzer_cpp.so is compiled beforehand on lxplus; outputs go to the job scratch dir
FILES=${1//+/,}
# p-Pb samples: correlations in the data N_trk classes, N_trk^offline = 0.844 x gen N_trk (see genV0_analyzer.cpp)
OPT=""
if [[ $2 == dijetppb* ]]; then OPT=",-1,true,0.844"; fi
root -l -b -q "genV0_analyzer.cpp+(\"$FILES\",\"$_CONDOR_SCRATCH_DIR/$2\"$OPT)"
