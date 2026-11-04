#!/bin/bash
# Condor wrapper for correlations/genV0_analyzer.cpp: corr_job.sh <inputs joined by +> <output prefix>
set -e
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc13-opt/setup.sh
cd /afs/cern.ch/work/e/eprifti/PythiaCR/correlations
# genV0_analyzer_cpp.so is compiled beforehand on lxplus; outputs go to the job scratch dir
FILES=${1//+/,}
# p-Pb samples: correlations in the data N_trk classes, N_trk^offline = 0.844 x gen N_trk (see genV0_analyzer.cpp)
OPT=""
# pp variant <sample>_mult (same events): also filled in the data N_trk classes, for the per-class comparison with data
if [[ $2 == dijetppb* || $2 == *_mult_corr* ]]; then OPT=",-1,true,0.844"; else OPT=",-1,false,1.0"; fi
# variant <sample>_pthat<N>: events with pTHat < N skipped
if [[ $2 =~ _pthat([0-9]+)_corr ]]; then OPT="$OPT,${BASH_REMATCH[1]}"; fi
# pp variant <sample>_jpct (same events): classes with the same fraction of jets as in data, gen N_ch edges from
# condor/multmap_<sample>.txt (macros/jpct_edges.py), mapped onto the data N_trk classes
if [[ $2 =~ ^(.*)_jpct_corr ]]; then OPT=",-1,true,1.0,0,\"$(cat /afs/cern.ch/work/e/eprifti/PythiaCR/condor/multmap_${BASH_REMATCH[1]}.txt)\""; fi
root -l -b -q "genV0_analyzer.cpp+(\"$FILES\",\"$_CONDOR_SCRATCH_DIR/$2\"$OPT)"
