#!/bin/bash
# Condor wrapper for macros/genCompare.C: compare_job.sh <in.root> <out.root> <scale> <pthatLo> <pthatHi>
set -e
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc13-opt/setup.sh
cd /afs/cern.ch/work/e/eprifti/PythiaCR
# genCompare_C.so is compiled beforehand on lxplus; ACLiC reuses it (the job writes the output in its scratch dir)
root -l -b -q "macros/genCompare.C+(\"$1\",\"$_CONDOR_SCRATCH_DIR/$2\",$3,$4,$5)"
