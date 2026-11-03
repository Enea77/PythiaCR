#!/bin/bash
# Condor wrapper for bin/aliceCheck (ALICE JHEP07(2023)136 reproduction): alice_job.sh <proc> <tune> <nEvents> <seed>
# proc = pp13_minbias (inclusive) | pp13_hardqcd (in-jet); output: histograms only
set -e
PROC=$1; TUNE=$2; NEV=$3; SEED=$4
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc13-opt/setup.sh
D=/afs/cern.ch/work/e/eprifti/PythiaCR
CARDS="$D/cards/common.cmnd $D/cards/${PROC}.cmnd"
if [[ $TUNE == crblc_* ]]; then CARDS="$CARDS $D/cards/crblc_common.cmnd"; fi
$D/bin/aliceCheck alice_${PROC}_${TUNE}_${SEED}.root $NEV $SEED $CARDS $D/cards/${TUNE}.cmnd
