#!/bin/bash
# Condor wrapper: run_job.sh <process> <tune> <nEvents> <seed>
# process = dijet | minbias | dijetppb (p-Pb, Angantyr) | dijetppbcentral (p-Pb, central-enriched), tune = cuetp8m1 | monash | crblc_mode0 | crblc_mode2 | crblc_mode3
set -e
PROC=$1; TUNE=$2; NEV=$3; SEED=$4
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc13-opt/setup.sh
D=/afs/cern.ch/work/e/eprifti/PythiaCR
if [[ $PROC == dijetppb* ]]; then   # p-Pb (Angantyr): dijet signal in one NN sub-collision, NN cross-section fit per tune
  CARDS="$D/cards/common.cmnd $D/cards/dijet.cmnd $D/cards/ppb.cmnd"
  if [[ $PROC == dijetppbcentral ]]; then CARDS="$CARDS $D/cards/central.cmnd"; fi   # central-enriched (impact-parameter width)
else
  CARDS="$D/cards/common.cmnd $D/cards/${PROC}.cmnd"
fi
if [[ $TUNE == crblc_* ]]; then CARDS="$CARDS $D/cards/crblc_common.cmnd"; fi
CARDS="$CARDS $D/cards/${TUNE}.cmnd"
if [[ $PROC == dijetppb* ]]; then CARDS="$CARDS $D/cards/sigfit/${TUNE}.cmnd"; fi
# Write to the job's local scratch; condor copies the file to EOS on exit
$D/bin/runPythia ${PROC}_${TUNE}_${SEED}.root $NEV $SEED $CARDS
