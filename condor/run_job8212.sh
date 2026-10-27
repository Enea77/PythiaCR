#!/bin/bash
# As run_job.sh, with PYTHIA 8.212 (bin/runPythia8212, LCG_84 CentOS 7); runs in the el7 container (MY.WantOS = "el7")
set -e
PROC=$1; TUNE=$2; NEV=$3; SEED=$4
D=/afs/cern.ch/work/e/eprifti/PythiaCR
source $D/setup8212.sh
CARDS="$D/cards/common.cmnd $D/cards/${PROC}.cmnd"
if [[ $TUNE == crblc_* ]]; then CARDS="$CARDS $D/cards/crblc_common.cmnd"; fi
CARDS="$CARDS $D/cards/${TUNE}.cmnd"
$D/bin/runPythia8212 ${PROC}_${TUNE}_${SEED}.root $NEV $SEED $CARDS
