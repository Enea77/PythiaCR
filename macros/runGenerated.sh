#!/bin/bash
# Gen-level histograms (macros/genCompare.C) of all runPythia outputs of one sample on EOS, then hadd.
#   cd /afs/cern.ch/work/e/eprifti/PythiaCR && macros/runGenerated.sh dijet_cuetp8m1
#   -> compare/dijet_cuetp8m1.root
set -e
S=$1
EOSDIR=/eos/user/e/eprifti/PythiaCR/$S
XRD=root://eosuser.cern.ch/$EOSDIR
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc13-opt/setup.sh
mkdir -p compare/$S
root -l -b -q -e '.L macros/genCompare.C+' > /dev/null
i=0
for f in $(ls $EOSDIR | grep '\.root$'); do
  root -l -b -q "macros/genCompare.C+(\"$XRD/$f\",\"compare/$S/$f\")" > compare/$S/${f%.root}.log 2>&1 &
  i=$((i + 1))
  (( i % 8 == 0 )) && wait
done
wait
echo "$i files; $(grep -L 'events' compare/$S/*.log | wc -l) failed"
hadd -f compare/$S.root compare/$S/*.root > /dev/null
echo "wrote compare/$S.root"
