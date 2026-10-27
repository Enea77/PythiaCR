#!/bin/bash
# Lines "<sample> <10 runPythia files joined by +> <prefix>" for corr.sub (condor splits queue items on commas).
for s in "$@"; do
  ls /eos/user/e/eprifti/PythiaCR/$s | grep '\.root$' | sort -V | \
    awk -v s=$s -v u=root://eosuser.cern.ch//eos/user/e/eprifti/PythiaCR/$s '{f[NR]=u"/"$0} END {for (i=1;i<=NR;i+=10) {l=f[i]; for (j=i+1;j<i+10&&j<=NR;j++) l=l"+"f[j]; printf "%s %s %s_corr%02d\n", s, l, s, int((i-1)/10)}}'
done
