# PYTHIA 8.212 + ROOT + FastJet from LCG_84 (CentOS 7): source inside /cvmfs/cms.cern.ch/common/cmssw-el7
source /cvmfs/sft.cern.ch/lcg/views/LCG_84/x86_64-centos7-gcc49-opt/setup.sh
# the view has no pythia8-config: point PYTHIA8DATA at the 8.212 install (Makefile.8212 uses the same path)
export PYTHIA8DATA=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/pythia8/212-90044/x86_64-centos7-gcc49-opt/share/Pythia8/xmldoc
