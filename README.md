# Intructions to run the analyzer

source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_2_9

cd CMSSW_10_2_9/src/

cmsenv

voms-proxy-init -voms cms -valid 192:00

mkdir -p myAnalyzers

cd myAnalyzers

git clone git@github.com:mobassirameen/BspiPAT.git

scram b -j12

cd BspiPAT/test/

cmsRun BspiRootupler.py

