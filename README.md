# Intructions to run the analyzer

source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_2_9

cd CMSSW_10_2_9/src/

cmsenv

voms-proxy-init -voms cms -valid 192:00

git clone git@github.com:mobassirameen/Bc2BsPi-Analysis.git

scram b -j12

cd myAnalyzers/BsPiPAT/test/

cmsRun BsPiRootupler.py

