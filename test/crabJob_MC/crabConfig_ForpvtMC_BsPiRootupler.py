from CRABClient.UserUtilities import config
import datetime
import time

config = config()

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d-%H-%M')

channel = 'BcToBsPi'
year = '2018'
step = 'FlatNtuplesForpvtMC'+year
myname = step+'-'+channel

config.section_('General')
config.General.requestName = step+'-'+channel+'-'+st
config.General.workArea = 'crab_'+step+'-'+channel+'-'+st
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../BspiRootupler.py'
config.JobType.allowUndistributedCMSSW =True
#config.JobType.maxMemoryMB = 3000

config.section_('Data')
config.Data.inputDBS = 'phys03'
#config.Data.inputDataset = '/PrivateMC-2018-BcToBsPi/moameen-BcToBsPi_13TeV_bcvegpy_pythia8_evtgen-381607a57a83c31d7571aa76fa4a79fc/USER'
config.Data.userInputFiles = open('input.txt').readlines()
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.allowNonValidInputDataset = True
#config.Data.publication = True
#config.Data.outLFNDirBase = '/eos/cms/store/group/phys_bphys/moameen/MCSamples/'  #This path will not work.
config.Data.outputPrimaryDataset = myname
config.Data.outLFNDirBase = '/store/user/moameen/FlatNtuplesForBsPi-28July2022'

config.section_('Site')
#config.Site.storageSite = 'T2_CH_CERN'
#config.Data.outLFNDirBase = '/store/user/moameen/'
config.Site.storageSite = 'T3_CH_CERNBOX'
