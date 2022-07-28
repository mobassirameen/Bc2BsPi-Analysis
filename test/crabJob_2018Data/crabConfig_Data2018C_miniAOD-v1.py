from CRABClient.UserUtilities import config
import datetime
import time

config = config()

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d-%H-%M')

channel = 'BcToBsPi'
year = '2018C'
step = 'Data'+year
myname = step+'-'+channel

config.General.requestName = step+'-'+channel+'-'+st
config.General.workArea = 'crab_'+step+'-'+channel+'-'+st
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../BspiRootupler_2018data.py'

config.Data.inputDataset = '/Charmonium/Run2018C-17Sep2018-v1/MINIAOD'

config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 50

config.JobType.allowUndistributedCMSSW = True

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_MuonPhys.txt'

#config.Data.runRange = '315252-316995' #Era A
#config.Data.runRange = '317080-319310' #Era B
config.Data.runRange = '319337-320065' #era C
#config.Data.runRange = '320673-325175' #EraD
                  

config.Data.outLFNDirBase = '/store/user/moameen/Bc2BsPi_Data2018C_08Jun2022'
config.Data.publication = False

config.Site.storageSite = 'T3_CH_CERNBOX'
