import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

#--------------- Data --------
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#-----------------------------
#--------------- MC ----------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#-----------------------------
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data') # for 2018, GT for DATA
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v14', '') # for 2018, GT for MC
# GT for MC
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v21', '')# for 2018
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_mc2017_realistic_v8', '')# for 2017
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_mcRun2_asymptotic_v8', '')# for 2016

## Message Logger and Event range
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1200))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

        #MiniAOD Private MC2018
        'gsiftp://eosuserftp.cern.ch/eos/user/m/moameen/MCsamplesForBc2BsPi/PrivateMC-2018-BcToBsPi/BcToBsPi_13TeV_bcvegpy_pythia8_evtgen/220310_180019/0000/BcToBsPi_miniAOD_1.root',
        #'gsiftp://eosuserftp.cern.ch/eos/user/m/moameen/MCsamplesForBc2BsPi/PrivateMC-2018-BcToBsPi/BcToBsPi_13TeV_bcvegpy_pythia8_evtgen/220310_180019/0000/BcToBsPi_miniAOD_10.root',

        #'/store/data/Run2018A/Charmonium/MINIAOD/17Sep2018-v1/00000/01A4FC1E-4A66-544A-B703-31CE3EFF6A28.root',
	#duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
             
 )
)

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring('HLT_Dimuon25_Jpsi_v*',
                                                                        'HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*',
                                                                        'HLT_DoubleMu4_JpsiTrk_Displaced_v*'                                   
                                                                       ),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.load("myAnalyzers.BspiPAT.slimmedMuonsTriggerMatcher_cfi") 

process.rootuple = cms.EDAnalyzer('Bspi',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          Trak_lowpt = cms.InputTag("lostTracks"),        
                          #GenParticles = cms.InputTag("genParticles"),
                          GenParticles = cms.InputTag("prunedGenParticles"),
                          #prunedGenParticles = cms.InputTag("prunedGenParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          OnlyBest = cms.bool(False),
                          #isMC = cms.bool(False), # False for Data
                          isMC = cms.bool(True), # True for MC 
                          OnlyGen = cms.bool(False),
                          Trkmass           = cms.double(0.493677),
                          #Trkmass           = cms.double(0.13957018),
                          TrkTrkMasscut     = cms.vdouble(0.970,1.070),        
                          BarebMasscut      = cms.vdouble(4.4,6.3),
                          bMasscut          = cms.vdouble(5.0,6.0),        
                          )
process.TFileService = cms.Service("TFileService",
        #fileName = cms.string('PvtMC-BsPi_genlevel.root'),
        #fileName = cms.string('PvtMC-BsPi_FlatNtuple_07May2022.root'),
        #fileName = cms.string('PvtMC-BsPi_FlatNtuple_17May2022_v0.root'),
        #fileName = cms.string('PvtMC-BsPi_FlatNtuple_08Jun2022.root'),
        fileName = cms.string('PvtMC-BsPi_FlatNtuple_27Jun2022.root'),
        #fileName = cms.string('BctoBsPi-Data2018A_miniAOD.root'),
        #fileName = cms.string('BctoBsPi-BcmassPeak_Testing.root'),
)


process.mySequence = cms.Sequence(
                                   process.triggerSelection *
				   process.slimmedMuonsWithTriggerSequence *
                                   process.rootuple
				   )
#process.p = cms.Path(process.mySequence)

#process.p = cms.Path(process.triggerSelection*process.rootuple)
#process.p = cms.Path(process.rootuple)
process.p = cms.Path(process.triggerSelection*process.slimmedMuonsWithTriggerSequence*process.rootuple)
