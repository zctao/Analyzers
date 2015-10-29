from CRABClient.UserUtilities import config
config = config()

#config.General.requestName = 'ttHTTNtuple_signal'
#config.General.requestName = 'ttHTTNtuple_TTJets'
#config.General.requestName = 'ttHTTNtuple_TTWToLNu'
config.General.requestName = 'ttHTTNtuple_TTWToQQ'

config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/uscms/home/ztao/work/CU_ttH_WD/CMSSW_7_4_7/src/Analyzers/ttH_analyzer/python/CU_ttH_EDA_cfg.py'

#config.Data.inputDataset = '/ttHToTT_M125_13TeV_powheg_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
#config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v3/MINIAODSIM'
#config.Data.inputDataset = '/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
config.Data.inputDataset = '/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.ignoreLocality = True
config.Data.publication = False

#config.Data.publishDataName = 'ttHToTauTau_Ntuple_signal'
#config.Data.publishDataName = 'ttHToTauTau_Ntuple_TTJets'
#config.Data.publishDataName = 'ttHToTauTau_Ntuple_TTWToLNu'
config.Data.publishDataName = 'ttHToTauTau_Ntuple_TTWToQQ'

config.Data.outLFNDirBase = '/store/user/ztao'
config.Site.storageSite = 'T3_US_FNALLPC'
