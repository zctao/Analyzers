from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'data_E2015D'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ttH_2lss_tau_cfg.py'
config.JobType.pyCfgParams = ['doSync=False',
                              'isData=True',
                              'SampleName=data_obs',
                              'doScale=False',
                              'doSystematics=False']

#config.Data.inputDataset = '/DoubleMuon/Run2015C_25ns-16Dec2015-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2015D-16Dec2015-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2015C_25ns-16Dec2015-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleEG/Run2015D-04Dec2015-v2/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2015C_25ns-16Dec2015-v1/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2015D-16Dec2015-v1/MINIAOD'
#config.Data.inputDataset = '/SingleMuon/Run2015C_25ns-16Dec2015-v1/MINIAOD'
#config.Data.inputDataset = '/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD'
#config.Data.inputDataset = '/SingleElectron/Run2015C_25ns-16Dec2015-v1/MINIAOD'
config.Data.inputDataset = '/SingleElectron/Run2015D-04Dec2015-v2/MINIAOD'

config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 30
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'
#config.Data.runRange = '254227-254914'    # Run2015C
config.Data.runRange = '256630-260627'    # Run2015D
config.Data.ignoreLocality = True
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/ztao/ttH_76X'

config.Site.storageSite = 'T3_US_FNALLPC'
