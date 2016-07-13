from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ttH_htt'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ttH_2lss_tau_cfg.py'
config.JobType.pyCfgParams = ['doSync=False',
                              'isData=False',
                              'SampleName=ttH_htt',
                              'IntLumi=2.26',   #fb^-1
                              'CrossSection=258.6', #fb
                              'sysJECType=NA',  # NA, JESUp, JESDown
                              'doSystematics=True'] # True only for NA
#config.JobType.inputFiles

config.Data.inputDataset = '/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 25000
#config.Data.lumiMask
#config.Data.runRange
config.Data.ignoreLocality = True
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/ztao/ttH_76X'

config.Site.storageSite = 'T3_US_FNALLPC'
