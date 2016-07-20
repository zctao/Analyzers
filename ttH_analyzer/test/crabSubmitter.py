# submit multiple crab jobs
# modified from Jorge's script

import os

string = '''from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '%(name)s'
config.General.transferLogs = True
config.General.workArea = 'crab'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ttH_2lss_tau_cfg.py'
config.JobType.pyCfgParams = %(cfgparams)s

config.Data.inputDataset = %(dataset)s
config.Data.splitting = '%(splitting)s'
config.Data.unitsPerJob = %(unit)s
%(lumimask)s
%(runrange)s
config.Data.ignoreLocality = True
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/ztao/ttH_76X'

config.Site.storageSite = 'T3_US_FNALLPC'
'''

channels = ['ttH_htt', 'ttH_hww', 'ttH_hzz', 'ttW', 'ttZ', 'rares_TTTT', 'flips_data_dimu_2015c', 'fakes_data_mu_2015d', 'data_obs_dimu_2015d']

for ch in channels:
    with open("SampleList.txt") as f:
        for line in f:
            if ch in line:
                sample = f.next().strip()
                if 'data' in ch:
                    run = f.next().strip()
                pset = f.next().strip()
                perjob = f.next().strip()
                
    vd = locals()
    vd['name'] = ch
    vd['dataset'] = sample
    vd['cfgparams'] = pset
    vd['unit'] = perjob
    if 'data' in ch:
        vd['lumimask'] = "config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'"
        vd['runrange'] = 'config.Data.runRange = '+ run
        vd['splitting'] = 'LumiBased'
    else:
        vd['lumimask'] = ''
        vd['runrange'] = ''
        vd['splitting'] = 'EventAwareLumiBased'

    open('crabConfig_'+ch+'.py','wt').write(string % vd)
    os.system('crab submit -c crabConfig_'+ch+'.py')
