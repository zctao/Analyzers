# submit multiple crab jobs
# modified from Jorge's script

import os

string = '''from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '%(name)s'
config.General.transferLogs = True
config.General.workArea = '/uscms/home/ztao/nobackup/crab'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'CU_ttH_EDA_cfg.py'
config.JobType.pyCfgParams = %(cfgparams)s
config.JobType.sendExternalFolder = True

config.Data.inputDataset = %(dataset)s
config.Data.splitting = '%(splitting)s'
config.Data.unitsPerJob = %(unit)s
%(lumimask)s
%(runrange)s
config.Data.ignoreLocality = True
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/ztao/ttH_80X'

config.Site.storageSite = 'T3_US_FNALLPC'
'''

channels = [#'ttH', 'ttH_jesup', 'ttH_jesdown',
            #'TTW', 'TTW_ext',
            #'TTW_jesup', 'TTW_ext_jesup', 'TTW_jesdown', 'TTW_ext_jesdown',
            #'TTZ', 'TTZ_jesup', 'TTZ_jesdown',
            #'WZ', 'WZ_jesup', 'WZ_jesdown',
            #'ZZ', 'ZZ_jesup', 'ZZ_jesdown',
            #'WW', 'WW_jesup', 'WW_jesdown',
            #'WJets', 'WJets_jesup', 'WJets_jesdown',
            #'WG', 'WG_jesup', 'WG_jesdown',
            #'ZG', 'ZG_jesup', 'ZG_jesdown',
            #'DYJets_M10to50', 'DYJets_M10to50_jesup', 'DYJets_M10to50_jesdown',
            #'DYJets_M50', 'DYJets_M50_jesup', 'DYJets_M50_jesdown',
            #'WWds', 'WWds_jesup', 'WWds_jesdown',
            #'WpWp', 'WpWp_jesup', 'WpWp_jesdown',
            #'WZZ', 'WZZ_jesup', 'WZZ_jesdown',
            #'WWZ', 'WWZ_jesup', 'WWZ_jesdown',
            #'WWW', 'WWW_jesup', 'WWW_jesdown',
            #'ZZZ', 'ZZZ_jesup', 'ZZZ_jesdown',
            #'TTTT', 'TTTT_jesup', 'TTTT_jesdown',
            #'TTJets_DiLep', 'TTJets_DiLep_jesup', 'TTJets_DiLep_jesdown',
            #'TTJets_DiLep_ext', 'TTJets_DiLep_ext_jesup', 'TTJets_DiLep_ext_jesdown',
            #'TTJets_LepT', 'TTJets_LepT_jesup', 'TTJets_LepT_jesdown',
            #'TTJets_LepT_ext', 'TTJets_LepT_ext_jesup', 'TTJets_LepT_ext_jesdown',
            #'TTJets_LepTbar', 'TTJets_LepTbar_jesup', 'TTJets_LepTbar_jesdown',
            #'TTJets_LepTbar_ext', 'TTJets_LepTbar_ext_jesup', 'TTJets_LepTbar_ext_jesdown',
            #'TTGJets', 'TTGJets_jesup', 'TTGJets_jesdown',
            #'TTGJets_ext', 'TTGJets_ext_jesup', 'TTGJets_ext_jesdown',
            #'ST_sLep', 'ST_sLep_jesup', 'ST_sLep_jesdown',
            #'ST_tT', 'ST_tT_jesup', 'ST_tT_jesdown',
            #'ST_tTbar', 'ST_tTbar_jesup', 'ST_tTbar_jesdown',
            #'ST_tWT', 'ST_tWT_jesup', 'ST_tWT_jesdown',
            #'ST_tWTbar', 'ST_tWTbar_jesup', 'ST_tWTbar_jesdown',
            #'TGJets', 'TGJets_jesup', 'TGJets_jesdown',
            #'tZq', 'tZq_jesup', 'tZq_jesdown',
            #'flips_data_dimu_2016b','flips_data_dimu_2016c','flips_data_dimu_2016d',
            #'flips_data_dimu_2016e','flips_data_dimu_2016f','flips_data_dimu_2016g',
            #'flips_data_mu_2016b','flips_data_mu_2016c','flips_data_mu_2016d',
            #'flips_data_mu_2016e','flips_data_mu_2016f','flips_data_mu_2016g',
            #'flips_data_mueg_2016b','flips_data_mueg_2016c','flips_data_mueg_2016d',
            #'flips_data_mueg_2016e','flips_data_mueg_2016f','flips_data_mueg_2016g',
            #'flips_data_e_2016b','flips_data_e_2016c','flips_data_e_2016d',
            #'flips_data_e_2016e','flips_data_e_2016f','flips_data_e_2016g',
            #'flips_data_dieg_2016b','flips_data_dieg_2016c','flips_data_dieg_2016d',
            #'flips_data_dieg_2016e','flips_data_dieg_2016f','flips_data_dieg_2016g',
            #'fakes_data_dimu_2016b','fakes_data_dimu_2016c','fakes_data_dimu_2016d',
            #'fakes_data_dimu_2016e','fakes_data_dimu_2016f','fakes_data_dimu_2016g',
            #'fakes_data_mu_2016b','fakes_data_mu_2016c','fakes_data_mu_2016d',
            #'fakes_data_mu_2016e','fakes_data_mu_2016f','fakes_data_mu_2016g',
            #'fakes_data_mueg_2016b','fakes_data_mueg_2016c','fakes_data_mueg_2016d',
            #'fakes_data_mueg_2016e','fakes_data_mueg_2016f','fakes_data_mueg_2016g',
            #'fakes_data_e_2016b','fakes_data_e_2016c','fakes_data_e_2016d',
            #'fakes_data_e_2016e','fakes_data_e_2016f','fakes_data_e_2016g',
            #'fakes_data_dieg_2016b','fakes_data_dieg_2016c','fakes_data_dieg_2016d',
            #'fakes_data_dieg_2016e','fakes_data_dieg_2016f','fakes_data_dieg_2016g',
            #'data_obs_dimu_2016b', 'data_obs_dimu_2016c', 'data_obs_dimu_2016d',
            #'data_obs_dimu_2016e', 'data_obs_dimu_2016f', 'data_obs_dimu_2016g',
            #'data_obs_mu_2016b', 'data_obs_mu_2016c', 'data_obs_mu_2016d',
            #'data_obs_mu_2016e', 'data_obs_mu_2016f', 'data_obs_mu_2016g',
            #'data_obs_mueg_2016b', 'data_obs_mueg_2016c', 'data_obs_mueg_2016d',
            #'data_obs_mueg_2016e', 'data_obs_mueg_2016f', 'data_obs_mueg_2016g',
            #'data_obs_e_2016b', 'data_obs_e_2016c', 'data_obs_e_2016d',
            #'data_obs_e_2016e', 'data_obs_e_2016f', 'data_obs_e_2016g',
            #'data_obs_dieg_2016b', 'data_obs_dieg_2016c', 'data_obs_dieg_2016d',
            #'data_obs_dieg_2016e', 'data_obs_dieg_2016f', 'data_obs_dieg_2016g',
            ]

def remove_prefix(text, prefix):
    return text[text.startswith(prefix) and len(prefix):]

for ch in channels:
    with open("../data/SampleList_Moriond17.txt") as f:
        for line in f:
            if not 'data' in ch:
                if line.strip()==ch.strip("_jesup") or line.strip()==ch.strip("_jesdown"):
                    sample = f.next().strip()
                    pset = f.next().strip()
                    if '_jesup' in ch:
                        pset=pset.replace("doSystematics=True","doSystematics=False")
                        pset=pset.replace("JECType=NA","JECType=JESUp")
                    if '_jesdown' in ch:
                        pset=pset.replace("doSystematics=True","doSystematics=False")
                        pset=pset.replace("JECType=NA","JECType=JESDown")

                    perjob = f.next().strip()
                    break
                
            else:
                if line.strip()==remove_prefix(ch,"data_obs") or line.strip()==remove_prefix(ch,"fakes_data") or line.strip()==remove_prefix(ch,"flips_data"):
                    sample = f.next().strip()
                    run = f.next().strip()
                    perjob = f.next().strip()
                    pset="['isData=True','SampleName=','SelectionRegion=','TurnOffHLTCut=True']"
                    if 'fakes_' in ch:
                        pset=pset.replace("SampleName=","SampleName=fakes_data")
                        pset=pset.replace("SelectionRegion=",
                                          "SelectionRegion=control_1lfakeable")
                    elif 'flips_' in ch:
                        pset=pset.replace("SampleName=","SampleName=flips_data")
                        pset=pset.replace("SelectionRegion=",
                                          "SelectionRegion=control_2los1tau")
                    elif 'obs_' in ch:
                        pset=pset.replace("SampleName=","SampleName=data_obs")
                        pset=pset.replace("SelectionRegion=",
                                          "SelectionRegion=signal_2lss1tau")
                    else:
                        print 'WARNING: channel name is not valid!!!'
                        exit()
                    break
                         
    vd = locals()
    vd['name'] = ch
    vd['dataset'] = sample
    vd['cfgparams'] = pset
    vd['unit'] = perjob
    if 'data' in ch:
        vd['lumimask'] = "config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'"
        vd['runrange'] = 'config.Data.runRange = '+ run
        vd['splitting'] = 'LumiBased'
    else:
        vd['lumimask'] = ''
        vd['runrange'] = ''
        vd['splitting'] = 'EventAwareLumiBased'

    open('crab/crabConfig_'+ch+'.py','wt').write(string % vd)
    os.system('crab submit -c crab/crabConfig_'+ch+'.py')
