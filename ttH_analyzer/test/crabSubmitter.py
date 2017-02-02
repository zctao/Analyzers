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

channels = [#'ttH','ttH_jesup','ttH_jesdown','ttH_tesup','ttH_tesdown',
            #'TTW', 'TTW_jesup','TTW_jesdown','TTW_tesup','TTW_tesdown',
            #'TTW_ext','TTW_ext_jesup','TTW_ext_jesdown','TTW_ext_tesup','TTW_ext_tesdown',
            #'TTZ','TTZ_jesup','TTZ_jesdown','TTZ_tesup','TTZ_tesdown',
            #'TTGJets','TTGJets_jesup','TTGJets_jesdown','TTGJets_tesup','TTGJets_tesdown',
            #'TTGJets_ext','TTGJets_ext_jesup','TTGJets_ext_jesdown','TTGJets_ext_tesup','TTGJets_ext_tesdown',
            #'TGJets','TGJets_jesup','TGJets_jesdown','TGJets_tesup','TGJets_tesdown',
            #'WG','WG_jesup','WG_jesdown','WG_tesup','WG_tesdown',
            #'ZG','ZG_jesup','ZG_jesdown','ZG_tesup','ZG_tesdown',
            #'WZ','WZ_jesup','WZ_jesdown','WZ_tesup','WZ_tesdown',
            #'ZZ','ZZ_jesup','ZZ_jesdown','ZZ_tesup','ZZ_tesdown',
            #'WW','WW_jesup','WW_jesdown','WW_tesup','WW_tesdown',
            #'WWds','WWds_jesup','WWds_jesdown','WWds_tesup','WWds_tesdown',
            #'WpWp','WpWp_jesup','WpWp_jesdown','WpWp_tesup','WpWp_tesdown',
            #'WZZ','WZZ_jesup','WZZ_jesdown','WZZ_tesup','WZZ_tesdown',
            #'WWZ','WWZ_jesup','WWZ_jesdown','WWZ_tesup','WWZ_tesdown',
            #'WWW','WWW_jesup','WWW_jesdown','WWW_tesup','WWW_tesdown',
            #'ZZZ','ZZZ_jesup','ZZZ_jesdown','ZZZ_tesup','ZZZ_tesdown',
            #'tZq','tZq_jesup','tZq_jesdown','tZq_tesup','tZq_tesdown',
            #'TTTT','TTTT_jesup','TTTT_jesdown','TTTT_tesup','TTTT_tesdown',
            #
            #'TTJets_DiLep',
            #'TTJets_DiLep_jesup','TTJets_DiLep_jesdown','TTJets_DiLep_tesup','TTJets_DiLep_tesdown',
            #'TTJets_DiLep_ext',
            #'TTJets_DiLep_ext_jesup','TTJets_DiLep_ext_jesdown','TTJets_DiLep_ext_tesup','TTJets_DiLep_ext_tesdown',
            #'TTJets_LepT',
            #'TTJets_LepT_jesup','TTJets_LepT_jesdown','TTJets_LepT_tesup','TTJets_LepT_tesdown',
            #'TTJets_LepT_ext',
            #'TTJets_LepT_ext_jesup','TTJets_LepT_ext_jesdown','TTJets_LepT_ext_tesup','TTJets_LepT_ext_tesdown',
            #'TTJets_LepTbar',
            #'TTJets_LepTbar_jesup','TTJets_LepTbar_jesdown','TTJets_LepTbar_tesup','TTJets_LepTbar_tesdown',
            #'TTJets_LepTbar_ext',
            #'TTJets_LepTbar_ext_jesup','TTJets_LepTbar_ext_jesdown','TTJets_LepTbar_ext_tesup','TTJets_LepTbar_ext_tesdown',
            #'ST_sLep',
            #'ST_sLep_jesup','ST_sLep_jesdown','ST_sLep_tesup', 'ST_sLep_tesdown',
            #'ST_tT',
            #'ST_tT_jesup','ST_tT_jesdown','ST_tT_tesup','ST_tT_tesdown',
            #'ST_tTbar',
            #'ST_tTbar_jesup','ST_tTbar_jesdown','ST_tTbar_tesup','ST_tTbar_tesdown',
            #'ST_tWT',
            #'ST_tWT_jesup','ST_tWT_jesdown','ST_tWT_tesup','ST_tWT_tesdown',
            #'ST_tWTbar',
            #'ST_tWTbar_jesup','ST_tWTbar_jesdown','ST_tWTbar_tesup','ST_tWTbar_tesdown',
            #'WJets',
            #'WJets_jesup','WJets_jesdown','WJets_tesup','WJets_tesdown',
            #'DYJets_M10to50',
            #'DYJets_M10to50_jesup','DYJets_M10to50_jesdown','DYJets_M10to50_tesup','DYJets_M10to50_tesdown',
            #'DYJets_M50',
            #'DYJets_M50_jesup','DYJets_M50_jesdown','DYJets_M50_tesup','DYJets_M50_tesdown',
            #
            #'flips_data_dimu_2016b','flips_data_dimu_2016c','flips_data_dimu_2016d',
            #'flips_data_dimu_2016e','flips_data_dimu_2016f','flips_data_dimu_2016g',
            #'flips_data_dimu_2016h_v2', 'flips_data_dimu_2016h_v3',
            #'flips_data_mu_2016b','flips_data_mu_2016c','flips_data_mu_2016d',
            #'flips_data_mu_2016e','flips_data_mu_2016f','flips_data_mu_2016g',
            #'flips_data_mu_2016h_v2', 'flips_data_mu_2016h_v3',
            #'flips_data_mueg_2016b','flips_data_mueg_2016c','flips_data_mueg_2016d',
            #'flips_data_mueg_2016e','flips_data_mueg_2016f','flips_data_mueg_2016g',
            #'flips_data_mueg_2016h_v2', 'flips_data_mueg_2016h_v3',
            #'flips_data_e_2016b','flips_data_e_2016c','flips_data_e_2016d',
            #'flips_data_e_2016e','flips_data_e_2016f','flips_data_e_2016g',
            #'flips_data_e_2016h_v2', 'flips_data_e_2016h_v3',
            #'flips_data_dieg_2016b','flips_data_dieg_2016c','flips_data_dieg_2016d',
            #'flips_data_dieg_2016e','flips_data_dieg_2016f','flips_data_dieg_2016g',
            #'flips_data_dieg_2016h_v2', 'flips_data_dieg_2016h_v3',
            #'fakes_data_dimu_2016b','fakes_data_dimu_2016c','fakes_data_dimu_2016d',
            #'fakes_data_dimu_2016e','fakes_data_dimu_2016f','fakes_data_dimu_2016g',
            #'fakes_data_dimu_2016h_v2', 'fakes_data_dimu_2016h_v3',
            #'fakes_data_mu_2016b','fakes_data_mu_2016c','fakes_data_mu_2016d',
            #'fakes_data_mu_2016e','fakes_data_mu_2016f','fakes_data_mu_2016g',
            #'fakes_data_mu_2016h_v2', 'fakes_data_mu_2016h_v3',
            #'fakes_data_mueg_2016b','fakes_data_mueg_2016c','fakes_data_mueg_2016d',
            #'fakes_data_mueg_2016e','fakes_data_mueg_2016f','fakes_data_mueg_2016g',
            #'fakes_data_mueg_2016h_v2', 'fakes_data_mueg_2016h_v3',
            #'fakes_data_e_2016b','fakes_data_e_2016c','fakes_data_e_2016d',
            #'fakes_data_e_2016e','fakes_data_e_2016f','fakes_data_e_2016g',
            #'fakes_data_e_2016h_v2', 'fakes_data_e_2016h_v3',
            #'fakes_data_dieg_2016b','fakes_data_dieg_2016c','fakes_data_dieg_2016d',
            #'fakes_data_dieg_2016e','fakes_data_dieg_2016f','fakes_data_dieg_2016g',
            #'fakes_data_dieg_2016h_v2', 'fakes_data_dieg_2016h_v3',
            #'data_obs_dimu_2016b', 'data_obs_dimu_2016c', 'data_obs_dimu_2016d',
            #'data_obs_dimu_2016e', 'data_obs_dimu_2016f', 'data_obs_dimu_2016g',
            #'data_obs_dimu_2016h_v2', 'data_obs_dimu_2016h_v3',
            #'data_obs_mu_2016b', 'data_obs_mu_2016c', 'data_obs_mu_2016d',
            #'data_obs_mu_2016e', 'data_obs_mu_2016f', 'data_obs_mu_2016g',
            #'data_obs_mu_2016h_v2', 'data_obs_mu_2016h_v3',
            #'data_obs_mueg_2016b', 'data_obs_mueg_2016c', 'data_obs_mueg_2016d',
            #'data_obs_mueg_2016e', 'data_obs_mueg_2016f', 'data_obs_mueg_2016g',
            #'data_obs_mueg_2016h_v2', 'data_obs_mueg_2016h_v3',
            #'data_obs_e_2016b', 'data_obs_e_2016c', 'data_obs_e_2016d',
            #'data_obs_e_2016e', 'data_obs_e_2016f', 'data_obs_e_2016g',
            #'data_obs_e_2016h_v2', 'data_obs_e_2016h_v3',
            #'data_obs_dieg_2016b', 'data_obs_dieg_2016c', 'data_obs_dieg_2016d',
            #'data_obs_dieg_2016e', 'data_obs_dieg_2016f', 'data_obs_dieg_2016g',
            #'data_obs_dieg_2016h_v2', 'data_obs_dieg_2016h_v3',
            ]

def remove_prefix(text, prefix):
    return text[text.startswith(prefix) and len(prefix):]

sample=""
pset=""

for ch in channels:
    
    del sample
    del pset

    print ch
    
    with open("../data/SampleList_Moriond17.txt") as f:
        for line in f:
            if not 'data' in ch:
                if line.strip()==ch.replace("_jesup","") or line.strip()==ch.replace("_jesdown","") or line.strip()==ch.replace("_tesup","") or line.strip()==ch.replace("_tesdown",""):
                    sample = f.next().strip()
                    pset = f.next().strip()
                    if '_jesup' in ch:
                        pset=pset.replace("doSystematics=True","doSystematics=False")
                        pset=pset.replace("JECType=NA","JECType=JESUp")
                    if '_jesdown' in ch:
                        pset=pset.replace("doSystematics=True","doSystematics=False")
                        pset=pset.replace("JECType=NA","JECType=JESDown")

                    if '_tesup' in ch:
                        pset=pset.replace("doSystematics=True","doSystematics=False")
                        pset=pset.replace("TauESType=NA","TauESType=tauESUp")

                    if '_tesdown' in ch:
                        pset=pset.replace("doSystematics=True","doSystematics=False")
                        pset=pset.replace("TauESType=NA","TauESType=tauESDown")
                        
                    perjob = f.next().strip()
                    break
                
            else:
                if line.strip()==remove_prefix(ch,"data_obs") or line.strip()==remove_prefix(ch,"fakes_data") or line.strip()==remove_prefix(ch,"flips_data"):
                    sample = f.next().strip()
                    run = f.next().strip()
                    perjob = f.next().strip()
                    pset="['isData=True','SampleName=','SelectionRegion=','TurnOffHLTCut=True', 'HIPSafeMediumMuon=True', 'Is2016H=False']"
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

                    if '_2016g' in ch or '_2016h' in ch:
                        pset=pset.replace("HIPSafeMediumMuon=True","HIPSafeMediumMuon=False")
                    if '_2016h' in ch:
                        pset=pset.replace("Is2016H=False","Is2016H=True")
                        
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
