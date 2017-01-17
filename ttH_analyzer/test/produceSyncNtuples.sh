#!/bin/bash
# sync ntuple for signal region
echo 'start producing sync ntuple for signal region ...'
cmsRun CU_ttH_EDA_cfg.py evtSelOff=False TurnOffHLTCut=True SampleName=sync_event_sr doSystematics=True SelectionRegion=signal_2lss1tau reHLT=False GridMode=False doJERSmearing=False

# sync ntuple for fake extrapolation region
echo 'start producing sync ntuple for fake extrapolation region ...'
cmsRun CU_ttH_EDA_cfg.py evtSelOff=False TurnOffHLTCut=True SampleName=sync_event_fake doSystematics=False SelectionRegion=control_1lfakeable reHLT=False GridMode=False doJERSmearing=False

# sync ntuple for charge flip control region
echo 'start producing sync ntuple for charge flip control region ...'
cmsRun CU_ttH_EDA_cfg.py evtSelOff=False TurnOffHLTCut=True SampleName=sync_event_flip doSystematics=False SelectionRegion=control_2los1tau reHLT=False GridMode=False doJERSmearing=False
