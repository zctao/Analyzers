#! /bin/sh
# compare object selection level ntuples
echo 'Generate comparison plots for object selection level ntuples...'
root -b 'syncNtupleComparer.C(\
"", \
"~/nobackup/ttHTT_syncNtuple/80X/syncNtuple.root", "syncTree", \
"/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple_ttH_80X.root", "syncTree", \
"", "", \
"", "", \
"~/public_html/syncObjSel/"
)'
# compare event selection level ntuples
echo 'Generate comparison plots for event selection level ntuples...'
# Signal Region
echo 'Signal Region:'
root -b 'syncNtupleComparer.C(\
"2lSS1tau_SR/",
"~/nobackup/ttHTT_syncNtuple/80X/syncNtuple_event.root", "syncTree_2lSS1tau_SR", \
"/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple_event_ttH_80X.root", "syncTree_2lSS1tau_SR", \
"", "", \
"", "", \
"~/public_html/syncEvtSel/"
)'
# Fake
echo 'Fake Lepton Contral Region:'
root -b 'syncNtupleComparer.C(\
"2lSS1tau_Fake/",
"~/nobackup/ttHTT_syncNtuple/80X/syncNtuple_event.root", "syncTree_2lSS1tau_Fake", \
"/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple_event_ttH_80X.root", "syncTree_2lSS1tau_Fake", \
"", "", \
"", "", \
"~/public_html/syncEvtSel/"
)'
# Charge Flip
echo 'Charge Flip Control Region:'
root -b 'syncNtupleComparer.C(\
"2lSS1tau_Flip/",
"~/nobackup/ttHTT_syncNtuple/80X/syncNtuple_event.root", "syncTree_2lSS1tau_Flip", \
"/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple_event_ttH_80X.root", "syncTree_2lSS1tau_Flip", \
"", "", \
"", "", \
"~/public_html/syncEvtSel/"
)'
