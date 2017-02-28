#!/bin/bash
echo "hadding dataset: DoubleEG"
eval ./haddDataCrabOutputs2016.sh /eos/uscms/store/user/ztao/ttH_80X/DoubleEG/ [b,c,d,e,f,g,h_v2,h_v3] 24feb2017 170224

echo "hadding dataset: DoubleMuon"
eval ./haddDataCrabOutputs2016.sh /eos/uscms/store/user/ztao/ttH_80X/DoubleMuon/ [b,c,d,e,f,g,h_v2,h_v3] 24feb2017 170224

echo "hadding dataset: MuonEG"
eval ./haddDataCrabOutputs2016.sh /eos/uscms/store/user/ztao/ttH_80X/MuonEG/ [b,c,d,e,f,g,h_v2,h_v3] 24feb2017 170224

echo "hadding dataset: SingleElectron"
eval ./haddDataCrabOutputs2016.sh /eos/uscms/store/user/ztao/ttH_80X/SingleElectron/ [b,c,d,e,f,g,h_v2,h_v3] 24feb2017 170224

echo "hadding dataset: SingleMuon"
eval ./haddDataCrabOutputs2016.sh /eos/uscms/store/user/ztao/ttH_80X/SingleMuon/ [b,c,d,e,f,g,h_v2,h_v3] 24feb2017 170224
