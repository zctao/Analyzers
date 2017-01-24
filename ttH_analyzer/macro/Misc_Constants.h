#ifndef Misc_Constants_h
#define Misc_Constants_h

#include "TString.h"

#include <map>
#include <vector>

/////////////////////////////////////////////
// Names of b-tagging systematics
const TString BTagSysts[16] =
	{"LFUp","LFDown","HFUp","HFDown",
	 "HFStats1Up","HFStats1Down","HFStats2Up","HFStats2Down",
	 "LFStats1Up","LFStats1Down","LFStats2Up","LFStats2Down",
	 "cErr1Up","cErr1Down","cErr2Up","cErr2Down"};

/////////////////////////////////////////////
// Names of theoretical uncertainties
const TString ThSysts[4] = {"x1Up", "x1Down", "y1Up", "y1Down"};

/////////////////////////////////////////////
// Higgs decay modes in signal sample
const TString HiggsDecayMode[3] = {"htt", "hww", "hzz"};

/////////////////////////////////////////////
// systematic histogram name
const TString sys_coname = "CMS_ttHl_";

/////////////////////////////////////////////
// data samples
const std::vector<TString> DataSamples =
	{"DoubleMuon","SingleMuon","DoubleEG","SingleElectron","MuonEG"};

/////////////////////////////////////////////
// Map of samples used for certain channel
const std::map<TString,std::vector<TString>> SamplesInChannel =
	{ {"ttH", {"ttH"}},
	  {"TTW", {"TTW"}},
	  {"TTZ", {"TTZ"}},
	  {"EWK", {"WZ", "ZZ", "WW", "WG", "ZG"}},
	  {"Rare", {"TTTT", "tqZ", "WZZ"/*, "WWZ", "WWW", "ZZZ"*/}},
	  {"fakes_data", DataSamples},
	  {"flips_data", DataSamples},
	  {"data_obs", DataSamples}
	};

/////////////////////////////////////////////
// Crab job output directories for different samples
const std::map<TString, TString> dir_map =
	{{"ttH", "/eos/uscms/store/user/ztao/ttH_80X/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/"},
	 {"TTZ", "/eos/uscms/store/user/ztao/ttH_80X/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/"},
	 {"TTW", "/eos/uscms/store/user/ztao/ttH_80X/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/"},
	 {"WZ", "/eos/uscms/store/user/ztao/ttH_80X/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/"},
	 {"ZZ", "/eos/uscms/store/user/ztao/ttH_80X/ZZTo4L_13TeV-amcatnloFXFX-pythia8/"},
	 {"WW", "/eos/uscms/store/user/ztao/ttH_80X/WWTo2L2Nu_13TeV-powheg/"},
	 {"WG", "/eos/uscms/store/user/ztao/ttH_80X/WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/"},
	 {"ZG", "/eos/uscms/store/user/ztao/ttH_80X/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/"},
	 {"TTTT", "/eos/uscms/store/user/ztao/ttH_80X/TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/"},
	 //{"tZq", "/"},
	 {"DoubleMuon", "/eos/uscms/store/user/ztao/ttH_80X/DoubleMuon/"},
	 {"SingleMuon", "/eos/uscms/store/user/ztao/ttH_80X/SingleMuon/"},
	 {"DoubleEG", "/eos/uscms/store/user/ztao/ttH_80X/DoubleEG/"},
	 {"SingleElectron", "/eos/uscms/store/user/ztao/ttH_80X/SingleElectron/"},
	 {"MuonEG", "/eos/uscms/store/user/ztao/ttH_80X/MuonEG/"}
	};

#endif
