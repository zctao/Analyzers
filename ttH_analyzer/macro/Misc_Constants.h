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
TString batch = "jan2017/";
TString eos_dir = "/eos/uscms/store/user/ztao/ttH_80X/"
const std::map<TString, TString> dir_map =
	{{"ttH", eos_dir+"ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/"+batch},
	 {"TTW", eos_dir+"TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/"+batch},
	 {"TTZ", eos_dir+"TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/"+batch},
	 {"TTGJets", eos_dir+"TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/"+batch},
	 {"TGJets", eos_dir+"TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/"+batch},
	 {"WG", eos_dir+"WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/"+batch},
	 {"ZG", eos_dir+"ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/"+batch},
	 {"WW", eos_dir+"WWTo2L2Nu_13TeV-powheg/"+batch},
	 {"WpWp", eos_dir+"WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/"+batch},
	 {"WZ", eos_dir+"WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/"+batch},
	 {"ZZ", eos_dir+"ZZTo4L_13TeV-amcatnloFXFX-pythia8/"+batch},
	 {"WWW", eos_dir+""+batch},
	 {"WWZ", eos_dir+""+batch},
	 {"WZZ", eos_dir+""+batch},
	 {"ZZZ", eos_dir+""+batch},
	 {"tZq", eos_dir+"tZq_ll_4f_13TeV-amcatnlo-pythia8/"+batch},
	 {"TTTT", eos_dir+"TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/"+batch},
	 {"TTJets_ll",},
	 {"TTJets_lt",},
	 {"TTJets_ltbar",},
	 {"DYJets_M10to50", eos_dir+"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/"+batch},
	 {"DYJets_M50", eos_dir+"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/"+batch},
	 {"WJets",eos_dir+"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/"+batch},
	 {"ST_sLep",},
	 {"ST_tT",},
	 {"ST_tTbar",},
	 {"ST_tWT",},
	 {"ST_tWTbar",},
	 {"WWds", eos_dir+""+batch},
	 
	 {"DoubleMuon", "/eos/uscms/store/user/ztao/ttH_80X/DoubleMuon/"+batch},
	 {"SingleMuon", "/eos/uscms/store/user/ztao/ttH_80X/SingleMuon/"+batch},
	 {"DoubleEG", "/eos/uscms/store/user/ztao/ttH_80X/DoubleEG/"+batch},
	 {"SingleElectron", "/eos/uscms/store/user/ztao/ttH_80X/SingleElectron/"+batch},
	 {"MuonEG", "/eos/uscms/store/user/ztao/ttH_80X/MuonEG/"+batch}
	};

#endif
