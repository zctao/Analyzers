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
	  {"TTZ", {"TTZ", "TTGJets"}},
	  {"EWK", {"WZ", "ZZ", "WW", "WG", "ZG"}},
	  {"Rares", {"TTTT", "tZq", "WZZ", "WWZ", "WWW", "ZZZ"}},
	  {"fakes_data", DataSamples},
	  {"flips_data", DataSamples},
	  {"data_obs", DataSamples}
	};

/////////////////////////////////////////////
// Crab job output directories for different samples
TString batch = "feb2017/";
//TString data_lumi = "12_9fb/";
TString data_lumi = "36_8fb/";
TString eos_dir = "/eos/uscms/store/user/ztao/ttH_80X/";
const std::map<TString, TString> dir_map =
	{
		//{"ttH", eos_dir+"ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/"+batch},
		{"ttH", eos_dir+"ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix/"+batch},
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
		{"WWW", eos_dir+"WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/"+batch},
		{"WWZ", eos_dir+"WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/"+batch},
		{"WZZ", eos_dir+"WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/"+batch},
		{"ZZZ", eos_dir+"ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/"+batch},
		{"tZq", eos_dir+"tZq_ll_4f_13TeV-amcatnlo-pythia8/"+batch},
		{"TTTT", eos_dir+"TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/"+batch},
		{"TTJets_ll",eos_dir+""+batch},
		{"TTJets_lt",eos_dir+""+batch},
		{"TTJets_ltbar",eos_dir+""+batch},
		{"DYJets_M10to50", eos_dir+"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/"+batch},
		{"DYJets_M50", eos_dir+"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/"+batch},
		{"WJets",eos_dir+"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/"+batch},
		{"ST_sLep",eos_dir+""+batch},
		{"ST_tT",eos_dir+""+batch},
		{"ST_tTbar",eos_dir+""+batch},
		{"ST_tWT",eos_dir+""+batch},
		{"ST_tWTbar",eos_dir+""+batch},
		{"WWds", eos_dir+""+batch}, 
		{"DoubleMuon", "/eos/uscms/store/user/ztao/ttH_80X/DoubleMuon/"+batch+data_lumi},
		{"SingleMuon", "/eos/uscms/store/user/ztao/ttH_80X/SingleMuon/"+batch+data_lumi},
		{"DoubleEG", "/eos/uscms/store/user/ztao/ttH_80X/DoubleEG/"+batch+data_lumi},
		{"SingleElectron", "/eos/uscms/store/user/ztao/ttH_80X/SingleElectron/"+batch+data_lumi},
		{"MuonEG", "/eos/uscms/store/user/ztao/ttH_80X/MuonEG/"+batch+data_lumi}
	};

#endif
