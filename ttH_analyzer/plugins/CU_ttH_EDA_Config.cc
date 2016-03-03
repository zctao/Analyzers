#ifndef CU_ttH_EDA_Config_cc
#define CU_ttH_EDA_Config_cc

/// Includes
#include "CU_ttH_EDA.h"

void CU_ttH_EDA::Load_configuration(string config_filename)
{
	const string error_message = "A problem in "
								 "CU_ttH_EDA::Load_configuration(..) occured "
								 "while reading a config.";

	YAML::Node config = YAML::LoadFile(config_filename);

	/// Setting up analysis type
	if (!config["Generic"]["analysis_type"])
		throw std::runtime_error(error_message);
	Load_configuration_set_type(
		config["Generic"]["analysis_type"].as<string>());

	/// Setting up generic parameters
	if (!config["Generic"]["HLT_config_tag"])
		throw std::runtime_error(error_message);
	hltTag = config["Generic"]["HLT_config_tag"].as<string>();

	if (!config["Generic"]["filter_config_tag"])
		throw std::runtime_error(error_message);
	filterTag = config["Generic"]["filter_config_tag"].as<string>();

	if (!config["Generic"]["verbosity"])
		throw std::runtime_error(error_message);
	verbose_ = config["Generic"]["verbosity"].as<bool>();

	if (!config["Generic"]["print_HLT_event_path"])
		throw std::runtime_error(error_message);
	dumpHLT_ = config["Generic"]["print_HLT_event_path"].as<bool>();

	/// Setting up HLT triggers
	if (!config["Triggers"]["collect_trigger_stats"])
		throw std::runtime_error(error_message);
	trigger_stats = config["Triggers"]["collect_trigger_stats"].as<bool>();

	if (!config["Triggers"]["HLT_electron_triggers"])
		throw std::runtime_error(error_message);
	trigger_on_HLT_e =
		config["Triggers"]["HLT_electron_triggers"].as<vector<string>>();

	if (!config["Triggers"]["HLT_muon_triggers"])
		throw std::runtime_error(error_message);
	trigger_on_HLT_mu =
		config["Triggers"]["HLT_muon_triggers"].as<vector<string>>();

	if (!config["Triggers"]["HLT_electron_electron_triggers"])
		throw std::runtime_error(error_message);
	trigger_on_HLT_ee = config["Triggers"]["HLT_electron_electron_triggers"]
							.as<vector<string>>();

	if (!config["Triggers"]["HLT_electron_muon_triggers"])
		throw std::runtime_error(error_message);
	trigger_on_HLT_emu =
		config["Triggers"]["HLT_electron_muon_triggers"].as<vector<string>>();

	if (!config["Triggers"]["HLT_muon_muon_triggers"])
		throw std::runtime_error(error_message);
	trigger_on_HLT_mumu =
		config["Triggers"]["HLT_muon_muon_triggers"].as<vector<string>>();

	/// Setting up cuts
	if (!config["Cuts"]["min_tight_lepton_pT"])
		throw std::runtime_error(error_message);
	min_tight_lepton_pT = config["Cuts"]["min_tight_lepton_pT"].as<float>();

	if (!config["Cuts"]["min_tight_tau_pT"])
		throw std::runtime_error(error_message);
	min_tight_tau_pT = config["Cuts"]["min_tight_tau_pT"].as<float>();
	
	if (!config["Cuts"]["min_jet_pT"])
		throw std::runtime_error(error_message);
	min_jet_pT = config["Cuts"]["min_jet_pT"].as<float>();

	if (!config["Cuts"]["min_bjet_pT"])
		throw std::runtime_error(error_message);
	min_bjet_pT = config["Cuts"]["min_bjet_pT"].as<float>();

	if (!config["Cuts"]["max_jet_eta"])
		throw std::runtime_error(error_message);
	max_jet_eta = config["Cuts"]["max_jet_eta"].as<float>();

	if (!config["Cuts"]["max_bjet_eta"])
		throw std::runtime_error(error_message);
	max_bjet_eta = config["Cuts"]["max_bjet_eta"].as<float>();

	if (!config["Cuts"]["min_njets"])
		throw std::runtime_error(error_message);
	min_njets = config["Cuts"]["min_njets"].as<int>();

	if (!config["Cuts"]["min_nbtags"])
		throw std::runtime_error(error_message);
	min_nbtags = config["Cuts"]["min_nbtags"].as<int>();
	
	/// Setting up jets
	if (!config["Jets"]["jet_corrector"])
		throw std::runtime_error(error_message);
	jet_corrector = config["Jets"]["jet_corrector"].as<string>();

	/// Setting up miniAODhelper
	if (!config["miniAODhelper_parameters"]["using_real_data"])
		throw std::runtime_error(error_message);
	Load_configuration_MAODH(
		config["miniAODhelper_parameters"]["using_real_data"].as<bool>());

	if (!config["miniAODhelper_parameters"]["b_tag_strength"])
		throw std::runtime_error(error_message);
	MAODHelper_b_tag_strength =
		config["miniAODhelper_parameters"]["b_tag_strength"].as<char>();
}

void CU_ttH_EDA::Load_configuration_set_type(const string &conf_analysis_type)
{
	if (conf_analysis_type == "lepton+jet") {
		analysis_type = Analyze_lepton_jet;
		return;
	}

	if (conf_analysis_type == "dilepton") {
		analysis_type = Analyze_dilepton;
		return;
	}

	if (conf_analysis_type == "taus_lepton+jet") {
		analysis_type = Analyze_ditaus_lepton;
		return;
	}

	if (conf_analysis_type == "taus_dilepton") {
		analysis_type = Analyze_tau_ssleptons;
		return;
	}
}

void CU_ttH_EDA::Load_configuration_MAODH(bool data)
{
	if (analysis_type == Analyze_lepton_jet) {
		miniAODhelper.SetUp(MAODHelper_era, MAODHelper_sample_nr,
							analysisType::LJ, // LJ, DIL, TauLJ, TauDIL
							data			  // is data
							);
		return;
	}

	if (analysis_type == Analyze_dilepton) {
		miniAODhelper.SetUp(MAODHelper_era, MAODHelper_sample_nr,
							analysisType::DIL, // LJ, DIL, TauLJ, TauDIL
							data			   // is data
							);
		return;
	}

	if (analysis_type == Analyze_ditaus_lepton) {
		miniAODhelper.SetUp(MAODHelper_era, MAODHelper_sample_nr,
							analysisType::TauLJ, // LJ, DIL, TauLJ, TauDIL
							data				 // is data
							);
		return;
	}

	if (analysis_type == Analyze_tau_ssleptons) {
		miniAODhelper.SetUp(MAODHelper_era, MAODHelper_sample_nr,
							analysisType::TauDIL, // LJ, DIL, TauLJ, TauDIL
							data				  // is data
							);
		return;
	}
}

#endif
