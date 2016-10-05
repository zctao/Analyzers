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
	
	/// Setting up jets
	//if (!config["Jets"]["jet_corrector"])
	//	throw std::runtime_error(error_message);
	//jet_corrector = config["Jets"]["jet_corrector"].as<string>();

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

	if (conf_analysis_type == "1l2tau") {
		analysis_type = Analyze_1l2tau;
		return;
	}

	if (conf_analysis_type == "2lss1tau") {
		analysis_type = Analyze_2lss1tau;
		return;
	}

	if (conf_analysis_type == "3leptons") {
		analysis_type = Analyze_3l;
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

	if (analysis_type == Analyze_1l2tau) {
		miniAODhelper.SetUp(MAODHelper_era, MAODHelper_sample_nr,
							analysisType::TauLJ, // LJ, DIL, TauLJ, TauDIL
							data				 // is data
							);
		return;
	}

	if (analysis_type == Analyze_2lss1tau) {
		miniAODhelper.SetUp(MAODHelper_era, MAODHelper_sample_nr,
							analysisType::TauDIL, // LJ, DIL, TauLJ, TauDIL
							data				  // is data
							);
		return;
	}
}

#endif
