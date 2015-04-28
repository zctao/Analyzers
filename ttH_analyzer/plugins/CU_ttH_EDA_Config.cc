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
	if (!config["Triggers"]["HLT_electron_triggers"])
		throw std::runtime_error(error_message);
	trigger_on_HLT_electrons =
		config["Triggers"]["HLT_electron_triggers"].as<string>();

	if (!config["Triggers"]["HLT_muon_triggers"])
		throw std::runtime_error(error_message);
	trigger_on_HLT_muons = config["Triggers"]["HLT_muon_triggers"].as<string>();

	/// Setting up cuts
	if (!config["Cuts"]["min_tight_lepton_pT"])
		throw std::runtime_error(error_message);
	min_tight_lepton_pT = config["Cuts"]["min_tight_lepton_pT"].as<double>();

	/// Setting up jets
	if (!config["Jets"]["jet_corrector"])
		throw std::runtime_error(error_message);
	jet_corrector = config["Jets"]["jet_corrector"].as<string>();

	/// Setting up miniAODhelper
	if (!config["miniAODhelper_parameters"]["analysis_type"] ||
		!config["miniAODhelper_parameters"]["using_real_data"])
		throw std::runtime_error(error_message);
	Load_configuration_MAODH(
		config["miniAODhelper_parameters"]["analysis_type"].as<string>(),
		config["miniAODhelper_parameters"]["using_real_data"].as<bool>());
}

void CU_ttH_EDA::Load_configuration_MAODH(const string &analysis_type,
										  bool data)
{
	if (analysis_type == "lepton+jet") {
		miniAODhelper.SetUp(MAODHelper_era, MAODHelper_sample_nr,
							analysisType::LJ, // LJ, DIL, TauLJ, TauDIL
							data			  // is data
							);
		return;
	}

	if (analysis_type == "dilepton") {
		miniAODhelper.SetUp(MAODHelper_era, MAODHelper_sample_nr,
							analysisType::DIL, // LJ, DIL, TauLJ, TauDIL
							data			   // is data
							);
		return;
	}

	if (analysis_type == "taus_lepton+jet") {
		miniAODhelper.SetUp(MAODHelper_era, MAODHelper_sample_nr,
							analysisType::TauLJ, // LJ, DIL, TauLJ, TauDIL
							data				 // is data
							);
		return;
	}

	if (analysis_type == "taus_dilepton") {
		miniAODhelper.SetUp(MAODHelper_era, MAODHelper_sample_nr,
							analysisType::TauDIL, // LJ, DIL, TauLJ, TauDIL
							data				  // is data
							);
		return;
	}
}

#endif