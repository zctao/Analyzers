#ifndef kinMVA_ttV_2lss_cc
#define kinMVA_ttV_2lss_cc

#include "Analyzers/ttH_analyzer/interface/kinMVA_ttV_2lss.h"

// Constructor
kinMVA_ttV_2lss::kinMVA_ttV_2lss()
{
	reader = new TMVA::Reader("!Color:!Silent");

	Set_up_Reader(reader);
}

// Destructor
kinMVA_ttV_2lss::~kinMVA_ttV_2lss()
{
	delete reader;
}

// Member functions
void kinMVA_ttV_2lss::Set_up_Reader(TMVA::Reader *reader)
{
	reader->AddVariable("max_Lep_eta := max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))", &max_lep_eta);
	reader->AddVariable("MT_met_lep1", &MT_met_lep1);
	reader->AddVariable("numJets_float := nJet25_Recl", &nJet25);
	reader->AddVariable("mindr_lep1_jet", &mindr_lep1_jet);
	reader->AddVariable("mindr_lep2_jet", &mindr_lep2_jet);
	reader->AddVariable("LepGood_conePt[iF_Recl[0]]", &lep1_conePt);
	reader->AddVariable("LepGood_conePt[iF_Recl[1]]", &lep2_conePt);

	// Spectators
	float iF0 = 0;
	float iF1 = 1;
	float iF2 = 2;
    reader->AddSpectator("iF_Recl[0]", &iF0);
	reader->AddSpectator("iF_Recl[1]", &iF1);
	reader->AddSpectator("iF_Recl[2]", &iF2);

	const std::string base =
		std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/";

	reader->BookMVA("BDTG method", base + "2lss_ttV_BDTG.weights.xml");
}

void kinMVA_ttV_2lss::Calculate_mvaVars(const CU_ttH_EDA_event_vars& event,
										int sysType = 0)
// call Calculate_MVAvars function every event after selection
{
	float lep1_eta, lep2_eta;
	float lep1_phi, lep2_phi;

	vector<pat::Jet> vjets;
	if (sysType == 1) {        // JESUp
		vjets = event.jets_selected_sorted_jesup;
	}
	else if (sysType == -1) {  // JESDown
		vjets = event.jets_selected_sorted_jesdown;
	}
	else {                     // NA
		vjets = event.jets_selected_sorted;
	}
	
	assign_lep_kinVars(event, lep1_conePt, lep2_conePt, lep1_eta,
					   lep2_eta, lep1_phi, lep2_phi);
	
	Set_max_lep_eta(lep1_eta,lep2_eta);
	
	Set_MT_met_lep1(lep1_conePt, lep1_phi, event.pfMET);
	
	Set_mindr_lep_jet(mindr_lep1_jet, lep1_eta, lep1_phi, vjets);
	
	Set_mindr_lep_jet(mindr_lep2_jet, lep2_eta, lep2_phi, vjets);
	
	Set_nJet25(vjets);
}


#endif
