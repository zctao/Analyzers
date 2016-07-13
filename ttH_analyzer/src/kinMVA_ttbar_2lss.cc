#ifndef kinMVA_ttbar_2lss_cc
#define kinMVA_ttbar_2lss_cc

#include "Analyzers/ttH_analyzer/interface/kinMVA_ttbar_2lss.h"

// Constructor
kinMVA_ttbar_2lss::kinMVA_ttbar_2lss()
{
	reader = new TMVA::Reader("!Color:!Silent");
	
	Set_up_Reader(reader);
}

// Destructor
kinMVA_ttbar_2lss::~kinMVA_ttbar_2lss()
{
	delete reader;
}

// Member functions
void kinMVA_ttbar_2lss::Set_up_Reader(TMVA::Reader * reader)
// Need to set up reader only once in the constructor
{
	reader->AddVariable("max_Lep_eta := max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))", &max_lep_eta);
	reader->AddVariable("numJets_float := nJet25_Recl", &nJet25);
	reader->AddVariable("mindr_lep1_jet", &mindr_lep1_jet);
	reader->AddVariable("mindr_lep2_jet", &mindr_lep2_jet);
	reader->AddVariable("met := min(met_pt,400)", &MET);
	reader->AddVariable("avg_dr_jet", &avg_dr_jet);
	reader->AddVariable("MT_met_lep1", &MT_met_lep1);

	// Spectators
	float iF0 = 0;
	float iF1 = 1;
	float iF2 = 2;
	reader->AddSpectator("iF_Recl[0]", &iF0);
	reader->AddSpectator("iF_Recl[1]", &iF1);
	reader->AddSpectator("iF_Recl[2]", &iF2);

	const std::string base =
		std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/";

	reader->BookMVA("BDTG method", base + "2lss_ttbar_BDTG.weights.xml");
}

void kinMVA_ttbar_2lss::Calculate_mvaVars(const CU_ttH_EDA_event_vars& event)
// call Calculate_MVAvars function every event after selection
{
	float lep1_eta, lep2_eta;
	float lep1_phi, lep2_phi;

	assign_lep_kinVars(event, lep1_conePt, lep2_conePt, lep1_eta,
					   lep2_eta, lep1_phi, lep2_phi);

	Set_max_lep_eta(lep1_eta,lep2_eta);
	
	Set_MT_met_lep1(lep1_conePt, lep1_phi, event.pfMET);
	
	Set_mindr_lep_jet(mindr_lep1_jet, lep1_eta, lep1_phi,
					  event.jets_selected_sorted);
	
	Set_mindr_lep_jet(mindr_lep2_jet, lep2_eta, lep2_phi,
					  event.jets_selected_sorted);
	
	Set_nJet25(event.jets_selected_sorted);
	
	Set_MET(event.pfMET.pt());
	
	Set_avg_dr_jet(event.jets_selected_sorted);
	
}

void kinMVA_ttbar_2lss::Set_avg_dr_jet(const vector<pat::Jet>& jets)
{
	assert(jets.size()>=2);

	float sum_dr_jet = 0.;
	int ncomb = 0;

	for (auto i = jets.begin(); i != jets.end()-1; ++i) {
		for (auto j = i+1; j != jets.end(); ++j) {
			++ncomb;
			sum_dr_jet += DeltaR(i->eta(), i->phi(), j->eta(), j->phi());
		}
	}

	avg_dr_jet = sum_dr_jet / ncomb;
}

#endif
