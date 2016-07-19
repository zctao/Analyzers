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

void kinMVA_ttV_2lss::Calculate_mvaVars(const std::vector<miniLepton>& leptons,
										const std::vector<pat::Tau>& taus,
										const std::vector<pat::Jet>& jets,
										const pat::MET& MET)
// call Calculate_MVAvars function every event after selection
{
	// make sure there two and only two selected leptons
	assert(leptons.size() == 2);
	// input leptons should already be sorted by conept
	assert(leptons[0].conePt() > leptons[1].conePt());

	Set_max_lep_eta(leptons[0].eta(),leptons[1].eta());

	Set_MT_met_lep1(leptons[0].conePt(), leptons[0].phi(), MET);

	mindr_lep1_jet = mindr_lep_jet(leptons[0].eta(), leptons[0].phi(), jets);
	
	mindr_lep2_jet = mindr_lep_jet(leptons[1].eta(), leptons[1].phi(), jets);
	
	Set_nJet25(jets);

	lep1_conePt = leptons[0].conePt();
	lep2_conePt = leptons[1].conePt();
	
}


#endif
