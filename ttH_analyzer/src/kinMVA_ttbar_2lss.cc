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
	reader->AddVariable("met := min(met_pt,400)", &met);
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

void kinMVA_ttbar_2lss::Calculate_mvaVars(const std::vector<miniLepton>& leptons,
										  const std::vector<pat::Tau>& taus,
										  const std::vector<pat::Jet>& jets,
										  const pat::MET& MET)
{
	// initialize all variables to be calculated
	MT_met_lep1 = -9999.;
	mindr_lep1_jet = -9999.;
	mindr_lep2_jet = -9999.;
	max_lep_eta = -9999.;
	met = -9999.;
	nJet25 = -9999.;
	avg_dr_jet = -9999.;
	allVarSet = false;

	// calculation
	if (leptons.size() > 0) {
		MT_met_lep1 = Set_MT_met_lep1(leptons[0].conePt(), leptons[0].phi(), MET);
		mindr_lep1_jet = mindr_lep_jet(leptons[0].eta(), leptons[0].phi(), jets);
	}

	if (leptons.size() > 1) {
		// input leptons should already be sorted by conept
		assert(leptons[0].conePt() > leptons[1].conePt());
		
		mindr_lep2_jet = mindr_lep_jet(leptons[1].eta(), leptons[1].phi(), jets);
		max_lep_eta = Set_max_lep_eta(leptons[0].eta(),leptons[1].eta());
	}

	met = Set_MET(MET.pt());
	nJet25 = Set_nJet25(jets);
	avg_dr_jet = Set_avg_dr_jet(jets);

	// check if all variables are set
	allVarSet =
		MT_met_lep1!=-9999. and mindr_lep1_jet != -9999. and
		mindr_lep2_jet != -9999. and max_lep_eta != -9999. and
		met != -9999. and nJet25 != -9999. and avg_dr_jet != -9999.;
	
}

float kinMVA_ttbar_2lss::Set_avg_dr_jet(const vector<pat::Jet>& jets)
{
	if (jets.size()<2) return -9999.;

	float sum_dr_jet = 0.;
	int ncomb = 0;

	for (auto i = jets.begin(); i != jets.end()-1; ++i) {
		for (auto j = i+1; j != jets.end(); ++j) {
			++ncomb;
			sum_dr_jet += reco::deltaR(i->eta(), i->phi(), j->eta(), j->phi());
		}
	}

	return sum_dr_jet / ncomb;
}

#endif
