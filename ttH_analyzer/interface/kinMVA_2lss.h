#ifndef kinMVA_2lss_h
#define kinMVA_2lss_h

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "Analyzers/ttH_analyzer/interface/CU_ttH_EDA_event_vars.h"
#include "Analyzers/ttH_analyzer/interface/miniLepton.h"
#include "TMVA/Reader.h"
#include "TLorentzVector.h"

#include <cmath>
#include <algorithm>

using namespace std;

class kinMVA_2lss
{
 public:

	kinMVA_2lss() {}
	~kinMVA_2lss() {}
	
	float Get_max_lep_eta() const {return max_lep_eta;}
	float Get_MT_met_lep1() const {return MT_met_lep1;}
	float Get_mindr_lep1_jet() const {return mindr_lep1_jet;}
	float Get_mindr_lep2_jet() const {return mindr_lep2_jet;}
	float Get_nJet25() const {return nJet25;}
	
	float DeltaR(double eta1, double phi1, double eta2, double phi2)
	{return sqrt( (eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2) );}

	virtual void Set_up_Reader(TMVA::Reader *) = 0;
	virtual void Calculate_mvaVars(const std::vector<miniLepton>&,
								   const std::vector<pat::Tau>&,
								   const std::vector<pat::Jet>&,
								   const pat::MET&) = 0;

	double Get_mvaScore()
	{
		if (allVarSet)
			return reader->EvaluateMVA("BDTG method");
		else
			return -9999.;
	}
	
 protected:

	TMVA::Reader *reader;
	
	float max_lep_eta;
	float MT_met_lep1;
	float mindr_lep1_jet;
	float mindr_lep2_jet;
	float nJet25;
	bool allVarSet;  // a flag indicates if all MVA variables are set
	
	float Set_max_lep_eta(double, double);
	float Set_MT_met_lep1(double, double, const pat::MET&);
	float mindr_lep_jet(float, float, const vector<pat::Jet>&);
	float Set_nJet25(const vector<pat::Jet>&);
	
};

#endif
