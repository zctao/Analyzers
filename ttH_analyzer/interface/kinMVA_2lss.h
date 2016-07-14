#ifndef kinMVA_2lss_h
#define kinMVA_2lss_h

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "Analyzers/ttH_analyzer/interface/CU_ttH_EDA_event_vars.h"
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
	
	float Get_max_lep_eta() {return max_lep_eta;}
	float Get_MT_met_lep1() {return MT_met_lep1;}
	float Get_mindr_lep1_jet() {return mindr_lep1_jet;}
	float Get_mindr_lep2_jet() {return mindr_lep2_jet;}
	float Get_nJet25() {return nJet25;}
	float Get_lep1_conePt() {return lep1_conePt;}
	float Get_lep2_conePt() {return lep2_conePt;}
	
	float DeltaR(double eta1, double phi1, double eta2, double phi2)
	{return sqrt( (eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2) );}

	//template <typename T> float Compute_lep_conePt(const T&);
	float Compute_lep_conePt(const pat::Muon&);
	float Compute_lep_conePt(const pat::Electron&);

	virtual void Set_up_Reader(TMVA::Reader *) = 0;
	virtual void Calculate_mvaVars(const CU_ttH_EDA_event_vars&, int) = 0;

	double Get_mvaScore() {return reader->EvaluateMVA("BDTG method");}
	
 protected:

	TMVA::Reader *reader;
	
	float max_lep_eta;
	float MT_met_lep1;
	float mindr_lep1_jet;
	float mindr_lep2_jet;
	float nJet25;
	float lep1_conePt, lep2_conePt;
	float lep1_eta, lep2_eta;
	float lep1_phi, lep2_phi;

	void Set_max_lep_eta(double, double);
	void Set_MT_met_lep1(double, double, const pat::MET&);
	void Set_mindr_lep_jet(float&, float, float, const vector<pat::Jet>&);
	void Set_nJet25(const vector<pat::Jet>&);

	// Get the kinematic variables of the leading and trailing leptons
	void assign_lep_kinVars(const CU_ttH_EDA_event_vars&,
							float&,float&,float&,float&,float&,float&);
	
};

#endif
