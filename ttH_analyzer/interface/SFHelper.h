#ifndef SFHelper_h
#define SFHelper_h

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "Analyzers/ttH_analyzer/interface/miniLepton.h"
#endif

#include "Analyzers/ttH_analyzer/interface/Types_enum.h"
// Root
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"

#include <iostream>
#include <assert.h>

class SFHelper
{
 public:
	// constructor and destructor
	SFHelper(Analysis_types, Selection_types, bool);
	~SFHelper();

	// member functions
	float Get_HLTSF(int);
	float Get_LeptonIDSF(float,float,bool,bool,bool);	
	float Get_PUWeight(int);	
	float Get_TauIDSF(float,float,bool);
	//float Get_MCWeight();
	float Get_FakeRate(float,float,bool,bool); // for ele or mu
	float Get_FakeRate(float,float);  // for tau
	//float Get_ChargeFlipWeight();
	float Get_EleChargeMisIDProb(float,float,int,int);
	
#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
	float Get_LeptonIDSF(const miniLepton&);
	float Get_EvtCSVWeight(std::vector<pat::Jet> &, const std::string &);
	float Get_TauIDSF(const pat::Tau&, bool);
	float Get_FakeRate(const miniLepton&);
	float Get_FakeRate(const pat::Tau&);
	float Get_EleChargeMisIDProb(const miniLepton&, int);
#endif
		
	// utilities
	float read2DHist(TH2*, float, float);
	float evalTGraph(TGraphAsymmErrors*, float);
	float readTGraph(TGraphAsymmErrors*, float);
	float readTF(TF1*, float);
	
 private:

	bool _isdata;
	Analysis_types _analysis;
	Selection_types _selection;

	//float _hlt_sf;      // trigger scale factor
	//float _lepIDEff_sf; // lepton ID efficiency scale factor
	//float _csv_weight;  // btag csv re-weighting
	//float _pu_weight;   // pileup weight
	//float _tauID_sf;    // tau ID scale factor
	//float _mc_weight;   // gen weight
	//float _FR_weight;   // fake rate weight;
	//float _QF_weight;   // charge flip weight;

	// Fake lepton rate 
	TFile *file_fr_lep;
	TH2F *h_fakerate_el;
	TH2F *h_fakerate_mu;

	// Jet to tau fake rate
	TFile* file_fr_tau;
	TGraphAsymmErrors *g_fakerate_tau_mvaM_etaL_mc;
	TGraphAsymmErrors *g_fakerate_tau_mvaM_etaH_mc;
	TF1 *f_fakerate_tau_mvaM_etaL_ratio;
	TF1 *f_fakerate_tau_mvaM_etaH_ratio;

	// PU weight
	TFile* file_puweight;
	TH1F *h_puweight;
	
	// Lepton ID scale factor lookup tables
	TFile* file_recoToLoose_leptonSF_mu1_b;
	TFile* file_recoToLoose_leptonSF_mu1_e;
	TFile* file_recoToLoose_leptonSF_mu2;
	TFile* file_recoToLoose_leptonSF_mu3;
	TFile* file_recoToLoose_leptonSF_el;
	TFile* file_recoToLoose_leptonSF_gsf;
	TFile* file_looseToTight_leptonSF_mu_2lss;
	TFile* file_looseToTight_leptonSF_el_2lss;
	TFile* file_looseToTight_leptonSF_mu_3l;
	TFile* file_looseToTight_leptonSF_el_3l;

	TGraphAsymmErrors *h_recoToLoose_leptonSF_mu1_b;
	TGraphAsymmErrors *h_recoToLoose_leptonSF_mu1_e;
	TH2F *h_recoToLoose_leptonSF_mu2;
	TGraphAsymmErrors *h_recoToLoose_leptonSF_mu3;
	TH2F *h_recoToLoose_leptonSF_el1;
	TH2F *h_recoToLoose_leptonSF_el2;
	TH2F *h_recoToLoose_leptonSF_el3;
	TH2F *h_recoToLoose_leptonSF_gsf;
	TH2F *h_looseToTight_leptonSF_mu_2lss;
	TH2F *h_looseToTight_leptonSF_el_2lss;
	TH2F *h_looseToTight_leptonSF_mu_3l;
	TH2F *h_looseToTight_leptonSF_el_3l;

	// Electron Charge MisID
	TFile *file_eleMisCharge;
	TH2F *h_chargeMisId;
	
	// CSV
#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
	BTagCalibrationReader* BTagCaliReader;
	//std::string Btag_sysList[16] =
	//	{"LFUp","LFDown","HFUp","HFDown",
	//	 "HFStats1Up","HFStats1Down","HFStats2Up","HFStats2Down",
	//	 "LFStats1Up","LFStats1Down","LFStats2Up","LFStats2Down",
	//	 "cErr1Up","cErr1Down","cErr2Up","cErr2Down"};
	void Set_up_BTagCalibration_Readers();
	void Delete_BTagCalibration_Readers();
#endif
	
	void Set_up_FakeRate_Lut();
	void Set_up_TauSF_Lut();
	void Set_up_ChargeMisID_Lut();
	void Set_up_PUWeight_hist();
	void Set_up_LeptonSF_Lut();

	void Delete_FakeRate_Lut();
	void Delete_TauSF_Lut();
	void Delete_ChargeMisID_Lut();
	void Delete_PUWeight_hist();
	void Delete_LeptonSF_Lut();


	float Get_LeptonSF_loose(float,float,bool,bool);
	float Get_LeptonSF_tight_vs_loose(float,float,bool,bool);
#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
	float Get_LeptonSF_loose(const miniLepton&);
	float Get_LeptonSF_tight_vs_loose(const miniLepton&);
	float Get_JetCSVWeight(pat::Jet&, std::string);
#endif
	
};

#endif
