#if defined(__ROOTCLING__) || defined(__ACLIC__)

#ifndef TreeAnalyzer_h
#define TreeAnalyzer_h

#include "TROOT.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TString.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "Analyzers/ttH_analyzer/interface/Types_enum.h"
#include "Analyzers/ttH_analyzer/interface/SFHelper.h"
#include "Analyzers/ttH_analyzer/interface/eventNtuple.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

using namespace std;

class TreeAnalyzer
{
 public:
	
	// constructor and destructor
	TreeAnalyzer(TTree*, Analysis_types, Selection_types, bool, int verbosity=0);
	~TreeAnalyzer();
	
	// member function
	void fill_Datacards_MC(std::map<TString, TH1D*>&);
	void fill_Datacards_Data(vector<TH1D*>&, vector<vector<unsigned long long>>&);
	vector<TH1D*> makeHistograms(bool, vector<vector<unsigned long long>>&);
	void dump_Events(TString);
	void dump_Events(TString, vector<vector<unsigned long long>>&);
	
 private:

	void getWeights();
	void updateWeights();
    void buildFourVectors();
	bool checkBits(unsigned int, unsigned int);
	bool passTriggers();
	bool passFilters();
	bool passAdditionalSelection(bool cr = false);
	void triggerBits_Decoder(unsigned int);
	void filterBits_Decoder(unsigned int);
	
	TTree* _tree;
	SFHelper* _sfhelper;
	Analysis_types _AnaType;
	Selection_types _SelType;
	bool _isdata;
	int _verbosity;  // 0=unverbose; 1=show passed events; 2=show failed events;

	eventNtuple _ntuple;
	// updated weights and scale factors
	float _event_weight;
	float _PU_weight;
	float _triggerSF_weight;
	float _tauSF_weight;
	float _tauSF_weight_FRjt_normUp;
	float _tauSF_weight_FRjt_normDown;
	float _tauSF_weight_FRjt_shapeUp;
	float _tauSF_weight_FRjt_shapeDown;
	float _leptonSF_weight;
	float _MC_weight;
	float _MC_weight_scale_muF0p5;
	float _MC_weight_scale_muF2;
	float _MC_weight_scale_muR0p5;
	float _MC_weight_scale_muR2;
	float _bTagSF_weight;
	float _btagSF_weight_LFUp;
	float _btagSF_weight_LFDown;
	float _btagSF_weight_HFUp;
	float _btagSF_weight_HFDown;
	float _btagSF_weight_HFStats1Up;
	float _btagSF_weight_HFStats1Down;
	float _btagSF_weight_HFStats2Up;
	float _btagSF_weight_HFStats2Down;
	float _btagSF_weight_LFStats1Up;
	float _btagSF_weight_LFStats1Down;
	float _btagSF_weight_LFStats2Up;
	float _btagSF_weight_LFStats2Down;
	float _btagSF_weight_cErr1Up;
	float _btagSF_weight_cErr1Down;
	float _btagSF_weight_cErr2Up;
	float _btagSF_weight_cErr2Down;

	float _fr_weight_FRe_normUp;
	float _fr_weight_FRe_normDown;
	float _fr_weight_FRe_ptUp;
	float _fr_weight_FRe_ptDown;
	float _fr_weight_FRe_bUp;
	float _fr_weight_FRe_bDown;
	float _fr_weight_FRe_ecUp;
	float _fr_weight_FRe_ecDown;
	float _fr_weight_FRm_normUp;
	float _fr_weight_FRm_normDown;
	float _fr_weight_FRm_ptUp;
	float _fr_weight_FRm_ptDown;
	float _fr_weight_FRm_bUp;
	float _fr_weight_FRm_bDown;
	float _fr_weight_FRm_ecUp;
	float _fr_weight_FRm_ecDown;

	// 4 vectors
	TLorentzVector _lep0;
	TLorentzVector _lep1;
	TLorentzVector _tau;
	TLorentzVector _bjet0;  // two jets with the highest csv values
	TLorentzVector _bjet1;
	std::vector<TLorentzVector> _untag_jets;

	float _leps_conept[2];
	int _leps_id[2];
	int _leps_charge[2];
	int _leps_istight[2];

	int _tau_charge;

	bool _fourvectorsbuilt;
};

#endif

#endif
