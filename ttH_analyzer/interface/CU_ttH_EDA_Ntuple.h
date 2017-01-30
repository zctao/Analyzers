#ifndef CU_ttH_EDA_Ntuple_h
#define CU_ttH_EDA_Ntuple_h

#include <map>
#include <algorithm>

#include "TTree.h"
#include "TClass.h"
#include "TLorentzVector.h"

#include "Analyzers/ttH_analyzer/interface/CU_ttH_EDA_event_vars.h"

/*
 *
 * Ntuple class
 *
 */
//#ifdef __CINT__
//#pragma link C++ class CU_ttH_EDA_Ntuple+;
//#endif

class CU_ttH_EDA_Ntuple //: public TClass
{
	
 private:
	
	// private member functions
	void fill_ntuple_electrons(const std::vector<pat::Electron> &);
	void fill_ntuple_muons(const std::vector<pat::Muon> &);
	void fill_ntuple_taus(const std::vector<pat::Tau> &);
	void fill_ntuple_jets(const std::vector<pat::Jet> &);
	//void fill_ntuple_met(const pat::MET &);
		
 public:
	/// function member
    CU_ttH_EDA_Ntuple();
	~CU_ttH_EDA_Ntuple();
	
	void initialize();
	void set_up_branches(TTree *);
	void write_ntuple(const CU_ttH_EDA_event_vars &);
	
	/// variables
	// event variables
	unsigned long long nEvent;
	int ls;   // luminosity section number
	int run;  // run number
	float event_weight;
	float PU_weight;
	float MC_weight;
	float bTagSF_weight; //csv_weight;
	float leptonSF_weight;
	float tauSF_weight;
	float triggerSF_weight;//hltSF;
	float FR_weight;

	/////////////////////////////
	// systematics
	float MC_weight_scale_muF0p5;
	float MC_weight_scale_muF2;
	float MC_weight_scale_muR0p5;
	float MC_weight_scale_muR2;
	float btagSF_weight_LFUp;
	float btagSF_weight_LFDown;
	float btagSF_weight_HFUp;
	float btagSF_weight_HFDown;
	float btagSF_weight_HFStats1Up;
	float btagSF_weight_HFStats1Down;
	float btagSF_weight_HFStats2Up;
	float btagSF_weight_HFStats2Down;
	float btagSF_weight_LFStats1Up;
	float btagSF_weight_LFStats1Down;
	float btagSF_weight_LFStats2Up;
	float btagSF_weight_LFStats2Down;
	float btagSF_weight_cErr1Up;
	float btagSF_weight_cErr1Down;
	float btagSF_weight_cErr2Up;
	float btagSF_weight_cErr2Down;
	/////////////////////////////
	
	int isGenMatched;
	int HiggsDecayType;   // Higgs decay product pdgId

	int lepCategory;   // 0: mumu; 1: ee; 2: emu
	int btagCategory;  // 0: loose; 1: medium (>=2 medium btags)
	
	float npuTrue;
	float npuInTime;
	
	int pass_single_mu;
	int pass_single_e;
	int pass_double_mu;
	int pass_double_e;
	int pass_elemu;
	int matchHLTPath;
	// save hlt paths
	int HLT_Ele27_WPTight_Gsf;
	int HLT_IsoMu24;
	int HLT_IsoTkMu24;
	int HLT_IsoMu22_eta2p1;
	int HLT_IsoTkMu22_eta2p1;
	int HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
	int HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL;
	int HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
	int HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
	
	// event level MVA
	float MVA_2lss_ttV;
	float MVA_2lss_ttbar;
	float MT_met_lep0;
	int    n_jet25_recl;
	float mindr_lep0_jet;
	float mindr_lep1_jet;
	float lep0_conept;
	float lep1_conept;
	float avg_dr_jet;
	float max_lep_eta;

	int ibin;  // bin index in 1D BDT shape template 
	
	int n_presel_mu;
	int n_mvasel_mu;
	int n_fakeablesel_mu;
	int n_presel_ele;
	int n_mvasel_ele;
	int n_fakeablesel_ele;
	int n_presel_tau;
	int n_tau;
	int n_presel_jet;
	
	// muons
	float mu0_pt;
	float mu0_conept;
	float mu0_eta;
	float mu0_phi;
	float mu0_E;
	int    mu0_charge;
	int    mu0_jetNDauChargedMVASel;
	float mu0_miniRelIso;
	float mu0_miniIsoCharged;
	float mu0_miniIsoNeutral;
	float mu0_jetPtRel;
	float mu0_jetPtRatio;
	float mu0_jetCSV;
	float mu0_sip3D;
	float mu0_dxy;
	float mu0_dz;
	float mu0_segmentCompatibility;
	float mu0_leptonMVA;
	float mu0_mediumID;
	float mu0_dpt_div_pt;
	int    mu0_ismvasel;
	int    mu0_isfakeablesel;
	int mu0_mcMatchType;
	float mu1_pt;
	float mu1_conept;
	float mu1_eta;
	float mu1_phi;
	float mu1_E;
	int    mu1_charge;
	int    mu1_jetNDauChargedMVASel;
	float mu1_miniRelIso;
	float mu1_miniIsoCharged;
	float mu1_miniIsoNeutral;
	float mu1_jetPtRel;
	float mu1_jetPtRatio;
	float mu1_jetCSV;
	float mu1_sip3D;
	float mu1_dxy;
	float mu1_dz;
	float mu1_segmentCompatibility;
	float mu1_leptonMVA;
	float mu1_mediumID;
	float mu1_dpt_div_pt;
	int    mu1_ismvasel;
	int    mu1_isfakeablesel;
	int mu1_mcMatchType;
	
	// electrons
	float ele0_pt;
	float ele0_conept;
	float ele0_eta;
	float ele0_phi;
	float ele0_E;
	int    ele0_charge;
	int    ele0_jetNDauChargedMVASel;
	float ele0_miniRelIso;
	float ele0_miniIsoCharged;
	float ele0_miniIsoNeutral;
	float ele0_jetPtRel;
	float ele0_jetPtRatio;
	float ele0_jetCSV;
	float ele0_sip3D;
	float ele0_dxy;
	float ele0_dz;
	float ele0_ntMVAeleID;
	float ele0_leptonMVA;
	int    ele0_isChargeConsistent;
	int    ele0_passesConversionVeto;
	int    ele0_nMissingHits;
	int    ele0_ismvasel;
	int    ele0_isfakeablesel;
	int ele0_mcMatchType;
	float ele1_pt;
	float ele1_conept;
	float ele1_eta;
	float ele1_phi;
	float ele1_E;
	int    ele1_charge;
	int    ele1_jetNDauChargedMVASel;
	float ele1_miniRelIso;
	float ele1_miniIsoCharged;
	float ele1_miniIsoNeutral;
	float ele1_jetPtRel;
	float ele1_jetPtRatio;
	float ele1_jetCSV;
	float ele1_sip3D;
	float ele1_dxy;
	float ele1_dz;
	float ele1_ntMVAeleID;
	float ele1_leptonMVA;
	int    ele1_isChargeConsistent;
	int    ele1_passesConversionVeto;
	int    ele1_nMissingHits;
	int    ele1_ismvasel;
	int    ele1_isfakeablesel;
	int ele1_mcMatchType;
	
	// taus
	float tau0_pt;
	float tau0_eta;
	float tau0_phi;
	float tau0_E;
	int    tau0_charge;
	float tau0_dxy;
	float tau0_dz;
	int    tau0_decayMode;
	int    tau0_decayModeFindingOldDMs;
	int    tau0_decayModeFindingNewDMs;
	float tau0_byCombinedIsolationDeltaBetaCorr3Hits;
	int    tau0_byLooseCombinedIsolationDeltaBetaCorr3Hits;
	int    tau0_byMediumCombinedIsolationDeltaBetaCorr3Hits;
	int    tau0_byTightCombinedIsolationDeltaBetaCorr3Hits;
	int    tau0_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03;
	int    tau0_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03;
	int    tau0_byTightCombinedIsolationDeltaBetaCorr3HitsdR03;
	int    tau0_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
	int    tau0_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
	int    tau0_byTightIsolationMVArun2v1DBdR03oldDMwLT;
	int    tau0_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
	int    tau0_againstMuonLoose3;
	int    tau0_againstMuonTight3;
	int    tau0_againstElectronVLooseMVA6;
	int    tau0_againstElectronLooseMVA6;
	int    tau0_againstElectronMediumMVA6;
	int    tau0_againstElectronTightMVA6;
	int tau0_mcMatchType;
	float tau1_pt;
	float tau1_eta;
	float tau1_phi;
	float tau1_E;
	int    tau1_charge;
	float tau1_dxy;
	float tau1_dz;
	int    tau1_decayMode;
	int    tau1_decayModeFindingOldDMs;
	int    tau1_decayModeFindingNewDMs;
	float tau1_byCombinedIsolationDeltaBetaCorr3Hits;
	int    tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits;
	int    tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits;
	int    tau1_byTightCombinedIsolationDeltaBetaCorr3Hits;
	int    tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03;
	int    tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03;
	int    tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03;
	int    tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
	int    tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
	int    tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT;
	int    tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
	int    tau1_againstMuonLoose3;
	int    tau1_againstMuonTight3;
	int    tau1_againstElectronVLooseMVA6;
	int    tau1_againstElectronLooseMVA6;
	int    tau1_againstElectronMediumMVA6;
	int    tau1_againstElectronTightMVA6;
	int tau1_mcMatchType;
	// jets
	float jet0_pt;
	float jet0_eta;
	float jet0_phi;
	float jet0_E;
	float jet0_CSV;
	float jet1_pt;
	float jet1_eta;
	float jet1_phi;
	float jet1_E;
	float jet1_CSV;
	float jet2_pt;
	float jet2_eta;
	float jet2_phi;
	float jet2_E;
	float jet2_CSV;
	float jet3_pt;
	float jet3_eta;
	float jet3_phi;
	float jet3_E;
	float jet3_CSV;
	// MET
	float PFMET;
	float PFMETphi;
	float MHT;
	float metLD;
	
	//ClassDef(CU_ttH_EDA_Ntuple,1);
	
};

#endif
	
