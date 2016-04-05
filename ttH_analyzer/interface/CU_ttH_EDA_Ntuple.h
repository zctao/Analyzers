#ifndef CU_ttH_EDA_Ntuple_h
#define CU_ttH_EDA_Ntuple_h

#include <map>

#include "TTree.h"
#include "TClass.h"

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
	int nEvent;
	int ls;   // luminosity section number
	int run;  // run number
	int n_presel_mu;
	int n_cutsel_mu;
	int n_mvasel_mu;
	int n_fakeablesel_mu;
	int n_presel_ele;
	int n_cutsel_ele;
	int n_mvasel_ele;
	int n_fakeablesel_ele;
	int n_presel_tau;
	int n_presel_jet;
	// event level MVA
	double MVA_2lss_ttV;
	double MVA_2lss_ttbar;
	double MT_met_lep0;
	int    n_jet25_recl;
	double mindr_lep0_jet;
	double mindr_lep1_jet;
	double lep0_conePt;
	double lep1_conePt;
	double avg_dr_jet;
	// muons
	double mu0_pt;
	double mu0_eta;
	double mu0_phi;
	double mu0_E;
	int    mu0_charge;
	int    mu0_jetNDauChargedMVASel;
	double mu0_miniRelIso;
	double mu0_miniIsoCharged;
	double mu0_miniIsoNeutral;
	double mu0_jetPtRel;
	double mu0_jetPtRatio;
	double mu0_jetCSV;
	double mu0_sip3D;
	double mu0_dxy;
	double mu0_dz;
	double mu0_segmentCompatibility;
	double mu0_leptonMVA;
	double mu0_mediumID;
	double mu0_dpt_div_pt;
	int    mu0_iscutsel;
	int    mu0_ismvasel;
	int    mu0_isfakeablesel;
	double mu1_pt;
	double mu1_eta;
	double mu1_phi;
	double mu1_E;
	int    mu1_charge;
	int    mu1_jetNDauChargedMVASel;
	double mu1_miniRelIso;
	double mu1_miniIsoCharged;
	double mu1_miniIsoNeutral;
	double mu1_jetPtRel;
	double mu1_jetPtRatio;
	double mu1_jetCSV;
	double mu1_sip3D;
	double mu1_dxy;
	double mu1_dz;
	double mu1_segmentCompatibility;
	double mu1_leptonMVA;
	double mu1_mediumID;
	double mu1_dpt_div_pt;
	int    mu1_iscutsel;
	int    mu1_ismvasel;
	int    mu1_isfakeablesel;
	
	// electrons
	double ele0_pt;
	double ele0_eta;
	double ele0_phi;
	double ele0_E;
	int    ele0_charge;
	int    ele0_jetNDauChargedMVASel;
	double ele0_miniRelIso;
	double ele0_miniIsoCharged;
	double ele0_miniIsoNeutral;
	double ele0_jetPtRel;
	double ele0_jetPtRatio;
	double ele0_jetCSV;
	double ele0_sip3D;
	double ele0_dxy;
	double ele0_dz;
	double ele0_ntMVAeleID;
	double ele0_leptonMVA;
	int    ele0_isChargeConsistent;
	int    ele0_passesConversionVeto;
	int    ele0_nMissingHits;
	int    ele0_iscutsel;
	int    ele0_ismvasel;
	int    ele0_isfakeablesel;
	double ele1_pt;
	double ele1_eta;
	double ele1_phi;
	double ele1_E;
	int    ele1_charge;
	int    ele1_jetNDauChargedMVASel;
	double ele1_miniRelIso;
	double ele1_miniIsoCharged;
	double ele1_miniIsoNeutral;
	double ele1_jetPtRel;
	double ele1_jetPtRatio;
	double ele1_jetCSV;
	double ele1_sip3D;
	double ele1_dxy;
	double ele1_dz;
	double ele1_ntMVAeleID;
	double ele1_leptonMVA;
	int    ele1_isChargeConsistent;
	int    ele1_passesConversionVeto;
	int    ele1_nMissingHits;
	int    ele1_iscutsel;
	int    ele1_ismvasel;
	int    ele1_isfakeablesel;
	
	// taus
	double tau0_pt;
	double tau0_eta;
	double tau0_phi;
	double tau0_E;
	int    tau0_charge;
	double tau0_dxy;
	double tau0_dz;
	int    tau0_decayModeFindingOldDMs;
	int    tau0_decayModeFindingNewDMs;
	double tau0_byCombinedIsolationDeltaBetaCorr3Hits;
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
	double tau1_pt;
	double tau1_eta;
	double tau1_phi;
	double tau1_E;
	int    tau1_charge;
	double tau1_dxy;
	double tau1_dz;
	int    tau1_decayModeFindingOldDMs;
	int    tau1_decayModeFindingNewDMs;
	double tau1_byCombinedIsolationDeltaBetaCorr3Hits;
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
	// jets
	double jet0_pt;
	double jet0_eta;
	double jet0_phi;
	double jet0_E;
	double jet0_CSV;
	double jet1_pt;
	double jet1_eta;
	double jet1_phi;
	double jet1_E;
	double jet1_CSV;
	double jet2_pt;
	double jet2_eta;
	double jet2_phi;
	double jet2_E;
	double jet2_CSV;
	double jet3_pt;
	double jet3_eta;
	double jet3_phi;
	double jet3_E;
	double jet3_CSV;
	// MET
	double PFMET;
	double PFMETphi;
	double MHT;
	double metLD;
	
	//ClassDef(CU_ttH_EDA_Ntuple,1);
	
};
	
#endif
	
