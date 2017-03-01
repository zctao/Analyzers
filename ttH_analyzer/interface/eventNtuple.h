#ifndef eventNtuple_h
#define eventNtuple_h

#include "TTree.h"

class eventNtuple
{
 public:
	
	eventNtuple(){};
	~eventNtuple(){};

	void set_branch_address(TTree*);
	
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
	// trigger and filter bits
	unsigned int triggerBits;
	unsigned int filterBits;

	int nBadMuons;
	
	// event level MVA
	float MVA_2lss_ttV;
	float MVA_2lss_ttbar;
	float MT_met_lep0;
	float mindr_lep0_jet;
	float mindr_lep1_jet;
	float lep0_conept;
	float lep1_conept;
	float avg_dr_jet;

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
	int mu0_isPFMuon;
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
	int mu1_isPFMuon;
	
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
	// vectors of jets
	std::vector<float> *jets_pt = 0;
	std::vector<float> *jets_eta = 0;
	std::vector<float> *jets_phi = 0;
	std::vector<float> *jets_E = 0;
	std::vector<float> *jets_csv = 0;
	std::vector<int> *jets_flavor = 0;
	// MET
	float PFMET;
	float PFMETphi;
	float MHT;
	float metLD;
	float METSignificance;
	float METCov00;
	float METCov10;
	float METCov01;
	float METCov11;
};

void eventNtuple::set_branch_address(TTree* tree)
{
	tree->SetBranchAddress("run", &run);
	tree->SetBranchAddress("ls", &ls);
	tree->SetBranchAddress("nEvent", &nEvent);
	tree->SetBranchAddress("event_weight", &event_weight);
	tree->SetBranchAddress("PU_weight", &PU_weight);
	tree->SetBranchAddress("MC_weight", &MC_weight);
	tree->SetBranchAddress("bTagSF_weight", &bTagSF_weight);
	tree->SetBranchAddress("leptonSF_weight", &leptonSF_weight);
	tree->SetBranchAddress("tauSF_weight", &tauSF_weight);
	tree->SetBranchAddress("triggerSF_weight", &triggerSF_weight);
	tree->SetBranchAddress("FR_weight", &FR_weight);
	tree->SetBranchAddress("MC_weight_scale_muF0p5", &MC_weight_scale_muF0p5);
	tree->SetBranchAddress("MC_weight_scale_muF2", &MC_weight_scale_muF2);
	tree->SetBranchAddress("MC_weight_scale_muR0p5", &MC_weight_scale_muR0p5);
	tree->SetBranchAddress("MC_weight_scale_muR2", &MC_weight_scale_muR2);
	tree->SetBranchAddress("btagSF_weight_LFUp", &btagSF_weight_LFUp);
	tree->SetBranchAddress("btagSF_weight_LFDown", &btagSF_weight_LFDown);
	tree->SetBranchAddress("btagSF_weight_HFUp", &btagSF_weight_HFUp);
	tree->SetBranchAddress("btagSF_weight_HFDown", &btagSF_weight_HFDown);
	tree->SetBranchAddress("btagSF_weight_HFStats1Up", &btagSF_weight_HFStats1Up);
	tree->SetBranchAddress("btagSF_weight_HFStats1Down", &btagSF_weight_HFStats1Down);
	tree->SetBranchAddress("btagSF_weight_HFStats2Up", &btagSF_weight_HFStats2Up);
	tree->SetBranchAddress("btagSF_weight_HFStats2Down", &btagSF_weight_HFStats2Down);
	tree->SetBranchAddress("btagSF_weight_LFStats1Up", &btagSF_weight_LFStats1Up);
	tree->SetBranchAddress("btagSF_weight_LFStats1Down", &btagSF_weight_LFStats1Down);
	tree->SetBranchAddress("btagSF_weight_LFStats2Up", &btagSF_weight_LFStats2Up);
	tree->SetBranchAddress("btagSF_weight_LFStats2Down", &btagSF_weight_LFStats2Down);
	tree->SetBranchAddress("btagSF_weight_cErr1Up", &btagSF_weight_cErr1Up);
	tree->SetBranchAddress("btagSF_weight_cErr1Down", &btagSF_weight_cErr1Down);
	tree->SetBranchAddress("btagSF_weight_cErr2Up", &btagSF_weight_cErr2Up);
	tree->SetBranchAddress("btagSF_weight_cErr2Down", &btagSF_weight_cErr2Down);
	tree->SetBranchAddress("isGenMatched", &isGenMatched);
	tree->SetBranchAddress("HiggsDecayType", &HiggsDecayType);
	tree->SetBranchAddress("lepCategory", &lepCategory);
	tree->SetBranchAddress("btagCategory", &btagCategory);
	tree->SetBranchAddress("npuTrue", &npuTrue);
	tree->SetBranchAddress("npuInTime", &npuInTime);
	tree->SetBranchAddress("pass_single_mu", &pass_single_mu);
	tree->SetBranchAddress("pass_single_e", &pass_single_e);
	tree->SetBranchAddress("pass_double_mu", &pass_double_mu);
	tree->SetBranchAddress("pass_double_e", &pass_double_e);
	tree->SetBranchAddress("pass_elemu", &pass_elemu);
	tree->SetBranchAddress("matchHLTPath", &matchHLTPath);
	tree->SetBranchAddress("triggerBits", &triggerBits);
	tree->SetBranchAddress("filterBits", &filterBits);
	tree->SetBranchAddress("nBadMuons", &nBadMuons);	
	tree->SetBranchAddress("MVA_2lss_ttbar", &MVA_2lss_ttbar);
	tree->SetBranchAddress("MVA_2lss_ttV", &MVA_2lss_ttbar);
	tree->SetBranchAddress("MT_met_lep0", &MT_met_lep0);
	tree->SetBranchAddress("mindr_lep0_jet", &mindr_lep0_jet);
	tree->SetBranchAddress("mindr_lep1_jet", &mindr_lep1_jet);
	tree->SetBranchAddress("lep0_conept", &lep0_conept);
	tree->SetBranchAddress("lep1_conept", &lep1_conept);
	tree->SetBranchAddress("avg_dr_jet", &avg_dr_jet);
	tree->SetBranchAddress("ibin", &ibin);
	tree->SetBranchAddress("n_presel_mu", &n_presel_mu);
	tree->SetBranchAddress("n_mvasel_mu", &n_mvasel_mu);
	tree->SetBranchAddress("n_fakeablesel_mu", &n_fakeablesel_mu);
	tree->SetBranchAddress("n_presel_ele", &n_presel_ele);
	tree->SetBranchAddress("n_mvasel_ele", &n_mvasel_ele);
	tree->SetBranchAddress("n_fakeablesel_ele", &n_fakeablesel_ele);
	tree->SetBranchAddress("n_presel_tau", &n_presel_tau);
	tree->SetBranchAddress("n_tau", &n_tau);
	tree->SetBranchAddress("n_presel_jet", &n_presel_jet);
	// muons
	tree->SetBranchAddress("mu0_pt",                   &mu0_pt);
	tree->SetBranchAddress("mu0_conept",               &mu0_conept);
	tree->SetBranchAddress("mu0_eta",                  &mu0_eta);
	tree->SetBranchAddress("mu0_phi",                  &mu0_phi);
	tree->SetBranchAddress("mu0_E",                    &mu0_E);
	tree->SetBranchAddress("mu0_charge",               &mu0_charge);
	tree->SetBranchAddress("mu0_jetNDauChargedMVASel", &mu0_jetNDauChargedMVASel);
	tree->SetBranchAddress("mu0_miniRelIso",           &mu0_miniRelIso);
	tree->SetBranchAddress("mu0_miniIsoCharged",       &mu0_miniIsoCharged);
	tree->SetBranchAddress("mu0_miniIsoNeutral",       &mu0_miniIsoNeutral);
	tree->SetBranchAddress("mu0_jetPtRel",             &mu0_jetPtRel);
	tree->SetBranchAddress("mu0_jetPtRatio",           &mu0_jetPtRatio);
	tree->SetBranchAddress("mu0_jetCSV",               &mu0_jetCSV);
	tree->SetBranchAddress("mu0_sip3D",                &mu0_sip3D);
	tree->SetBranchAddress("mu0_dxy",                  &mu0_dxy);
	tree->SetBranchAddress("mu0_dz",                   &mu0_dz);
	tree->SetBranchAddress("mu0_segmentCompatibility", &mu0_segmentCompatibility);
	tree->SetBranchAddress("mu0_leptonMVA",            &mu0_leptonMVA);
	tree->SetBranchAddress("mu0_mediumID",             &mu0_mediumID);
	tree->SetBranchAddress("mu0_dpt_div_pt",           &mu0_dpt_div_pt);
	tree->SetBranchAddress("mu0_ismvasel",             &mu0_ismvasel);
	tree->SetBranchAddress("mu0_isfakeablesel",        &mu0_isfakeablesel);
	tree->SetBranchAddress("mu0_mcMatchType",          &mu0_mcMatchType);
	tree->SetBranchAddress("mu0_isPFMuon",             &mu0_isPFMuon);
	tree->SetBranchAddress("mu1_pt",                   &mu1_pt);
	tree->SetBranchAddress("mu1_conept",               &mu1_conept);
	tree->SetBranchAddress("mu1_eta",                  &mu1_eta);
	tree->SetBranchAddress("mu1_phi",                  &mu1_phi);
	tree->SetBranchAddress("mu1_E",                    &mu1_E);
	tree->SetBranchAddress("mu1_charge",               &mu1_charge);
	tree->SetBranchAddress("mu1_jetNDauChargedMVASel", &mu1_jetNDauChargedMVASel);
	tree->SetBranchAddress("mu1_miniRelIso",           &mu1_miniRelIso);
	tree->SetBranchAddress("mu1_miniIsoCharged",       &mu1_miniIsoCharged);
	tree->SetBranchAddress("mu1_miniIsoNeutral",       &mu1_miniIsoNeutral);
	tree->SetBranchAddress("mu1_jetPtRel",             &mu1_jetPtRel);
	tree->SetBranchAddress("mu1_jetPtRatio",           &mu1_jetPtRatio);
	tree->SetBranchAddress("mu1_jetCSV",               &mu1_jetCSV);
	tree->SetBranchAddress("mu1_sip3D",                &mu1_sip3D);
	tree->SetBranchAddress("mu1_dxy",                  &mu1_dxy);
	tree->SetBranchAddress("mu1_dz",                   &mu1_dz);
	tree->SetBranchAddress("mu1_segmentCompatibility", &mu1_segmentCompatibility);
	tree->SetBranchAddress("mu1_leptonMVA",            &mu1_leptonMVA);
	tree->SetBranchAddress("mu1_mediumID",             &mu1_mediumID);
	tree->SetBranchAddress("mu1_dpt_div_pt",           &mu1_dpt_div_pt);
	tree->SetBranchAddress("mu1_ismvasel",             &mu1_ismvasel);
	tree->SetBranchAddress("mu1_isfakeablesel",        &mu1_isfakeablesel);
	tree->SetBranchAddress("mu1_mcMatchType",          &mu1_mcMatchType);
	tree->SetBranchAddress("mu1_isPFMuon",             &mu1_isPFMuon);
	// electron
	tree->SetBranchAddress("ele0_pt",                   &ele0_pt);
	tree->SetBranchAddress("ele0_conept",               &ele0_conept);
	tree->SetBranchAddress("ele0_eta",                  &ele0_eta);
	tree->SetBranchAddress("ele0_phi",                  &ele0_phi);
	tree->SetBranchAddress("ele0_E",                    &ele0_E);
	tree->SetBranchAddress("ele0_charge",               &ele0_charge);
	tree->SetBranchAddress("ele0_jetNDauChargedMVASel", &ele0_jetNDauChargedMVASel);
	tree->SetBranchAddress("ele0_miniRelIso",           &ele0_miniRelIso);
	tree->SetBranchAddress("ele0_miniIsoCharged",       &ele0_miniIsoCharged);
	tree->SetBranchAddress("ele0_miniIsoNeutral",       &ele0_miniIsoNeutral);
	tree->SetBranchAddress("ele0_jetPtRel",             &ele0_jetPtRel);
	tree->SetBranchAddress("ele0_jetPtRatio",           &ele0_jetPtRatio);
	tree->SetBranchAddress("ele0_jetCSV",               &ele0_jetCSV);
	tree->SetBranchAddress("ele0_sip3D",                &ele0_sip3D);
	tree->SetBranchAddress("ele0_dxy",                  &ele0_dxy);
	tree->SetBranchAddress("ele0_dz",                   &ele0_dz);
	tree->SetBranchAddress("ele0_ntMVAeleID",           &ele0_ntMVAeleID);
	tree->SetBranchAddress("ele0_leptonMVA",            &ele0_leptonMVA);
	tree->SetBranchAddress("ele0_isChargeConsistent", &ele0_isChargeConsistent);
	tree->SetBranchAddress("ele0_passesConversionVeto", &ele0_passesConversionVeto);
	tree->SetBranchAddress("ele0_nMissingHits", &ele0_nMissingHits);
	tree->SetBranchAddress("ele0_ismvasel", &ele0_ismvasel);
	tree->SetBranchAddress("ele0_isfakeablesel", &ele0_isfakeablesel);
	tree->SetBranchAddress("ele0_mcMatchType", &ele0_mcMatchType);
	tree->SetBranchAddress("ele1_pt",                   &ele1_pt);
	tree->SetBranchAddress("ele1_conept",               &ele1_conept);
	tree->SetBranchAddress("ele1_eta",                  &ele1_eta);
	tree->SetBranchAddress("ele1_phi",                  &ele1_phi);
	tree->SetBranchAddress("ele1_E",                    &ele1_E);
	tree->SetBranchAddress("ele1_charge",               &ele1_charge);
	tree->SetBranchAddress("ele1_jetNDauChargedMVASel", &ele1_jetNDauChargedMVASel);
	tree->SetBranchAddress("ele1_miniRelIso",           &ele1_miniRelIso);
	tree->SetBranchAddress("ele1_miniIsoCharged",       &ele1_miniIsoCharged);
	tree->SetBranchAddress("ele1_miniIsoNeutral",       &ele1_miniIsoNeutral);
	tree->SetBranchAddress("ele1_jetPtRel",             &ele1_jetPtRel);
	tree->SetBranchAddress("ele1_jetPtRatio",           &ele1_jetPtRatio);
	tree->SetBranchAddress("ele1_jetCSV",               &ele1_jetCSV);
	tree->SetBranchAddress("ele1_sip3D",                &ele1_sip3D);
	tree->SetBranchAddress("ele1_dxy",                  &ele1_dxy);
	tree->SetBranchAddress("ele1_dz",                   &ele1_dz);
	tree->SetBranchAddress("ele1_ntMVAeleID",           &ele1_ntMVAeleID);
	tree->SetBranchAddress("ele1_leptonMVA",            &ele1_leptonMVA);
	tree->SetBranchAddress("ele1_isChargeConsistent", &ele1_isChargeConsistent);
	tree->SetBranchAddress("ele1_passesConversionVeto", &ele1_passesConversionVeto);
	tree->SetBranchAddress("ele1_nMissingHits", &ele1_nMissingHits);
	tree->SetBranchAddress("ele1_ismvasel", &ele1_ismvasel);
	tree->SetBranchAddress("ele1_isfakeablesel", &ele1_isfakeablesel);
	tree->SetBranchAddress("ele1_mcMatchType", &ele1_mcMatchType);
	// taus
	tree->SetBranchAddress("tau0_pt", &tau0_pt);
	tree->SetBranchAddress("tau0_eta", &tau0_eta);
	tree->SetBranchAddress("tau0_phi", &tau0_phi);
	tree->SetBranchAddress("tau0_E", &tau0_E);
	tree->SetBranchAddress("tau0_charge", &tau0_charge);
	tree->SetBranchAddress("tau0_dxy", &tau0_dxy);
	tree->SetBranchAddress("tau0_dz", &tau0_dz);
	tree->SetBranchAddress("tau0_decayMode", &tau0_decayMode);
	tree->SetBranchAddress("tau0_mcMatchType", &tau0_mcMatchType);
	// jets
    tree->SetBranchAddress("jets_pt", &jets_pt);
	tree->SetBranchAddress("jets_eta", &jets_eta);
	tree->SetBranchAddress("jets_phi", &jets_phi);
	tree->SetBranchAddress("jets_E", &jets_E);
	tree->SetBranchAddress("jets_csv", &jets_csv);
	tree->SetBranchAddress("jets_flavor", &jets_flavor);
	// MET
	tree->SetBranchAddress("PFMET", &PFMET);
	tree->SetBranchAddress("PFMETphi", &PFMETphi);
	tree->SetBranchAddress("METSignificance", &METSignificance);
	tree->SetBranchAddress("METCov00", &METCov00);
	tree->SetBranchAddress("METCov01", &METCov01);
	tree->SetBranchAddress("METCov10", &METCov10);
	tree->SetBranchAddress("METCov11", &METCov11);
	tree->SetBranchAddress("MHT", &MHT);
}

#endif
