#ifndef CU_ttH_EDA_Ntuple_cc
#define CU_ttH_EDA_Ntuple_cc

#include "Analyzers/ttH_analyzer/interface/CU_ttH_EDA_Ntuple.h"
#include "Analyzers/ttH_analyzer/interface/CU_ttH_EDA_event_vars.h"
#include <cmath>

// Constructor
CU_ttH_EDA_Ntuple::CU_ttH_EDA_Ntuple() {}

// Destructor
CU_ttH_EDA_Ntuple::~CU_ttH_EDA_Ntuple() {}

// Member functions

void CU_ttH_EDA_Ntuple::write_ntuple(const CU_ttH_EDA_event_vars &local)
{
	// Event variales
	nEvent = local.event_nr;
	ls = local.lumisection_nr;
	run = local.run_nr;
	event_weight = local.weight;
	genWeight = local.gen_weight;
	
	// Muons
	n_presel_mu = local.n_muons_loose;
	n_fakeablesel_mu = local.n_muons_fakeable;
	n_mvasel_mu = local.n_muons_tight;

	fill_ntuple_muons(local.mu_preselected_sorted);
	//fill_ntuple_muons(local.mu_fakeable);
	
	// Electrons
	n_presel_ele = local.n_electrons_loose;
	n_fakeablesel_ele = local.n_electrons_fakeable;
	n_mvasel_ele = local.n_electrons_tight;
	
	fill_ntuple_electrons(local.e_preselected_sorted);
	//fill_ntuple_electrons(local.e_fakeable);
		
	// Taus
	n_presel_tau = local.n_taus_pre;
	n_tau = local.n_taus;
	fill_ntuple_taus(local.tau_preselected_sorted);

	// Jets
	n_presel_jet = local.n_jets;
	fill_ntuple_jets(local.jets_selected_sorted);
	
	// MET/MHT
	PFMET = local.pfMET.pt();
	PFMETphi = local.pfMET.phi();
	MHT = local.MHT;
	metLD = local.metLD;
}

void CU_ttH_EDA_Ntuple::fill_ntuple_muons(const std::vector<pat::Muon>& muons)
{
	if (muons.size() > 0) {
		mu0_pt = muons[0].pt();
		mu0_eta = muons[0].eta();
		mu0_phi = muons[0].phi();
		mu0_E = muons[0].energy();
		mu0_charge = muons[0].charge();
		mu0_jetNDauChargedMVASel = muons[0].userFloat("nearestJetNDauCharged");
		mu0_miniRelIso = muons[0].userFloat("miniIso");
		mu0_miniIsoCharged = muons[0].userFloat("miniAbsIsoCharged");
		mu0_miniIsoNeutral = muons[0].userFloat("miniAbsIsoNeutralcorr");
		mu0_jetPtRel = muons[0].userFloat("nearestJetPtRel");
		mu0_jetPtRatio = muons[0].userFloat("nearestJetPtRatio");
		mu0_jetCSV = muons[0].userFloat("nearestJetCsv");
		mu0_sip3D = muons[0].userFloat("sip3D");
		mu0_dxy = muons[0].userFloat("dxy");
		mu0_dz = muons[0].userFloat("dz");
		mu0_segmentCompatibility = muons[0].segmentCompatibility();
		mu0_leptonMVA = muons[0].userFloat("leptonMVA");
		//mu0_mediumID = muons[0].isMediumMuon();
		////////////////////
		/// short term HIP safe muon mediumID
		bool goodGlb0 = (muons[0].isGlobalMuon() && muons[0].userFloat("normalizedChiSq") < 3 && muons[0].userFloat("localChiSq") < 12 && muons[0].userFloat("trackKink") < 20);
		mu0_mediumID = (muons[0].userFloat("validFraction") >= 0.8 && muons[0].segmentCompatibility() >= (goodGlb0 ? 0.303 : 0.451));
		////////////////////
		mu0_dpt_div_pt = muons[0].muonBestTrack()->ptError()/muons[0].muonBestTrack()->pt();
		mu0_isfakeablesel = muons[0].userFloat("idFakeable") > 0.5;
		mu0_ismvasel = muons[0].userFloat("idMVABased") > 0.5;
	}
	
	if (muons.size() > 1 ) {
		mu1_pt = muons[1].pt();
		mu1_eta = muons[1].eta();
		mu1_phi = muons[1].phi();
		mu1_E = muons[1].energy();
		mu1_charge = muons[1].charge();
		mu1_jetNDauChargedMVASel = muons[1].userFloat("nearestJetNDauCharged");
		mu1_miniRelIso = muons[1].userFloat("miniIso");
		mu1_miniIsoCharged = muons[1].userFloat("miniAbsIsoCharged");
		mu1_miniIsoNeutral = muons[1].userFloat("miniAbsIsoNeutralcorr");
		mu1_jetPtRel = muons[1].userFloat("nearestJetPtRel");
		mu1_jetPtRatio = muons[1].userFloat("nearestJetPtRatio");
		mu1_jetCSV = muons[1].userFloat("nearestJetCsv");
		mu1_sip3D = muons[1].userFloat("sip3D");
		mu1_dxy = muons[1].userFloat("dxy");
		mu1_dz = muons[1].userFloat("dz");
		mu1_segmentCompatibility = muons[1].segmentCompatibility();
		mu1_leptonMVA = muons[1].userFloat("leptonMVA");
		//mu1_mediumID = muons[1].isMediumMuon();
		////////////////////
		/// short term HIP safe muon mediumID
		bool goodGlb1 = (muons[1].isGlobalMuon() && muons[1].userFloat("normalizedChiSq") < 3 && muons[1].userFloat("localChiSq") < 12 && muons[1].userFloat("trackKink") < 20);
		mu1_mediumID = (muons[1].userFloat("validFraction") >= 0.8 && muons[1].segmentCompatibility() >= (goodGlb1 ? 0.303 : 0.451));
		////////////////////		
		mu1_dpt_div_pt = muons[1].muonBestTrack()->ptError()/muons[1].muonBestTrack()->pt();
		mu1_isfakeablesel = muons[1].userFloat("idFakeable") > 0.5;
		mu1_ismvasel = muons[1].userFloat("idMVABased") > 0.5;
	}
}

void CU_ttH_EDA_Ntuple::fill_ntuple_electrons(const std::vector<pat::Electron>& electrons)
{
	if (electrons.size() > 0) {
		ele0_pt = electrons[0].pt();
		ele0_eta = electrons[0].eta();
		ele0_phi = electrons[0].phi();
		ele0_E = electrons[0].energy();
		ele0_charge = electrons[0].charge();
		ele0_jetNDauChargedMVASel = electrons[0].userFloat("nearestJetNDauCharged");
		ele0_miniRelIso = electrons[0].userFloat("miniIso");
		ele0_miniIsoCharged = electrons[0].userFloat("miniAbsIsoCharged");
		ele0_miniIsoNeutral = electrons[0].userFloat("miniAbsIsoNeutralcorr");
		ele0_jetPtRel = electrons[0].userFloat("nearestJetPtRel");
		ele0_jetPtRatio = electrons[0].userFloat("nearestJetPtRatio");
		ele0_jetCSV = electrons[0].userFloat("nearestJetCsv");
		ele0_sip3D = electrons[0].userFloat("sip3D");
		ele0_dxy = electrons[0].userFloat("dxy");
		ele0_dz = electrons[0].userFloat("dz");
		ele0_ntMVAeleID = electrons[0].userFloat("eleMvaId"); // non-triggering mva ID
		ele0_leptonMVA = electrons[0].userFloat("leptonMVA");
		ele0_isChargeConsistent = (electrons[0].isGsfCtfScPixChargeConsistent() + electrons[0].isGsfScPixChargeConsistent()) > 1;
		ele0_passesConversionVeto = electrons[0].passConversionVeto();
		ele0_nMissingHits = electrons[0].userFloat("numMissingHits");
		ele0_ismvasel = electrons[0].userFloat("idMVABased") > 0.5;
		ele0_isfakeablesel = electrons[0].userFloat("idFakeable") > 0.5;
	}
	
	if (electrons.size() > 1 ) {
		ele1_pt = electrons[1].pt();
		ele1_eta = electrons[1].eta();
		ele1_phi = electrons[1].phi();
		ele1_E = electrons[1].energy();
		ele1_charge = electrons[1].charge();
		ele1_jetNDauChargedMVASel = electrons[1].userFloat("nearestJetNDauCharged");
		ele1_miniRelIso = electrons[1].userFloat("miniIso");
		ele1_miniIsoCharged = electrons[1].userFloat("miniAbsIsoCharged");
		ele1_miniIsoNeutral = electrons[1].userFloat("miniAbsIsoNeutralcorr");
		ele1_jetPtRel = electrons[1].userFloat("nearestJetPtRel");
		ele1_jetPtRatio = electrons[1].userFloat("nearestJetPtRatio");
		ele1_jetCSV = electrons[1].userFloat("nearestJetCsv");
		ele1_sip3D = electrons[1].userFloat("sip3D");
		ele1_dxy = electrons[1].userFloat("dxy");
		ele1_dz = electrons[1].userFloat("dz");
		ele1_ntMVAeleID = electrons[1].userFloat("eleMvaId"); // non-triggering mva ID
		ele1_leptonMVA = electrons[1].userFloat("leptonMVA");
		ele1_isChargeConsistent = (electrons[1].isGsfCtfScPixChargeConsistent() + electrons[1].isGsfScPixChargeConsistent()) > 1;
		ele1_passesConversionVeto = electrons[1].passConversionVeto();
		ele1_nMissingHits = electrons[1].userFloat("numMissingHits");
		ele1_ismvasel = electrons[1].userFloat("idMVABased") > 0.5;
		ele1_isfakeablesel = electrons[1].userFloat("idFakeable") > 0.5;
	}
}

void CU_ttH_EDA_Ntuple::fill_ntuple_taus(const std::vector<pat::Tau>& taus)
{
	if (taus.size() > 0) {
		tau0_pt = taus[0].pt();
		tau0_eta = taus[0].eta();
		tau0_phi = taus[0].phi();
		tau0_E = taus[0].energy();
		tau0_charge = taus[0].charge();
		tau0_dxy = taus[0].userFloat("dxy");
		tau0_dz = taus[0].userFloat("dz");
		tau0_decayMode = taus[0].decayMode();
		tau0_decayModeFindingOldDMs =  taus[0].tauID("decayModeFinding");  // decayModeFindingOldDMs
		tau0_decayModeFindingNewDMs =  taus[0].tauID("decayModeFindingNewDMs");
		
		tau0_byCombinedIsolationDeltaBetaCorr3Hits = taus[0].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
		tau0_byLooseCombinedIsolationDeltaBetaCorr3Hits = taus[0].tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
		tau0_byMediumCombinedIsolationDeltaBetaCorr3Hits = taus[0].tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
		tau0_byTightCombinedIsolationDeltaBetaCorr3Hits = taus[0].tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
		
		//tau0_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 = taus[0].tauID("byLooseCombinedIsolationDeltaBetaCorr3HitsdR03");
		//tau0_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = taus[0].tauID("byMediumCombinedIsolationDeltaBetaCorr3HitsdR03");
		//tau0_byTightCombinedIsolationDeltaBetaCorr3HitsdR03 = taus[0].tauID("byTightCombinedIsolationDeltaBetaCorr3HitsdR03");
		
		tau0_byLooseIsolationMVArun2v1DBdR03oldDMwLT = taus[0].tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT");
		tau0_byMediumIsolationMVArun2v1DBdR03oldDMwLT = taus[0].tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT");
		tau0_byTightIsolationMVArun2v1DBdR03oldDMwLT = taus[0].tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT");
		tau0_byVTightIsolationMVArun2v1DBdR03oldDMwLT = taus[0].tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT");

		tau0_againstMuonLoose3 = taus[0].tauID("againstMuonLoose3");
		tau0_againstMuonTight3 = taus[0].tauID("againstMuonTight3");

		tau0_againstElectronVLooseMVA6 = taus[0].tauID("againstElectronVLooseMVA6");
		tau0_againstElectronLooseMVA6 = taus[0].tauID("againstElectronLooseMVA6");
		tau0_againstElectronMediumMVA6 = taus[0].tauID("againstElectronMediumMVA6");
		tau0_againstElectronTightMVA6 = taus[0].tauID("againstElectronTightMVA6");
	}
	
	if (taus.size() > 1 ) {
		tau1_pt = taus[1].pt();
		tau1_eta = taus[1].eta();
		tau1_phi = taus[1].phi();
		tau1_E = taus[1].energy();
		tau1_charge = taus[1].charge();
		tau1_dxy = taus[1].userFloat("dxy");
		tau1_dz = taus[1].userFloat("dz");
		tau1_decayMode = taus[1].decayMode();
		tau1_decayModeFindingOldDMs =  taus[1].tauID("decayModeFinding");  // decayModeFindingOldDMs
		tau1_decayModeFindingNewDMs =  taus[1].tauID("decayModeFindingNewDMs");
		
		tau1_byCombinedIsolationDeltaBetaCorr3Hits = taus[1].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
		tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits = taus[1].tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
		tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits = taus[1].tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
		tau1_byTightCombinedIsolationDeltaBetaCorr3Hits = taus[1].tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
		
		//tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 = taus[1].tauID("byLooseCombinedIsolationDeltaBetaCorr3HitsdR03");
		//tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = taus[1].tauID("byMediumCombinedIsolationDeltaBetaCorr3HitsdR03");
		//tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03 = taus[1].tauID("byTightCombinedIsolationDeltaBetaCorr3HitsdR03");
		
		tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT = taus[1].tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT");
		tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT = taus[1].tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT");
		tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT = taus[1].tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT");
		tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT = taus[1].tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT");

		tau1_againstMuonLoose3 = taus[1].tauID("againstMuonLoose3");
		tau1_againstMuonTight3 = taus[1].tauID("againstMuonTight3");

		tau1_againstElectronVLooseMVA6 = taus[1].tauID("againstElectronVLooseMVA6");
		tau1_againstElectronLooseMVA6 = taus[1].tauID("againstElectronLooseMVA6");
		tau1_againstElectronMediumMVA6 = taus[1].tauID("againstElectronMediumMVA6");
		tau1_againstElectronTightMVA6 = taus[1].tauID("againstElectronTightMVA6");
	}
}

void CU_ttH_EDA_Ntuple::fill_ntuple_jets(const std::vector<pat::Jet>& jets)
{
	if (jets.size() > 0) {
		jet0_pt = jets[0].pt();
		jet0_eta = jets[0].eta();
		jet0_phi = jets[0].phi();
		jet0_E = jets[0].energy();
		jet0_CSV = jets[0].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	}
	
	if (jets.size() > 1 ) {
		jet1_pt = jets[1].pt();
		jet1_eta = jets[1].eta();
		jet1_phi = jets[1].phi();
		jet1_E = jets[1].energy();
		jet1_CSV = jets[1].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	}
	
	if (jets.size() > 2) {
		jet2_pt = jets[2].pt();
		jet2_eta = jets[2].eta();
		jet2_phi = jets[2].phi();
		jet2_E = jets[2].energy();
		jet2_CSV = jets[2].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 
	}
	
	if (jets.size() > 3 ) {
		jet3_pt = jets[3].pt();
		jet3_eta = jets[3].eta();
		jet3_phi = jets[3].phi();
		jet3_E = jets[3].energy();
		jet3_CSV = jets[3].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	}
}

void CU_ttH_EDA_Ntuple::initialize()
{
	// event variables
	nEvent = -9999;
	ls = -9999;
	run = -9999;
	event_weight = -9999.;
	genWeight = -9999.;
	n_presel_mu = -9999;
	n_mvasel_mu = -9999;
	n_fakeablesel_mu = -9999;
	n_presel_ele = -9999;
	n_mvasel_ele = -9999;
	n_fakeablesel_ele = -9999;
	n_presel_tau = -9999;
	n_tau = -9999;
	n_presel_jet = -9999;

	MVA_2lss_ttV = -9999.;
	MVA_2lss_ttbar = -9999.;
	MT_met_lep0 = -9999.;
	n_jet25_recl = -9999;
	mindr_lep0_jet = -9999.;
	mindr_lep1_jet = -9999.;
	lep0_conept = -9999.;
	lep1_conept = -9999.;
	avg_dr_jet = -9999.;
	max_lep_eta = -9999.;

	// muons
	mu0_pt = -9999.;
	//mu0_conept = -9999.;
	mu0_eta = -9999.;
	mu0_phi = -9999.;
	mu0_E = -9999.;
	mu0_charge = -9999;
	mu0_jetNDauChargedMVASel = -9999;
	mu0_miniRelIso = -9999.;
	mu0_miniIsoCharged = -9999.;
	mu0_miniIsoNeutral = -9999.;
	mu0_jetPtRel = -9999.;
	mu0_jetPtRatio = -9999.;
	mu0_jetCSV = -9999.;
	mu0_sip3D = -9999.;
	mu0_dxy = -9999.;
	mu0_dz = -9999.;
	mu0_segmentCompatibility = -9999.;
	mu0_leptonMVA = -9999.;
	mu0_mediumID = -9999.;
	mu0_dpt_div_pt = -9999.;
	mu0_ismvasel = 0;//-9999;
	mu0_isfakeablesel = 0;//-9999;
	mu1_pt = -9999.;
	//mu1_conept = -9999.;
	mu1_eta = -9999.;
	mu1_phi = -9999.;
	mu1_E = -9999.;
	mu1_charge = -9999;
	mu1_jetNDauChargedMVASel = -9999;
	mu1_miniRelIso = -9999.;
	mu1_miniIsoCharged = -9999.;
	mu1_miniIsoNeutral = -9999.;
	mu1_jetPtRel = -9999.;
	mu1_jetPtRatio = -9999.;
	mu1_jetCSV = -9999.;
	mu1_sip3D = -9999.;
	mu1_dxy = -9999.;
	mu1_dz = -9999.;
	mu1_segmentCompatibility = -9999.;
	mu1_leptonMVA = -9999.;
	mu1_mediumID = -9999.;
	mu1_dpt_div_pt = -9999.;
	mu1_ismvasel = 0;//-9999;
	mu1_isfakeablesel = 0;//-9999;
	
	// electrons
	ele0_pt = -9999.;
	//ele0_conept = -9999.;
	ele0_eta = -9999.;
	ele0_phi = -9999.;
	ele0_E = -9999.;
	ele0_charge = -9999;
	ele0_jetNDauChargedMVASel = -9999;
	ele0_miniRelIso = -9999.;
	ele0_miniIsoCharged = -9999.;
	ele0_miniIsoNeutral = -9999.;
	ele0_jetPtRel = -9999.;
	ele0_jetPtRatio = -9999.;
	ele0_jetCSV = -9999.;
	ele0_sip3D = -9999.;
	ele0_dxy = -9999.;
	ele0_dz = -9999.;
	ele0_ntMVAeleID = -9999.;
	ele0_leptonMVA = -9999.;
	ele0_isChargeConsistent = 0;//-9999;
	ele0_passesConversionVeto = 0;//-9999;
	ele0_nMissingHits = -9999;
	ele0_ismvasel = 0;//-9999;
	ele0_isfakeablesel = 0;//-9999;
	ele1_pt = -9999.;
	//ele1_conept = -9999.;
	ele1_eta = -9999.;
	ele1_phi = -9999.;
	ele1_E = -9999.;
	ele1_charge = -9999;
	ele1_jetNDauChargedMVASel = -9999;
	ele1_miniRelIso = -9999.;
	ele1_miniIsoCharged = -9999.;
	ele1_miniIsoNeutral = -9999.;
	ele1_jetPtRel = -9999.;
	ele1_jetPtRatio = -9999.;
	ele1_jetCSV = -9999.;
	ele1_sip3D = -9999.;
	ele1_dxy = -9999.;
	ele1_dz = -9999.;
	ele1_ntMVAeleID = -9999.;
	ele1_leptonMVA = -9999.;
	ele1_isChargeConsistent = 0;//-9999;
	ele1_passesConversionVeto = 0;//-9999;
	ele1_nMissingHits = -9999;
	ele1_ismvasel = 0;//-9999;
	ele1_isfakeablesel = 0;//-9999;
	
	// taus
	tau0_pt = -9999.;
	tau0_eta = -9999.;
	tau0_phi = -9999.;
	tau0_E = -9999.;
	tau0_charge = -9999;
	tau0_dxy = -9999.;
	tau0_dz = -9999.;
	tau0_decayMode = -9999;
	tau0_decayModeFindingOldDMs = -9999;
	tau0_decayModeFindingNewDMs = -9999;
	tau0_byCombinedIsolationDeltaBetaCorr3Hits = -9999.;
	tau0_byLooseCombinedIsolationDeltaBetaCorr3Hits = -9999;
	tau0_byMediumCombinedIsolationDeltaBetaCorr3Hits = -9999;
	tau0_byTightCombinedIsolationDeltaBetaCorr3Hits = -9999;
	tau0_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
	tau0_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
	tau0_byTightCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
	tau0_byLooseIsolationMVArun2v1DBdR03oldDMwLT = -9999;
	tau0_byMediumIsolationMVArun2v1DBdR03oldDMwLT = -9999;
	tau0_byTightIsolationMVArun2v1DBdR03oldDMwLT = -9999;
	tau0_byVTightIsolationMVArun2v1DBdR03oldDMwLT = -9999;
	tau0_againstMuonLoose3 = -9999;
	tau0_againstMuonTight3 = -9999;
	tau0_againstElectronVLooseMVA6 = -9999;
	tau0_againstElectronLooseMVA6 = -9999;
	tau0_againstElectronMediumMVA6 = -9999;
	tau0_againstElectronTightMVA6 = -9999;
	tau1_pt = -9999.;
	tau1_eta = -9999.;
	tau1_phi = -9999.;
	tau1_E = -9999.;
	tau1_charge = -9999;
	tau1_dxy = -9999.;
	tau1_dz = -9999.;
	tau1_decayMode = -9999;
	tau1_decayModeFindingOldDMs = -9999;
	tau1_decayModeFindingNewDMs = -9999;
	tau1_byCombinedIsolationDeltaBetaCorr3Hits = -9999.;
	tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits = -9999;
	tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits = -9999;
	tau1_byTightCombinedIsolationDeltaBetaCorr3Hits = -9999;
	tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
	tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
	tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
	tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT = -9999;
	tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT = -9999;
	tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT = -9999;
	tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT = -9999;
	tau1_againstMuonLoose3 = -9999;
	tau1_againstMuonTight3 = -9999;
	tau1_againstElectronVLooseMVA6 = -9999;
	tau1_againstElectronLooseMVA6 = -9999;
	tau1_againstElectronMediumMVA6 = -9999;
	tau1_againstElectronTightMVA6 = -9999;
		
	// jets
	jet0_pt = -9999.;
	jet0_eta = -9999.;
	jet0_phi = -9999.;
	jet0_E = -9999.;
	jet0_CSV = -9999.;
	jet1_pt = -9999.;
	jet1_eta = -9999.;
	jet1_phi = -9999.;
	jet1_E = -9999.;
	jet1_CSV = -9999.;
	jet2_pt = -9999.;
	jet2_eta = -9999.;
	jet2_phi = -9999.;
	jet2_E = -9999.;
	jet2_CSV = -9999.;
	jet3_pt = -9999.;
	jet3_eta = -9999.;
	jet3_phi = -9999.;
	jet3_E = -9999.;
	jet3_CSV = -9999.;
	// MET
	PFMET = -9999.;
	PFMETphi = -9999.;
	MHT = -9999.;
	metLD = -9999.;
}

void CU_ttH_EDA_Ntuple::set_up_branches(TTree *tree)
{
	// initialize ntuple
	this->initialize();

	// set up tree branches
	tree->Branch("nEvent", &nEvent);
	tree->Branch("ls", &ls);
	tree->Branch("run", &run);
	tree->Branch("event_weight", &event_weight);
	tree->Branch("genWeight", &genWeight);
	tree->Branch("n_presel_mu", &n_presel_mu);
	tree->Branch("n_mvasel_mu", &n_mvasel_mu);
	tree->Branch("n_fakeablesel_mu", &n_fakeablesel_mu);
	tree->Branch("n_presel_ele", &n_presel_ele);
	tree->Branch("n_mvasel_ele", &n_mvasel_ele);
	tree->Branch("n_fakeablesel_ele", &n_fakeablesel_ele);
	tree->Branch("n_presel_tau", &n_presel_tau);
	tree->Branch("n_tau", &n_tau);
	tree->Branch("n_presel_jet", &n_presel_jet);
	tree->Branch("MVA_2lss_ttV", &MVA_2lss_ttV);
	tree->Branch("MVA_2lss_ttbar", &MVA_2lss_ttbar);
	tree->Branch("MT_met_lep0", &MT_met_lep0);
	//tree->Branch("n_jet25_recl", &n_jet25_recl);
	tree->Branch("mindr_lep0_jet", &mindr_lep0_jet);
	tree->Branch("mindr_lep1_jet", &mindr_lep1_jet);
	tree->Branch("lep0_conept", &lep0_conept);
	tree->Branch("lep1_conept", &lep1_conept);
	tree->Branch("avg_dr_jet", &avg_dr_jet);
	// muons
	tree->Branch("mu0_pt",                   &mu0_pt);
	//tree->Branch("mu0_conept",               &mu0_conept);
	tree->Branch("mu0_eta",                  &mu0_eta);
	tree->Branch("mu0_phi",                  &mu0_phi);
	tree->Branch("mu0_E",                    &mu0_E);
	tree->Branch("mu0_charge",               &mu0_charge);
	tree->Branch("mu0_jetNDauChargedMVASel", &mu0_jetNDauChargedMVASel);
	tree->Branch("mu0_miniRelIso",           &mu0_miniRelIso);
	tree->Branch("mu0_miniIsoCharged",       &mu0_miniIsoCharged);
	tree->Branch("mu0_miniIsoNeutral",       &mu0_miniIsoNeutral);
	tree->Branch("mu0_jetPtRel",             &mu0_jetPtRel);
	tree->Branch("mu0_jetPtRatio",           &mu0_jetPtRatio);
	tree->Branch("mu0_jetCSV",               &mu0_jetCSV);
	tree->Branch("mu0_sip3D",                &mu0_sip3D);
	tree->Branch("mu0_dxy",                  &mu0_dxy);
	tree->Branch("mu0_dz",                   &mu0_dz);
	tree->Branch("mu0_segmentCompatibility", &mu0_segmentCompatibility);
	tree->Branch("mu0_leptonMVA",            &mu0_leptonMVA);
	tree->Branch("mu0_mediumID", &mu0_mediumID);
	tree->Branch("mu0_dpt_div_pt", &mu0_dpt_div_pt);
	tree->Branch("mu0_ismvasel", &mu0_ismvasel);
	tree->Branch("mu0_isfakeablesel", &mu0_isfakeablesel);
	tree->Branch("mu1_pt",                   &mu1_pt);
	//tree->Branch("mu1_conept",               &mu1_conept);
	tree->Branch("mu1_eta",                  &mu1_eta);
	tree->Branch("mu1_phi",                  &mu1_phi);
	tree->Branch("mu1_E",                    &mu1_E);
	tree->Branch("mu1_charge",               &mu1_charge);
	tree->Branch("mu1_jetNDauChargedMVASel", &mu1_jetNDauChargedMVASel);
	tree->Branch("mu1_miniRelIso",           &mu1_miniRelIso);
	tree->Branch("mu1_miniIsoCharged",       &mu1_miniIsoCharged);
	tree->Branch("mu1_miniIsoNeutral",       &mu1_miniIsoNeutral);
	tree->Branch("mu1_jetPtRel",             &mu1_jetPtRel);
	tree->Branch("mu1_jetPtRatio",           &mu1_jetPtRatio);
	tree->Branch("mu1_jetCSV",               &mu1_jetCSV);
	tree->Branch("mu1_sip3D",                &mu1_sip3D);
	tree->Branch("mu1_dxy",                  &mu1_dxy);
	tree->Branch("mu1_dz",                   &mu1_dz);
	tree->Branch("mu1_segmentCompatibility", &mu1_segmentCompatibility);
	tree->Branch("mu1_leptonMVA",            &mu1_leptonMVA);
	tree->Branch("mu1_mediumID", &mu1_mediumID);
	tree->Branch("mu1_dpt_div_pt", &mu1_dpt_div_pt);
	tree->Branch("mu1_ismvasel", &mu1_ismvasel);
	tree->Branch("mu1_isfakeablesel", &mu1_isfakeablesel);
	// electrons
	tree->Branch("ele0_pt",                   &ele0_pt);
	//tree->Branch("ele0_conept",               &ele0_conept);
	tree->Branch("ele0_eta",                  &ele0_eta);
	tree->Branch("ele0_phi",                  &ele0_phi);
	tree->Branch("ele0_E",                    &ele0_E);
	tree->Branch("ele0_charge",               &ele0_charge);
	tree->Branch("ele0_jetNDauChargedMVASel", &ele0_jetNDauChargedMVASel);
	tree->Branch("ele0_miniRelIso",           &ele0_miniRelIso);
	tree->Branch("ele0_miniIsoCharged",       &ele0_miniIsoCharged);
	tree->Branch("ele0_miniIsoNeutral",       &ele0_miniIsoNeutral);
	tree->Branch("ele0_jetPtRel",             &ele0_jetPtRel);
	tree->Branch("ele0_jetPtRatio",           &ele0_jetPtRatio);
	tree->Branch("ele0_jetCSV",               &ele0_jetCSV);
	tree->Branch("ele0_sip3D",                &ele0_sip3D);
	tree->Branch("ele0_dxy",                  &ele0_dxy);
	tree->Branch("ele0_dz",                   &ele0_dz);
	tree->Branch("ele0_ntMVAeleID",           &ele0_ntMVAeleID);
	tree->Branch("ele0_leptonMVA",            &ele0_leptonMVA);
	tree->Branch("ele0_isChargeConsistent", &ele0_isChargeConsistent);
	tree->Branch("ele0_passesConversionVeto", &ele0_passesConversionVeto);
	tree->Branch("ele0_nMissingHits", &ele0_nMissingHits);
	tree->Branch("ele0_ismvasel", &ele0_ismvasel);
	tree->Branch("ele0_isfakeablesel", &ele0_isfakeablesel);
	tree->Branch("ele1_pt",                   &ele1_pt);
	//tree->Branch("ele1_conept",               &ele1_conept);
	tree->Branch("ele1_eta",                  &ele1_eta);
	tree->Branch("ele1_phi",                  &ele1_phi);
	tree->Branch("ele1_E",                    &ele1_E);
	tree->Branch("ele1_charge",               &ele1_charge);
	tree->Branch("ele1_jetNDauChargedMVASel", &ele1_jetNDauChargedMVASel);
	tree->Branch("ele1_miniRelIso",           &ele1_miniRelIso);
	tree->Branch("ele1_miniIsoCharged",       &ele1_miniIsoCharged);
	tree->Branch("ele1_miniIsoNeutral",       &ele1_miniIsoNeutral);
	tree->Branch("ele1_jetPtRel",             &ele1_jetPtRel);
	tree->Branch("ele1_jetPtRatio",           &ele1_jetPtRatio);
	tree->Branch("ele1_jetCSV",               &ele1_jetCSV);
	tree->Branch("ele1_sip3D",                &ele1_sip3D);
	tree->Branch("ele1_dxy",                  &ele1_dxy);
	tree->Branch("ele1_dz",                   &ele1_dz);
	tree->Branch("ele1_ntMVAeleID",           &ele1_ntMVAeleID);
	tree->Branch("ele1_leptonMVA",            &ele1_leptonMVA);
	tree->Branch("ele1_isChargeConsistent", &ele1_isChargeConsistent);
	tree->Branch("ele1_passesConversionVeto", &ele1_passesConversionVeto);
	tree->Branch("ele1_nMissingHits", &ele1_nMissingHits);
	tree->Branch("ele1_ismvasel", &ele1_ismvasel);
	tree->Branch("ele1_isfakeablesel", &ele1_isfakeablesel);
	// taus
	tree->Branch("tau0_pt", &tau0_pt);
	tree->Branch("tau0_eta", &tau0_eta);
	tree->Branch("tau0_phi", &tau0_phi);
	tree->Branch("tau0_E", &tau0_E);
	tree->Branch("tau0_charge", &tau0_charge);
	tree->Branch("tau0_dxy", &tau0_dxy);
	tree->Branch("tau0_dz", &tau0_dz);
	tree->Branch("tau0_decayMode", &tau0_decayMode);
	tree->Branch("tau0_decayModeFindingOldDMs", &tau0_decayModeFindingOldDMs);
	tree->Branch("tau0_decayModeFindingNewDMs", &tau0_decayModeFindingNewDMs);
	tree->Branch("tau0_byCombinedIsolationDeltaBetaCorr3Hits", &tau0_byCombinedIsolationDeltaBetaCorr3Hits);
	tree->Branch("tau0_byLooseCombinedIsolationDeltaBetaCorr3Hits", &tau0_byLooseCombinedIsolationDeltaBetaCorr3Hits);
	tree->Branch("tau0_byMediumCombinedIsolationDeltaBetaCorr3Hits", &tau0_byMediumCombinedIsolationDeltaBetaCorr3Hits);
	tree->Branch("tau0_byTightCombinedIsolationDeltaBetaCorr3Hits", &tau0_byTightCombinedIsolationDeltaBetaCorr3Hits);
	tree->Branch("tau0_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03", &tau0_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03);
	tree->Branch("tau0_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03", &tau0_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03);
	tree->Branch("tau0_byTightCombinedIsolationDeltaBetaCorr3HitsdR03", &tau0_byTightCombinedIsolationDeltaBetaCorr3HitsdR03);
	tree->Branch("tau0_byLooseIsolationMVArun2v1DBdR03oldDMwLT", &tau0_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau0_byMediumIsolationMVArun2v1DBdR03oldDMwLT", &tau0_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau0_byTightIsolationMVArun2v1DBdR03oldDMwLT", &tau0_byTightIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau0_byVTightIsolationMVArun2v1DBdR03oldDMwLT", &tau0_byVTightIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau0_againstMuonLoose3", &tau0_againstMuonLoose3);
	tree->Branch("tau0_againstMuonTight3", &tau0_againstMuonTight3);
	tree->Branch("tau0_againstElectronVLooseMVA6", &tau0_againstElectronVLooseMVA6);
	tree->Branch("tau0_againstElectronLooseMVA6", &tau0_againstElectronLooseMVA6);
	tree->Branch("tau0_againstElectronMediumMVA6", &tau0_againstElectronMediumMVA6);
	tree->Branch("tau0_againstElectronTightMVA6", &tau0_againstElectronTightMVA6);
	tree->Branch("tau1_pt", &tau1_pt);
	tree->Branch("tau1_eta", &tau1_eta);
	tree->Branch("tau1_phi", &tau1_phi);
	tree->Branch("tau1_E", &tau1_E);
	tree->Branch("tau1_charge", &tau1_charge);
	tree->Branch("tau1_dxy", &tau1_dxy);
	tree->Branch("tau1_dz", &tau1_dz);
	tree->Branch("tau1_decayMode", &tau1_decayMode);
	tree->Branch("tau1_decayModeFindingOldDMs", &tau1_decayModeFindingOldDMs);
	tree->Branch("tau1_decayModeFindingNewDMs", &tau1_decayModeFindingNewDMs);
	tree->Branch("tau1_byCombinedIsolationDeltaBetaCorr3Hits", &tau1_byCombinedIsolationDeltaBetaCorr3Hits);
	tree->Branch("tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits", &tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits);
	tree->Branch("tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits", &tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits);
	tree->Branch("tau1_byTightCombinedIsolationDeltaBetaCorr3Hits", &tau1_byTightCombinedIsolationDeltaBetaCorr3Hits);
	tree->Branch("tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03", &tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03);
	tree->Branch("tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03", &tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03);
	tree->Branch("tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03", &tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03);
	tree->Branch("tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT", &tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT", &tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT", &tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT", &tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau1_againstMuonLoose3", &tau1_againstMuonLoose3);
	tree->Branch("tau1_againstMuonTight3", &tau1_againstMuonTight3);
	tree->Branch("tau1_againstElectronVLooseMVA6", &tau1_againstElectronVLooseMVA6);
	tree->Branch("tau1_againstElectronLooseMVA6", &tau1_againstElectronLooseMVA6);
	tree->Branch("tau1_againstElectronMediumMVA6", &tau1_againstElectronMediumMVA6);
	tree->Branch("tau1_againstElectronTightMVA6", &tau1_againstElectronTightMVA6);

	// jets
	tree->Branch("jet0_pt", &jet0_pt);
	tree->Branch("jet0_eta", &jet0_eta);
	tree->Branch("jet0_phi", &jet0_phi);
	tree->Branch("jet0_E", &jet0_E);
	tree->Branch("jet0_CSV", &jet0_CSV);
	tree->Branch("jet1_pt", &jet1_pt);
	tree->Branch("jet1_eta", &jet1_eta);
	tree->Branch("jet1_phi", &jet1_phi);
	tree->Branch("jet1_E", &jet1_E);
	tree->Branch("jet1_CSV", &jet1_CSV);
	tree->Branch("jet2_pt", &jet2_pt);
	tree->Branch("jet2_eta", &jet2_eta);
	tree->Branch("jet2_phi", &jet2_phi);
	tree->Branch("jet2_E", &jet2_E);
	tree->Branch("jet2_CSV", &jet2_CSV);
	tree->Branch("jet3_pt", &jet3_pt);
	tree->Branch("jet3_eta", &jet3_eta);
	tree->Branch("jet3_phi", &jet3_phi);
	tree->Branch("jet3_E", &jet3_E);
	tree->Branch("jet3_CSV", &jet3_CSV);
	// MET
	tree->Branch("PFMET", &PFMET);
	tree->Branch("PFMETphi", &PFMETphi);
	tree->Branch("MHT", &MHT);
	tree->Branch("metLD", &metLD);
}

//ClassImp(CU_ttH_EDA_Ntuple);

#endif
