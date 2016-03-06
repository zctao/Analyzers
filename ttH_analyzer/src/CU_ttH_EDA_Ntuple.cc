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
	nEvent = local.event_nr;  // Need to verify
	n_presel_mu = local.n_muons;
	n_presel_ele = local.n_electrons;
	n_presel_tau = local.n_taus;
	n_presel_jet = local.n_jets;
	
	fill_ntuple_electrons(local);
	fill_ntuple_muons(local);
	fill_ntuple_taus(local);
	fill_ntuple_jets(local);
	fill_ntuple_met(local);
	
}

void CU_ttH_EDA_Ntuple::fill_ntuple_muons(const CU_ttH_EDA_event_vars &local)
{
	if (local.mu_selected_sorted.size() > 0) {
		mu0_pt = local.mu_selected_sorted[0].pt();
		mu0_eta = local.mu_selected_sorted[0].eta();
		mu0_phi = local.mu_selected_sorted[0].phi();
		mu0_E = local.mu_selected_sorted[0].energy();
	}
	
	if (local.mu_selected_sorted.size() > 1 ) {
		mu1_pt = local.mu_selected_sorted[1].pt();
		mu1_eta = local.mu_selected_sorted[1].eta();
		mu1_phi = local.mu_selected_sorted[1].phi();
		mu1_E = local.mu_selected_sorted[1].energy();
	}
}

void CU_ttH_EDA_Ntuple::fill_ntuple_electrons(const CU_ttH_EDA_event_vars &local)
{
	if (local.e_selected_sorted.size() > 0) {
		ele0_pt = local.e_selected_sorted[0].pt();
		ele0_eta = local.e_selected_sorted[0].eta();
		ele0_phi = local.e_selected_sorted[0].phi();
		ele0_E = local.e_selected_sorted[0].energy();
	}
	
	if (local.e_selected_sorted.size() > 1 ) {
		ele1_pt = local.e_selected_sorted[1].pt();
		ele1_eta = local.e_selected_sorted[1].eta();
		ele1_phi = local.e_selected_sorted[1].phi();
		ele1_E = local.e_selected_sorted[1].energy();
	}
}

void CU_ttH_EDA_Ntuple::fill_ntuple_taus(const CU_ttH_EDA_event_vars &local)
{
	if (local.tau_selected_sorted.size() > 0) {
		tau0_pt = local.tau_selected_sorted[0].pt();
		tau0_eta = local.tau_selected_sorted[0].eta();
		tau0_phi = local.tau_selected_sorted[0].phi();
		tau0_E = local.tau_selected_sorted[0].energy();
	}
	
	if (local.tau_selected_sorted.size() > 1 ) {
		tau1_pt = local.tau_selected_sorted[1].pt();
		tau1_eta = local.tau_selected_sorted[1].eta();
		tau1_phi = local.tau_selected_sorted[1].phi();
		tau1_E = local.tau_selected_sorted[1].energy();
	}
}

void CU_ttH_EDA_Ntuple::fill_ntuple_jets(const CU_ttH_EDA_event_vars &local)
{
	if (local.jets_selected_sorted.size() > 0) {
		jet0_pt = local.jets_selected_sorted[0].pt();
		jet0_eta = local.jets_selected_sorted[0].eta();
		jet0_phi = local.jets_selected_sorted[0].phi();
		jet0_E = local.jets_selected_sorted[0].energy();
	}
	
	if (local.jets_selected_sorted.size() > 1 ) {
		jet1_pt = local.jets_selected_sorted[1].pt();
		jet1_eta = local.jets_selected_sorted[1].eta();
		jet1_phi = local.jets_selected_sorted[1].phi();
		jet1_E = local.jets_selected_sorted[1].energy();
	}
	
	if (local.jets_selected_sorted.size() > 2) {
		jet2_pt = local.jets_selected_sorted[2].pt();
		jet2_eta = local.jets_selected_sorted[2].eta();
		jet2_phi = local.jets_selected_sorted[2].phi();
		jet2_E = local.jets_selected_sorted[2].energy();
	}
	
	if (local.jets_selected_sorted.size() > 3 ) {
		jet3_pt = local.jets_selected_sorted[3].pt();
		jet3_eta = local.jets_selected_sorted[3].eta();
		jet3_phi = local.jets_selected_sorted[3].phi();
		jet3_E = local.jets_selected_sorted[3].energy();
	}
}

void CU_ttH_EDA_Ntuple::fill_ntuple_met(const CU_ttH_EDA_event_vars &local)
{
	double metx = local.MET_corrected.px();
	double mety = local.MET_corrected.py();
	PFMET = sqrt(metx*metx+mety*mety);
	}

void CU_ttH_EDA_Ntuple::initialize()
{
	// event variables
	nEvent = -9999;
	n_presel_mu = -9999;
	n_presel_ele = -9999;
	n_presel_tau = -9999;
	n_presel_jet = -9999;

	// muons
	mu0_pt = -9999;
	mu0_eta = -9999;
	mu0_phi = -9999;
	mu0_E = -9999;
	mu0_charge = -9999;
	mu0_jetNDauChargedMVASel = -9999;
	mu0_miniRelIso = -9999;
	mu0_miniIsoCharged = -9999;
	mu0_miniIsoNeutral = -9999;
	mu0_jetPtRel = -9999;
	mu0_jetPtRatio = -9999;
	mu0_jetCSV = -9999;
	mu0_sip3D = -9999;
	mu0_dxy = -9999;
	mu0_dz = -9999;
	mu0_segmentCompatibility = -9999;
	mu0_leptonMVA = -9999;
	mu1_pt = -9999;
	mu1_eta = -9999;
	mu1_phi = -9999;
	mu1_E = -9999;
	mu1_charge = -9999;
	mu1_jetNDauChargedMVASel = -9999;
	mu1_miniRelIso = -9999;
	mu1_miniIsoCharged = -9999;
	mu1_miniIsoNeutral = -9999;
	mu1_jetPtRel = -9999;
	mu1_jetPtRatio = -9999;
	mu1_jetCSV = -9999;
	mu1_sip3D = -9999;
	mu1_dxy = -9999;
	mu1_dz = -9999;
	mu1_segmentCompatibility = -9999;
	mu1_leptonMVA = -9999;
	
	// electrons
	ele0_pt = -9999;
	ele0_eta = -9999;
	ele0_phi = -9999;
	ele0_E = -9999;
	ele0_charge = -9999;
	ele0_jetNDauChargedMVASel = -9999;
	ele0_miniRelIso = -9999;
	ele0_miniIsoCharged = -9999;
	ele0_miniIsoNeutral = -9999;
	ele0_jetPtRel = -9999;
	ele0_jetPtRatio = -9999;
	ele0_jetCSV = -9999;
	ele0_sip3D = -9999;
	ele0_dxy = -9999;
	ele0_dz = -9999;
	ele0_ntMVAeleID = -9999;
	ele0_leptonMVA = -9999;
	ele1_pt = -9999;
	ele1_eta = -9999;
	ele1_phi = -9999;
	ele1_E = -9999;
	ele1_charge = -9999;
	ele1_jetNDauChargedMVASel = -9999;
	ele1_miniRelIso = -9999;
	ele1_miniIsoCharged = -9999;
	ele1_miniIsoNeutral = -9999;
	ele1_jetPtRel = -9999;
	ele1_jetPtRatio = -9999;
	ele1_jetCSV = -9999;
	ele1_sip3D = -9999;
	ele1_dxy = -9999;
	ele1_dz = -9999;
	ele1_ntMVAeleID = -9999;
	ele1_leptonMVA = -9999;	
	
	// taus
	tau0_pt = -9999;
	tau0_eta = -9999;
	tau0_phi = -9999;
	tau0_E = -9999;
	tau0_charge = -9999;
	tau0_dxy = -9999;
	tau0_dz = -9999;
	tau0_decayModeFindingOldDMs = -9999;
	tau0_decayModeFindingNewDMs = -9999;
	tau0_byCombinedIsolationDeltaBetaCorr3Hits = -9999;
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
	tau0_againstElectronVTightMVA6 = -9999;
	tau1_pt = -9999;
	tau1_eta = -9999;
	tau1_phi = -9999;
	tau1_E = -9999;
	tau1_charge = -9999;
	tau1_dxy = -9999;
	tau1_dz = -9999;
	tau1_decayModeFindingOldDMs = -9999;
	tau1_decayModeFindingNewDMs = -9999;
	tau1_byCombinedIsolationDeltaBetaCorr3Hits = -9999;
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
	tau1_againstElectronVTightMVA6 = -9999;
	
	// jets
	jet0_pt = -9999;
	jet0_eta = -9999;
	jet0_phi = -9999;
	jet0_E = -9999;
	jet0_CSV = -9999;
	jet1_pt = -9999;
	jet1_eta = -9999;
	jet1_phi = -9999;
	jet1_E = -9999;
	jet1_CSV = -9999;
	jet2_pt = -9999;
	jet2_eta = -9999;
	jet2_phi = -9999;
	jet2_E = -9999;
	jet2_CSV = -9999;
	jet3_pt = -9999;
	jet3_eta = -9999;
	jet3_phi = -9999;
	jet3_E = -9999;
	jet3_CSV = -9999;
	// MET
	PFMET = -9999;
	PFMETphi = -9999;
}

void CU_ttH_EDA_Ntuple::set_up_branches(TTree *tree)
{
	// initialize ntuple
	this->initialize();

	// set up tree branches
	tree->Branch("nEvent", &nEvent);
	tree->Branch("n_presel_mu", &n_presel_mu);
	tree->Branch("n_presel_ele", &n_presel_ele);
	tree->Branch("n_presel_tau", &n_presel_tau);
	tree->Branch("n_presel_jet", &n_presel_jet);
	// muons
	tree->Branch("mu0_pt",                   &mu0_pt);
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
	tree->Branch("mu1_pt",                   &mu1_pt);
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
	// electrons
	tree->Branch("ele0_pt",                   &ele0_pt);
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
	tree->Branch("ele1_pt",                   &ele1_pt);
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
	// taus
	tree->Branch("tau0_pt", &tau0_pt);
	tree->Branch("tau0_eta", &tau0_eta);
	tree->Branch("tau0_phi", &tau0_phi);
	tree->Branch("tau0_E", &tau0_E);
	tree->Branch("tau0_charge", &tau0_charge);
	tree->Branch("tau0_dxy", &tau0_dxy);
	tree->Branch("tau0_dz", &tau0_dz);
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
	tree->Branch("tau0_againstElectronVTightMVA6", &tau0_againstElectronVTightMVA6);
	tree->Branch("tau1_pt", &tau1_pt);
	tree->Branch("tau1_eta", &tau1_eta);
	tree->Branch("tau1_phi", &tau1_phi);
	tree->Branch("tau1_E", &tau1_E);
	tree->Branch("tau1_charge", &tau1_charge);
	tree->Branch("tau1_dxy", &tau1_dxy);
	tree->Branch("tau1_dz", &tau1_dz);
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
	tree->Branch("tau1_againstElectronVTightMVA6", &tau1_againstElectronVTightMVA6);
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
	// MET
	tree->Branch("PFMET", &PFMET);
	tree->Branch("PFMETphi", &PFMETphi);
}

//ClassImp(CU_ttH_EDA_Ntuple);

#endif
