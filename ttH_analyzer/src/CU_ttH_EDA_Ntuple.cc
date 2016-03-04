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

//ClassImp(CU_ttH_EDA_Ntuple);

#endif
