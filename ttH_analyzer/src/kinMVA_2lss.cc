#ifndef kinMVA_2lss_cc
#define kinMVA_2lss_cc

#include "Analyzers/ttH_analyzer/interface/kinMVA_2lss.h"

void kinMVA_2lss::Set_max_lep_eta(double eta1, double eta2)
{	
	max_lep_eta = max( abs(eta1), abs(eta2) );
}

void kinMVA_2lss::Set_MT_met_lep1(double l1_conePt, double l1_phi,
								  const pat::MET& met)
{
	MT_met_lep1 = sqrt( 2*l1_conePt*met.pt()*(1-cos(l1_phi-met.phi())) );	
}

void kinMVA_2lss::Set_mindr_lep_jet(float& mindr_LJ, float lep_eta,
									float lep_phi, const vector<pat::Jet>& jets)
{
	assert(jets.size()>=1);  // at least one jet
	
	double mindr = 666.;

	for (auto & j : jets) {
		double dr = DeltaR(lep_eta, lep_phi, j.eta(), j.phi());
		if ( dr < mindr)
			mindr = dr;		
	}

	mindr_LJ = mindr;
}

void kinMVA_2lss::Set_nJet25(const vector<pat::Jet>& jets)
{
	int nJ = 0;
	for (auto & j : jets) {
		if (j.pt() >= 25.)
			++nJ;
	}
	
	nJet25 = nJ;
}

/*
template <typename T>
float kinMVA_2lss::Compute_lep_conePt(const T& lep)
{
	// All input leptons should pass fakeable selection
	bool is_tight = lep.userFloat("idMVABased") > 0.5;
	
	if ( is_tight )
		return lep.pt();
	else
		return 0.85 * lep.pt() / lep.userFloat("nearestJetPtRatio");	
}

template float kinMVA_2lss::Compute_lep_conePt<pat::Muon>(const pat::Muon&, const vector<pat::Muon>&);
template float kinMVA_2lss::Compute_lep_conePt<pat::Electron>(const pat::Electron&, const vector<pat::Electron>&);
*/
float kinMVA_2lss::Compute_lep_conePt(const pat::Muon& mu)
{
	// All input leptons should pass fakeable selection
	bool is_tight = mu.userFloat("idMVABased") > 0.5;
	// extra selection for tight ? (CMSSW_7_6_X, AN321_v4)
	bool extra = false;
	if (mu.innerTrack().isAvailable()) {
		if (mu.innerTrack()->ptError()/mu.innerTrack()->pt() < 0.2)
			extra = true;
	}

	is_tight = is_tight and extra;


	
	if ( is_tight)// or mu.userFloat("nearestJetPtRatio")<0.)
		return mu.pt();
	else {
		assert(mu.userFloat("nearestJetPtRatio")>0.);
		// lepton could be categorized as 'fakeable' but no ptRatio calculated (default value -1.) if e.g. dR(lep, jet)>0.5   NEED TO FOLLOW UP
		return 0.85 * mu.pt() / mu.userFloat("nearestJetPtRatio");
	}
	
}

float kinMVA_2lss::Compute_lep_conePt(const pat::Electron& ele)
{
	// All input leptons should pass fakeable selection
	bool is_tight = ele.userFloat("idMVABased") > 0.5;
	// extra selection for tight ? (CMSSW_7_6_X, AN321_v4)
	bool extra = false;
	if ( ele.userFloat("numMissingHits") == 0 and ele.passConversionVeto() and
		ele.isGsfCtfScPixChargeConsistent() ) {
			extra = true;
	}

	is_tight = is_tight and extra;

	if ( is_tight)// or ele.userFloat("nearestJetPtRatio")<0.)
		return ele.pt();
	else {
		assert(ele.userFloat("nearestJetPtRatio")>0.);
		// lepton could be categorized as 'fakeable' but no ptRatio calculated (default value -1.) if e.g. dR(lep, jet)>0.5   NEED TO FOLLOW UP
		return 0.85 * ele.pt() / ele.userFloat("nearestJetPtRatio");
	}
	
}

void kinMVA_2lss::assign_lep_kinVars(const CU_ttH_EDA_event_vars& event,
						float& lep1_conePt, float& lep2_conePt,
						float& lep1_eta, float& lep2_eta,
						float& lep1_phi, float& lep2_phi)
{
	assert(event.mu_fakeable.size()+event.e_fakeable.size() >= 2);
	
	lep1_conePt = 0.;
	lep2_conePt = 0.;

	for (auto & mu : event.mu_fakeable) {
		//float conePt = Compute_lep_conePt<pat::Muon>(mu);
		float conePt = Compute_lep_conePt(mu);
		assert(conePt>0);
		if (conePt > lep1_conePt) {
			lep2_conePt = lep1_conePt;
			lep2_eta = lep1_eta;
			lep2_phi = lep1_phi;
			lep1_conePt = conePt;
			lep1_eta = mu.eta();
			lep1_phi = mu.phi();
		}
		else if (conePt > lep2_conePt) {
			lep2_conePt = conePt;
			lep2_eta = mu.eta();
			lep2_phi = mu.phi();
		}
	}

	for (auto & ele : event.e_fakeable) {
		//float conePt = Compute_lep_conePt<pat::Electron>(ele);
		float conePt = Compute_lep_conePt(ele);
		assert(conePt>0);
		if (conePt > lep1_conePt) {
			lep2_conePt = lep1_conePt;
			lep2_eta = lep1_eta;
			lep2_phi = lep1_phi;
			lep1_conePt = conePt;
			lep1_eta = ele.eta();
			lep1_phi = ele.phi();
		}
		else if (conePt > lep2_conePt) {
			lep2_conePt = conePt;
			lep2_eta = ele.eta();
			lep2_phi = ele.phi();
		}
	}

	assert(lep1_conePt >= lep2_conePt and lep2_conePt > 0.);
}

#endif
