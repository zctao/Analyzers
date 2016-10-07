#ifndef kinMVA_2lss_cc
#define kinMVA_2lss_cc

#include "Analyzers/ttH_analyzer/interface/kinMVA_2lss.h"

float kinMVA_2lss::Set_max_lep_eta(double eta1, double eta2)
{	
	return max( abs(eta1), abs(eta2) );
}

float kinMVA_2lss::Set_MT_met_lep1(double l1_conePt, double l1_phi,
								  const pat::MET& met)
{
    return sqrt( 2*l1_conePt*met.pt()*(1-cos(l1_phi-met.phi())) );	
}

float kinMVA_2lss::mindr_lep_jet(float lep_eta, float lep_phi,
								 const vector<pat::Jet>& jets)
{
	if (jets.size() < 1) return -9999.;  // need at least one jet
	
	float mindr = 666.;

	for (auto & j : jets) {
		float dr = DeltaR(lep_eta, lep_phi, j.eta(), j.phi());
		if ( dr < mindr)
			mindr = dr;		
	}

	return mindr;
}

float kinMVA_2lss::Set_nJet25(const vector<pat::Jet>& jets)
{
	int nJ = 0;
	for (auto & j : jets) {
		if (j.pt() >= 25.)
			++nJ;
	}
	
	return nJ;
}

#endif
