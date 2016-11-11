#ifndef miniLepton_cc
#define miniLepton_cc

#include "Analyzers/ttH_analyzer/interface/miniLepton.h"

miniLepton::miniLepton(const pat::Electron & ele)
{
	_pt = ele.pt();
	_eta = ele.eta();
	_phi = ele.phi();
	_mass = ele.mass();
	_charge = ele.charge();

	_type = LeptonType::kele;
	_lepID_set = false;
	_tight_charge = false;

	Set_LeptonID(ele);
	Set_tightCharge(ele);

	_no_missinghits = ele.userFloat("numMissingHits") == 0;
	_conversion_veto = ele.passConversionVeto();

	MCMatchType = 0;
}

miniLepton::miniLepton(const pat::Muon & mu)
{
	_pt = mu.pt();
	_eta = mu.eta();
	_phi = mu.phi();
	_mass = mu.mass();
	_charge = mu.charge();

	_type = LeptonType::kmu;
	_lepID_set = false;
	_tight_charge = false;

	Set_LeptonID(mu);
	Set_tightCharge(mu);

	//_conversion_veto = false;
	//_no_missinghits = false;

	MCMatchType = 0;
}

void miniLepton::Set_LeptonID(const pat::Electron & ele)
{

	_passloose = ele.userFloat("idPreselection") > 0.5;

	_passfakeable =
		ele.userFloat("idFakeable") > 0.5 and
		ele.userFloat("numMissingHits") == 0 and
		ele.passConversionVeto() and
		_passloose;

	_passtight =
		ele.userFloat("idMVABased") > 0.5 and
		//ele.userFloat("numMissingHits") == 0 and
		//ele.passConversionVeto() and
		//ele.isGsfCtfScPixChargeConsistent() and
		_passfakeable;
	
	_conept = -99;
	if (_passfakeable)
		_conept = (ele.userFloat("leptonMVA") > 0.75) ?
			ele.pt() : 0.85 * ele.pt() / ele.userFloat("nearestJetPtRatio");
	
	_lepID_set = true;
}

void miniLepton::Set_LeptonID(const pat::Muon & mu)
{
	_passloose = mu.userFloat("idPreselection") > 0.5;
	_passfakeable = mu.userFloat("idFakeable") > 0.5 and _passloose;
	_passtight = mu.userFloat("idMVABased") > 0.5 and _passfakeable;
   
	_conept = -99.;
	if (_passfakeable)
		_conept = (mu.userFloat("leptonMVA") > 0.75) ?
			mu.pt() : (0.85 * mu.pt() / mu.userFloat("nearestJetPtRatio"));
	
	_lepID_set = true;
}

void miniLepton::Set_tightCharge(const pat::Electron & e)
{
	_tight_charge =
		e.isGsfCtfScPixChargeConsistent() + e.isGsfScPixChargeConsistent() > 1;
}

void miniLepton::Set_tightCharge(const pat::Muon & mu)
{
	if (mu.innerTrack().isAvailable())
		_tight_charge = mu.innerTrack()->ptError()/mu.innerTrack()->pt() < 0.2;
}

float miniLepton::conePt() const
{
	assert(_lepID_set);
	return _conept;
}

TLorentzVector miniLepton::p4() const
{
	TLorentzVector l;
	l.SetPtEtaPhiM(_pt,_eta,_phi,_mass);
	return l;
}

#endif
