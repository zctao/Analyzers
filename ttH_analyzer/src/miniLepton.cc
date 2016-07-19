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

	Set_LeptonID(ele);
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

	Set_LeptonID(mu);
}

void miniLepton::Set_LeptonID(const pat::Electron & ele)
{

	_passloose = ele.userFloat("idPreselection") > 0.5;

	_passfakeable =
		ele.userFloat("idFakeable") > 0.5 and
		ele.userFloat("numMissingHits") == 0 and
		_passloose;

	_passtight =
		ele.userFloat("idMVABased") > 0.5 and
		ele.userFloat("numMissingHits") == 0 and
		ele.passConversionVeto() and
		ele.isGsfCtfScPixChargeConsistent() and
		_passfakeable;
	
    
	_conept = -99;
	if (_passfakeable and !_passtight)
		_conept = 0.85 * ele.pt() / ele.userFloat("nearestJetPtRatio");
	if (_passtight)
		_conept = ele.pt();
	 
	_lepID_set = true;
}

void miniLepton::Set_LeptonID(const pat::Muon & mu)
{
	_passloose = mu.userFloat("idPreselection") > 0.5;
	_passfakeable = mu.userFloat("idFakeable") > 0.5 and _passloose;

	_passtight = false;
	if (mu.userFloat("idMVABased") > 0.5 and
		mu.innerTrack().isAvailable()) {
		if (mu.innerTrack()->ptError()/mu.innerTrack()->pt() < 0.2)
			_passtight = true;
	}
	_passtight = _passtight and _passfakeable;
   
	_conept = -99;
	if (_passfakeable and !_passtight)
		_conept = 0.85 * mu.pt() / mu.userFloat("nearestJetPtRatio");
	if (_passtight)
		_conept = mu.pt();
	 
	_lepID_set = true;
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
