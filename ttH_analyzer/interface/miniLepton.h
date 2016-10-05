#ifndef miniLepton_h
#define miniLepton_h

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TLorentzVector.h"

//enum LeptonID{Loose, Fakeable, Tight};
enum LeptonType{kele, kmu};

class miniLepton
{
 public:

	// constructor and destructor
	miniLepton(const pat::Electron &);
	miniLepton(const pat::Muon &);
	~miniLepton(){};

	// member functions
	void SetLeptonID();
	void Set_pt(float ipt) {_pt = ipt;}
	void Set_eta(float ieta) {_eta = ieta;}
	void Set_phi(float iphi) {_phi = iphi;}
	void Set_mass(float imass) {_mass = imass;}
	void Set_charge(int icharge) {_charge = icharge;}

	float pt() const {return _pt;}
	float eta() const {return _eta;}	
	float phi() const {return _phi;}
	float mass() const {return _mass;}
	int charge() const {return _charge;}
	LeptonType Type() const {return _type;}
	TLorentzVector p4() const;
	
	//LeptonID LepID();
	bool passLooseSel() const {assert(_lepID_set); return _passloose;}
	bool passFakeableSel() const {assert(_lepID_set); return _passfakeable;}
	bool passTightSel() const {assert(_lepID_set); return _passtight;}

	bool tightCharge() const {return _tight_charge;}
	bool conversionVeto() const {assert(_type==LeptonType::kele); return _conversion_veto;}
	bool noMissingHits() const {assert(_type==LeptonType::kele); return _no_missinghits;}
	
	float conePt() const;
	
 private:
	
	float _pt;
	float _conept;
	float _eta;
	float _phi;
	float _mass;
	int _charge;
	
	bool _passloose;
	bool _passfakeable;
	bool _passtight;

	bool _tight_charge;
	bool _conversion_veto;
	bool _no_missinghits;
	
	LeptonType _type;

	bool _lepID_set;

	void Set_LeptonID(const pat::Electron &);
	void Set_LeptonID(const pat::Muon &);
	void Set_tightCharge(const pat::Electron &);
	void Set_tightCharge(const pat::Muon &);
	void Set_conePt();

};

#endif
