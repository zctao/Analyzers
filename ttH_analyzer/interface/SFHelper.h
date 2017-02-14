#ifndef SFHelper_h
#define SFHelper_h

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

/// BTag Calibration
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "Analyzers/ttH_analyzer/interface/NtupleSFHelper.h"
#include "Analyzers/ttH_analyzer/interface/miniLepton.h"

class SFHelper : public NtupleSFHelper
{
 public:
	// constructor and destructor
	SFHelper(Analysis_types, Selection_types, bool);
	~SFHelper();

	// member functions
	// overload base class functions
	//using NtupleSFHelper::Get_LeptonIDSF;
	using NtupleSFHelper::Get_TauIDSF;
	using NtupleSFHelper::Get_FakeRate;
	using NtupleSFHelper::Get_EleChargeMisIDProb;
	using NtupleSFHelper::Get_LeptonSF_loose;
	using NtupleSFHelper::Get_LeptonSF_tight_vs_loose;
		
	float Get_LeptonIDSF(const miniLepton&);	
	float Get_EvtCSVWeight(std::vector<pat::Jet> &, const std::string &);
	float Get_TauIDSF(const pat::Tau&, bool);
	float Get_FakeRate(const miniLepton&);
	float Get_FakeRate(const pat::Tau&);
	float Get_EleChargeMisIDProb(const miniLepton&, int);
	//float Get_MCWeight();
	//float Get_ChargeFlipWeight();
	
 private:

	bool _isdata;
	Analysis_types _analysis;
	Selection_types _selection;

	// CSV
	BTagCalibrationReader* BTagCaliReader;
	//std::string Btag_sysList[16] =
	//	{"LFUp","LFDown","HFUp","HFDown",
	//	 "HFStats1Up","HFStats1Down","HFStats2Up","HFStats2Down",
	//	 "LFStats1Up","LFStats1Down","LFStats2Up","LFStats2Down",
	//	 "cErr1Up","cErr1Down","cErr2Up","cErr2Down"};

	void Set_up_BTagCalibration_Readers();
	
	void Delete_BTagCalibration_Readers();

	float Get_LeptonSF_loose(const miniLepton&);
	float Get_LeptonSF_tight_vs_loose(const miniLepton&);
	float Get_JetCSVWeight(pat::Jet&, std::string);

};

#endif
