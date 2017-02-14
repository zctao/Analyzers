#ifndef SFHelper_cc
#define SFHelper_cc

#include "Analyzers/ttH_analyzer/interface/SFHelper.h"

// constructor
SFHelper::SFHelper(Analysis_types analysis, Selection_types selection, bool isdata)
	: NtupleSFHelper(analysis,selection,isdata)
{
	_analysis = analysis;
	_selection = selection;
	_isdata = isdata;

	if (not _isdata) {
		Set_up_TauSF_Lut();
		Set_up_PUWeight_hist();
		Set_up_LeptonSF_Lut();
		Set_up_BTagCalibration_Readers();		
	}

	if (_selection == Control_1lfakeable) {
		Set_up_FakeRate_Lut();
	}

	if (_selection == Control_2los1tau) {
		Set_up_ChargeMisID_Lut();
	}
}

SFHelper::~SFHelper()
{
	if (not _isdata) {
		Delete_TauSF_Lut();
		Delete_PUWeight_hist();
		Delete_LeptonSF_Lut();
		Delete_BTagCalibration_Readers();
	}

	if (_selection == Control_1lfakeable) {
		Delete_FakeRate_Lut();
	}

	if (_selection == Control_2los1tau) {
		Delete_ChargeMisID_Lut();
	}
}

void SFHelper::Set_up_BTagCalibration_Readers()
{
	const std::string base =
		std::string(getenv("CMSSW_BASE")) +  "/src/Analyzers/ttH_analyzer/data/";
	
	BTagCalibration calib_csvv2("csvv2", base + "CSVv2Moriond17_2017_1_26_BtoH.csv");
	
	BTagCaliReader = new BTagCalibrationReader(
	     BTagEntry::OP_RESHAPING, // operating point
		 "central",
		 {"up_jes", "down_jes", "up_lf", "down_lf", "up_hf", "down_hf",
				 "up_hfstats1", "down_hfstats1", "up_hfstats2", "down_hfstats2",
				 "up_lfstats1", "down_lfstats1", "up_lfstats2", "down_lfstats2",
				 "up_cferr1", "down_cferr1", "up_cferr2", "down_cferr2"}
											   );
	
	BTagCaliReader->load(calib_csvv2,BTagEntry::FLAV_B,"iterativefit");
	BTagCaliReader->load(calib_csvv2,BTagEntry::FLAV_C,"iterativefit");
	BTagCaliReader->load(calib_csvv2,BTagEntry::FLAV_UDSG,"iterativefit");
}



void SFHelper::Delete_BTagCalibration_Readers()
{
	delete BTagCaliReader;
}

float SFHelper::Get_LeptonSF_loose(const miniLepton& lepton)
{
	return Get_LeptonSF_loose(lepton.pt(), lepton.eta(),
							  lepton.Type()==LeptonType::kele,
							  lepton.Type()==LeptonType::kmu);
}


float SFHelper::Get_LeptonSF_tight_vs_loose(const miniLepton& lepton)
{
	assert(lepton.passTightSel());

	return Get_LeptonSF_tight_vs_loose(lepton.pt(), lepton.eta(),
									   lepton.Type()==LeptonType::kele,
									   lepton.Type()==LeptonType::kmu);
}

float SFHelper::Get_EvtCSVWeight(std::vector<pat::Jet>& jets,
								 const std::string& sys)
{
	assert(not _isdata);
	
	double weight_evt = 1.;
	
	for (auto & j : jets) {		
		double w = Get_JetCSVWeight(j, sys);
		weight_evt *= w;
	}
	
	return weight_evt;
}

float SFHelper::Get_JetCSVWeight(pat::Jet& jet, std::string sys/*pass by copy*/)
{
	double pt = jet.pt();
	double eta = jet.eta();
	double csv = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	int flavor = jet.hadronFlavour();
	
	BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;	
	if ( abs(flavor) == 5 )
		jf = BTagEntry::FLAV_B;
	else if ( abs(flavor) == 4 )
		jf = BTagEntry::FLAV_C;

	double weight_jet = 1.;

	if (sys == "JESUp")
		sys = "up_jes";
	else if (sys == "JESDown")
		sys = "down_jes";
	else if (sys == "LFUp" and jf == BTagEntry::FLAV_B)
		sys = "up_lf";
	else if (sys == "LFDown" and jf == BTagEntry::FLAV_B)
		sys = "down_lf";
	else if (sys == "HFStats1Up" and jf == BTagEntry::FLAV_B)
		sys = "up_hfstats1";
	else if (sys == "HFStats1Down" and jf == BTagEntry::FLAV_B)
		sys = "down_hfstats1";
	else if (sys == "HFStats2Up" and jf == BTagEntry::FLAV_B)
		sys = "up_hfstats2";
	else if (sys == "HFStats2Down" and jf == BTagEntry::FLAV_B)
		sys = "down_hfstats2";
	else if (sys == "HFUp" and jf == BTagEntry::FLAV_UDSG)
		sys = "up_hf";
	else if (sys == "HFDown" and jf == BTagEntry::FLAV_UDSG)
		sys = "down_hf";
	else if (sys == "LFStats1Up" and jf == BTagEntry::FLAV_UDSG)
		sys = "up_lfstats1";
	else if (sys == "LFStats1Down" and jf == BTagEntry::FLAV_UDSG)
		sys = "down_lfstats1";
	else if (sys == "LFStats2Up" and jf == BTagEntry::FLAV_UDSG)
		sys = "up_lfstats2";
	else if (sys == "LFStats2Down" and jf == BTagEntry::FLAV_UDSG)
		sys = "down_lfstats2";
	else if (sys == "cErr1Up" and jf == BTagEntry::FLAV_C)
		sys = "up_cferr1";
	else if (sys == "cErr1Down" and jf == BTagEntry::FLAV_C)
		sys = "down_cferr1";
	else if (sys == "cErr2Up" and jf == BTagEntry::FLAV_C)
		sys = "up_cferr2";
	else if (sys == "cErr2Down" and jf == BTagEntry::FLAV_C)
		sys = "down_cferr2";
	else
		sys = "central";

	weight_jet = BTagCaliReader->eval_auto_bounds(sys, jf, abs(eta), pt, csv);
	
	if (sys == "central") assert(weight_jet > 0.);

	return weight_jet;
}

float SFHelper::Get_TauIDSF(const pat::Tau& tau, bool isGenMatched)
{
	return Get_TauIDSF(tau.pt(), tau.eta(), isGenMatched);
}

float SFHelper::Get_LeptonIDSF(const miniLepton& lepton)
{
	assert(not _isdata);
	
	assert(lepton.passLooseSel());

	float sf = 1.;
	sf *= Get_LeptonSF_loose(lepton);
	// SF fakeable to loose is currently assumed to be 1.
	if (lepton.passTightSel())
		sf *= Get_LeptonSF_tight_vs_loose(lepton);

	return sf;
}

/*float SFHelper::Get_MCWeight()
{
	return 0.;
	}*/

float SFHelper::Get_FakeRate(const miniLepton& lepton)
{
	bool isEle = lepton.Type() == LeptonType::kele;
	bool isMu = lepton.Type() == LeptonType::kmu;

	return Get_FakeRate(lepton.conePt(), lepton.eta(), isEle, isMu);
}

float SFHelper::Get_FakeRate(const pat::Tau& tau)
{
	return Get_FakeRate(tau.pt(),tau.eta());
}

/*float SFHelper::Get_ChargeFlipWeight()
{

}*/

float SFHelper::Get_EleChargeMisIDProb(const miniLepton& lepton, int tauCharge)
{
	// muon
	if (lepton.Type() == LeptonType::kmu) return 0.;
	
	// electron
	assert(lepton.Type() == LeptonType::kele);

	return Get_EleChargeMisIDProb(lepton.pt(), lepton.eta(),
								  lepton.charge(), tauCharge);
}

#endif
