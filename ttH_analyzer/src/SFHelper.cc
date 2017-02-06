#ifndef SFHelper_cc
#define SFHelper_cc

#include "Analyzers/ttH_analyzer/interface/SFHelper.h"

// constructor
SFHelper::SFHelper(analysis_types analysis, Selection_types selection, bool isdata)
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

void SFHelper::Set_up_TauSF_Lut()
{
	file_fr_tau = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/FR_tau_2016.root").c_str(), "read");
	
	f_fakerate_tau_mvaM_etaL_ratio = (TF1*) file_fr_tau->Get("jetToTauFakeRate/dR03mvaMedium/absEtaLt1_5/fitFunction_data_div_mc_hadTaus_pt");
	f_fakerate_tau_mvaM_etaH_ratio = (TF1*) file_fr_tau->Get("jetToTauFakeRate/dR03mvaMedium/absEta1_5to9_9/fitFunction_data_div_mc_hadTaus_pt");
}

void SFHelper::Set_up_PUWeight_hist()
{
	file_puweight = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/PU_weights/PU_weights_2016_271036_284044.root").c_str(), "read");
	h_puweight = (TH1F*) file_puweight->Get("h_ratio_data_MC");
}

void SFHelper::Set_up_FakeRate_Lut()
{
	// electrons and muons
	file_fr_lep = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/FR_data_ttH_mva.root").c_str(),"read");
	
	h_fakerate_el = (TH2F*) file_fr_lep->Get("FR_mva075_el_data_comb");
	h_fakerate_mu = (TH2F*) file_fr_lep->Get("FR_mva075_mu_data_comb");
	
	//taus
	/*
	  file_fr_tau = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/FR_tau_2016.root").c_str(), "read");
	  
	  g_fakerate_tau_mvaM_etaL_mc = (TGraphAsymmErrors*) file_fr_tau->Get("jetToTauFakeRate/dR03mvaMedium/absEtaLt1_5/jetToTauFakeRate_mc_hadTaus_pt");
	  g_fakerate_tau_mvaM_etaH_mc = (TGraphAsymmErrors*) file_fr_tau->Get("jetToTauFakeRate/dR03mvaMedium/absEta1_5to9_9/jetToTauFakeRate_mc_hadTaus_pt");
	  f_fakerate_tau_mvaM_etaL_ratio = (TF1*) file_fr_tau->Get("jetToTauFakeRate/dR03mvaMedium/absEtaLt1_5/fitFunction_data_div_mc_hadTaus_pt");
	  f_fakerate_tau_mvaM_etaH_ratio = (TF1*) file_fr_tau->Get("jetToTauFakeRate/dR03mvaMedium/absEta1_5to9_9/fitFunction_data_div_mc_hadTaus_pt");
	*/
}

void SFHelper::Set_up_ChargeMisID_Lut()
{
	file_eleMisCharge = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/QF_data_el.root").c_str(), "read");
	h_chargeMisId = (TH2F*) file_eleMisCharge->Get("chargeMisId");
}

void SFHelper::Set_up_LeptonSF_Lut()
{
	//// loose vs reco
	/// for muons
	// loose vs reco
	file_recoToLoose_leptonSF_mu1_b = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/mu_ttH_presel_barrel.root").c_str(),"read");
	file_recoToLoose_leptonSF_mu1_e = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/mu_ttH_presel_endcap.root").c_str(),"read");
	file_recoToLoose_leptonSF_mu2 = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/MuonID_Z_RunBCD_prompt80X_7p65_looseID.root").c_str(),"read");
	file_recoToLoose_leptonSF_mu3 = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/ratios_HIP_trkEff.root").c_str(),"read");
	
	h_recoToLoose_leptonSF_mu1_b = (TGraphAsymmErrors*)(file_recoToLoose_leptonSF_mu1_b->Get("ratio"));
	h_recoToLoose_leptonSF_mu1_e = (TGraphAsymmErrors*)(file_recoToLoose_leptonSF_mu1_e->Get("ratio"));
	h_recoToLoose_leptonSF_mu2 = (TH2F*)(file_recoToLoose_leptonSF_mu2->Get("pt_abseta_ratio_MC_NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1"));
	h_recoToLoose_leptonSF_mu3 = (TGraphAsymmErrors*)(file_recoToLoose_leptonSF_mu3->Get("ratio_eta"));

	/// for electrons
	// loose vs reco
	file_recoToLoose_leptonSF_el = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/el_scaleFactors_20160724.root").c_str(),"read");
	file_recoToLoose_leptonSF_gsf = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/el_scaleFactors_gsf.root").c_str(),"read");

	h_recoToLoose_leptonSF_el1 = (TH2F*)(file_recoToLoose_leptonSF_el->Get("GsfElectronToFOID2D"));
	h_recoToLoose_leptonSF_el2 = (TH2F*)(file_recoToLoose_leptonSF_el->Get("MVAVLooseElectronToMini4"));
	h_recoToLoose_leptonSF_el3 = (TH2F*)(file_recoToLoose_leptonSF_el->Get("MVAVLooseElectronToConvIHit1"));
	h_recoToLoose_leptonSF_gsf = (TH2F*)(file_recoToLoose_leptonSF_gsf->Get("EGamma_SF2D"));

	//// tight vs loose
	if (_analysis == Analyze_2lss1tau) {
		/// for muon
		file_looseToTight_leptonSF_mu_2lss = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/lepMVAEffSF_m_2lss.root").c_str(),"read");
		h_looseToTight_leptonSF_mu_2lss = (TH2F*)(file_looseToTight_leptonSF_mu_2lss->Get("sf"));
		
		/// for electron
		file_looseToTight_leptonSF_el_2lss = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/lepMVAEffSF_e_2lss.root").c_str(),"read");
		h_looseToTight_leptonSF_el_2lss = (TH2F*)(file_looseToTight_leptonSF_el_2lss->Get("sf"));
	}
	
	if (_analysis == Analyze_3l) {
		/// for muon
		file_looseToTight_leptonSF_mu_3l = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/lepMVAEffSF_m_3l.root").c_str(),"read");
		h_looseToTight_leptonSF_mu_3l = (TH2F*)(file_looseToTight_leptonSF_mu_3l->Get("sf"));
		
		/// for electron
		file_looseToTight_leptonSF_el_3l = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/lepMVAEffSF_e_3l.root").c_str(),"read");
		h_looseToTight_leptonSF_el_3l = (TH2F*)(file_looseToTight_leptonSF_el_3l->Get("sf"));
	}
}

void SFHelper::Delete_FakeRate_Lut()
{
	file_fr_lep->Close();		
	delete file_fr_lep;
}

void SFHelper::Delete_TauSF_Lut()
{
	file_fr_tau->Close();		
	delete file_fr_tau;
}

void SFHelper::Delete_ChargeMisID_Lut()
{
	file_eleMisCharge->Close();
	delete file_eleMisCharge;
}

void SFHelper::Delete_PUWeight_hist()
{
	file_puweight->Close();
	delete file_puweight;
}

void SFHelper::Delete_LeptonSF_Lut()
{
	file_recoToLoose_leptonSF_mu1_b->Close();
	file_recoToLoose_leptonSF_mu1_e->Close();
	file_recoToLoose_leptonSF_mu2->Close();
	file_recoToLoose_leptonSF_mu3->Close();
	
	file_recoToLoose_leptonSF_el->Close();
	file_recoToLoose_leptonSF_gsf->Close();
	
	delete file_recoToLoose_leptonSF_mu1_b;
	delete file_recoToLoose_leptonSF_mu1_e;
	delete file_recoToLoose_leptonSF_mu2;
	delete file_recoToLoose_leptonSF_mu3;
	
	delete file_recoToLoose_leptonSF_el;
	delete file_recoToLoose_leptonSF_gsf;

	if (_analysis == Analyze_2lss1tau) {
		file_looseToTight_leptonSF_mu_2lss->Close();
		file_looseToTight_leptonSF_el_2lss->Close();
		
		delete file_looseToTight_leptonSF_mu_2lss;
		delete file_looseToTight_leptonSF_el_2lss;
	}
	if (_analysis == Analyze_3l) {
		file_looseToTight_leptonSF_mu_3l->Close();
		file_looseToTight_leptonSF_el_3l->Close();
		
		delete file_looseToTight_leptonSF_mu_3l;
		delete file_looseToTight_leptonSF_el_3l;
	}
}

void SFHelper::Delete_BTagCalibration_Readers()
{
	delete BTagCaliReader;
}

float SFHelper::Get_HLTSF(int lepCategory)
{
	// lepCategory: 0=mumu; 1=ee; 2=emu
	if (lepCategory == 0)
		return 1.01;
	else if (lepCategory == 1)
		return 1.02;
	else if (lepCategory == 2)
		return 1.02;

	std::cerr << "not valid lepton category !" << std::endl;
	assert(0);
	return 0.;
}

float SFHelper::Get_LeptonIDSF(const miniLepton& lepton)
{
	assert(not _isdata);
	
	float sf = 1.;
	assert(lepton.passLooseSel());

	sf *= Get_LeptonSF_loose(lepton);
	// SF fakeable to loose is currently assumed to be 1.
	if (lepton.passTightSel())
		sf *= Get_LeptonSF_tight_vs_loose(lepton);

	return sf;
}

float SFHelper::Get_LeptonSF_loose(const miniLepton& lepton)
{
	float sf =1.;

	if (lepton.Type()==LeptonType::kmu) {
		if (abs(lepton.eta())<1.2) {
			sf *= evalTGraph(h_recoToLoose_leptonSF_mu1_b, lepton.pt());
		}
		else {
			sf *= evalTGraph(h_recoToLoose_leptonSF_mu1_e, lepton.pt());
		}
		
		sf *= read2DHist(h_recoToLoose_leptonSF_mu2,
						 lepton.pt(), abs(lepton.eta()));
		
		sf *= evalTGraph(h_recoToLoose_leptonSF_mu3, lepton.eta());
	}
	else if (lepton.Type()==LeptonType::kele) {
		sf *= read2DHist(h_recoToLoose_leptonSF_el1,
						 lepton.pt(), abs(lepton.eta()));
		sf *= read2DHist(h_recoToLoose_leptonSF_el2,
						 lepton.pt(), abs(lepton.eta()));
		sf *= read2DHist(h_recoToLoose_leptonSF_el3,
						 lepton.pt(), abs(lepton.eta()));
		// !! different pt eta xaxis
		sf *= read2DHist(h_recoToLoose_leptonSF_gsf, lepton.eta(), lepton.pt());
	}

	return sf;
}

float SFHelper::Get_LeptonSF_tight_vs_loose(const miniLepton& lepton)
{
	float sf = 1.;
	
	assert(lepton.passTightSel());

	if (_analysis == Analyze_2lss1tau) {
		if (lepton.Type()==LeptonType::kmu) {
			sf = read2DHist(h_looseToTight_leptonSF_mu_2lss,
							lepton.pt(), abs(lepton.eta()));
		}
		else if (lepton.Type()==LeptonType::kele) {
			sf = read2DHist(h_looseToTight_leptonSF_el_2lss,
							lepton.pt(), abs(lepton.eta()));
		}
	}
	else if (_analysis == Analyze_3l) {
		if (lepton.Type()==LeptonType::kmu) {
			sf = read2DHist(h_looseToTight_leptonSF_mu_3l,
							lepton.pt(), abs(lepton.eta()));
		}
		else if (lepton.Type()==LeptonType::kele) {
			sf = read2DHist(h_looseToTight_leptonSF_el_3l,
							lepton.pt(), abs(lepton.eta()));
		}
	}

	return sf;
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

float SFHelper::Get_PUWeight(int nPU)
{
	assert(not _isdata);
	
	int xbin = h_puweight->FindBin(nPU);
	return h_puweight->GetBinContent(xbin);
}

float SFHelper::Get_TauIDSF(const pat::Tau& tau, bool isGenMatched)
{
	assert(not _isdata);

	// tau ID efficiency data/MC scale factor
	float tauEff_sf = 0.9;
	// tau ID fake rate data/MC scale factor
	float tauFR_sf = tau.eta() < 1.497 ?
	    readTF(f_fakerate_tau_mvaM_etaL_ratio, tau.pt()) :
		readTF(f_fakerate_tau_mvaM_etaH_ratio, tau.pt());

	return (isGenMatched ? tauEff_sf : tauFR_sf);
}

float SFHelper::Get_MCWeight()
{
	return 0.;
}

float SFHelper::Get_FakeRate(const miniLepton& lepton)
{
	float fakerate = 0;
	
	if (lepton.Type() == LeptonType::kele)
		fakerate = read2DHist(h_fakerate_el, lepton.conePt(), lepton.eta());
	else if (lepton.Type() == LeptonType::kmu)
		fakerate = read2DHist(h_fakerate_mu, lepton.conePt(), lepton.eta());
	
	if (lepton.conePt() < 10.) return 0.;
	
	return fakerate;
}

float SFHelper::Get_FakeRate(const pat::Tau& tau)
{
	float fr_mc = 0;
	float ratio = 0;

	if (abs(tau.eta()) < 1.479) {
		fr_mc = readTGraph(g_fakerate_tau_mvaM_etaL_mc, tau.pt());
		ratio = readTF(f_fakerate_tau_mvaM_etaL_ratio, tau.pt());
	}
	else {
		fr_mc = readTGraph(g_fakerate_tau_mvaM_etaH_mc, tau.pt());
		ratio = readTF(f_fakerate_tau_mvaM_etaH_ratio, tau.pt());
	}
	
	return fr_mc * ratio;
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

	// only apply the charge flip rate to the electron that is same sign as tau
	// due to the tau charge requirement in signal region
	if (lepton.charge() * tauCharge < 0) return 0;

	return read2DHist(h_chargeMisId, lepton.pt(), abs(lepton.eta()));
}

float SFHelper::read2DHist(TH2* h2d, float x, float y)
{
	TAxis* xaxis = h2d->GetXaxis();
	int nbinx = xaxis->GetNbins();
	int xbin = xaxis->FindBin(x);
    if (xbin < 1) xbin = 1;
	if (xbin > nbinx) xbin = nbinx;
	
	TAxis* yaxis = h2d->GetYaxis();
	int nbiny = yaxis->GetNbins();
	int ybin = yaxis->FindBin(abs(y));
    if (ybin < 1) ybin = 1;
	if (ybin > nbiny) ybin = nbiny;
	
	float result = h2d->GetBinContent(xbin, ybin);

	return result;
}

float SFHelper::evalTGraph(TGraphAsymmErrors* graph, float x)
{
	float x1 = std::max(float(graph->GetXaxis()->GetXmin()+1e-5),
						std::min(float(graph->GetXaxis()->GetXmax()-1e-5), x)
						);
	return graph->Eval(x1);
}

float SFHelper::readTGraph(TGraphAsymmErrors* graph, float x)
{
	int nPoints = graph -> GetN();
	double xp, yp = -1.;
	for (int i=0; i<nPoints; ++i) {
		int ip = graph -> GetPoint(i,xp,yp);
		float xerr_h = graph->GetErrorXhigh(i);
		float xerr_l = graph->GetErrorXlow(i);

		if (ip != -1) {
			if (x>xp-xerr_l and x<xp+xerr_h) return yp;
		}
	}
	return yp;
}

float SFHelper::readTF(TF1* f, float x)
{
	//float x1 = std::max(float(f->GetXaxis()->GetXmin()),
	//					std::min(float(f->GetXaxis()->GetXmax()), x)
	//					);
	return f->Eval(x);
}

#endif
