#ifndef NtupleSFHelper_cc
#define NtupleSFHelper_cc

#include "Analyzers/ttH_analyzer/interface/NtupleSFHelper.h"

NtupleSFHelper::NtupleSFHelper(Analysis_types analysis, Selection_types selection, bool isdata)
{
	_analysis = analysis;
	_selection = selection;
	_isdata = isdata;
	
	if (not _isdata) {
		Set_up_TauSF_Lut();
		Set_up_PUWeight_hist();
		Set_up_LeptonSF_Lut();
		//Set_up_BTagCalibration_Readers();		
	}

	if (_selection == Control_1lfakeable) {
		Set_up_FakeRate_Lut();
	}

	if (_selection == Control_2los1tau) {
		Set_up_ChargeMisID_Lut();
	}
}

NtupleSFHelper::~NtupleSFHelper()
{
	if (not _isdata) {
		Delete_TauSF_Lut();
		Delete_PUWeight_hist();
		Delete_LeptonSF_Lut();
		//Delete_BTagCalibration_Readers();
	}

	if (_selection == Control_1lfakeable) {
		Delete_FakeRate_Lut();
	}
	if (_selection == Control_2los1tau) {
		Delete_ChargeMisID_Lut();
	}
}

void NtupleSFHelper::Set_up_TauSF_Lut()
{
	file_fr_tau = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/FR_tau_2016.root").c_str(), "read");
	
	f_fakerate_tau_mvaM_etaL_ratio = (TF1*) file_fr_tau->Get("jetToTauFakeRate/dR03mvaMedium/absEtaLt1_5/fitFunction_data_div_mc_hadTaus_pt");
	f_fakerate_tau_mvaM_etaH_ratio = (TF1*) file_fr_tau->Get("jetToTauFakeRate/dR03mvaMedium/absEta1_5to9_9/fitFunction_data_div_mc_hadTaus_pt");
}

void NtupleSFHelper::Set_up_PUWeight_hist()
{
	file_puweight = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/PU_weights/PU_weights_2016_ReReco_MCMoriond17_271036_284044.root").c_str(), "read");
	h_puweight = (TH1F*) file_puweight->Get("h_ratio_data_MC");
}

void NtupleSFHelper::Set_up_FakeRate_Lut()
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

void NtupleSFHelper::Set_up_ChargeMisID_Lut()
{
	file_eleMisCharge = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/QF_data_el.root").c_str(), "read");
	h_chargeMisId = (TH2F*) file_eleMisCharge->Get("chargeMisId");
}

void NtupleSFHelper::Set_up_LeptonSF_Lut()
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
	/// for muon
	if (_analysis == Analyze_2lss1tau) {
		file_looseToTight_leptonSF_mu_2lss = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/lepMVAEffSF_m_2lss.root").c_str(),"read");
		h_looseToTight_leptonSF_mu_2lss = (TH2F*)(file_looseToTight_leptonSF_mu_2lss->Get("sf"));
	
		/// for electron
		file_looseToTight_leptonSF_el_2lss = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/lepMVAEffSF_e_2lss.root").c_str(),"read");
		h_looseToTight_leptonSF_el_2lss = (TH2F*)(file_looseToTight_leptonSF_el_2lss->Get("sf"));
	}
	else if (_analysis == Analyze_3l) {
		/// for muon
		file_looseToTight_leptonSF_mu_3l = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/lepMVAEffSF_m_3l.root").c_str(),"read");
		h_looseToTight_leptonSF_mu_3l = (TH2F*)(file_looseToTight_leptonSF_mu_3l->Get("sf"));
		
		/// for electron
		file_looseToTight_leptonSF_el_3l = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/lepMVAEffSF_e_3l.root").c_str(),"read");
		h_looseToTight_leptonSF_el_3l = (TH2F*)(file_looseToTight_leptonSF_el_3l->Get("sf"));
	}
}

void NtupleSFHelper::Delete_FakeRate_Lut()
{
	file_fr_lep->Close();		
	delete file_fr_lep;
}

void NtupleSFHelper::Delete_TauSF_Lut()
{
	file_fr_tau->Close();		
	delete file_fr_tau;
}

void NtupleSFHelper::Delete_ChargeMisID_Lut()
{
	file_eleMisCharge->Close();
	delete file_eleMisCharge;
}

void NtupleSFHelper::Delete_PUWeight_hist()
{
	file_puweight->Close();
	delete file_puweight;
}

void NtupleSFHelper::Delete_LeptonSF_Lut()
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
	else if (_analysis == Analyze_3l) {
		file_looseToTight_leptonSF_mu_3l->Close();
		file_looseToTight_leptonSF_el_3l->Close();
		
		delete file_looseToTight_leptonSF_mu_3l;
		delete file_looseToTight_leptonSF_el_3l;
	}
}

float NtupleSFHelper::Get_HLTSF(int lepCategory)
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

float NtupleSFHelper::Get_LeptonIDSF(float lepPt, float lepEta, bool isEle, bool isMu,
							   bool isTight)
{
	assert(not _isdata);
	
	float sf = 1.;
	sf *= Get_LeptonSF_loose(lepPt, lepEta, isEle, isMu);
	// SF fakeable to loose is currently assumed to be 1.
	if (isTight)
		sf *= Get_LeptonSF_tight_vs_loose(lepPt, lepEta, isEle, isMu);

	return sf;
}

float NtupleSFHelper::Get_LeptonSF_loose(float lepPt, float lepEta,
								   bool isEle, bool isMu)
{
	float sf =1.;

	if (isMu) {
		if (abs(lepEta)<1.2) {
			sf *= evalTGraph(h_recoToLoose_leptonSF_mu1_b, lepPt);
		}
		else {
			sf *= evalTGraph(h_recoToLoose_leptonSF_mu1_e, lepPt);
		}
		
		sf *= read2DHist(h_recoToLoose_leptonSF_mu2,
						 lepPt, abs(lepEta));
		
		sf *= evalTGraph(h_recoToLoose_leptonSF_mu3, lepEta);
	}
	else if (isEle) {
		sf *= read2DHist(h_recoToLoose_leptonSF_el1,
						 lepPt, abs(lepEta));
		sf *= read2DHist(h_recoToLoose_leptonSF_el2,
						 lepPt, abs(lepEta));
		sf *= read2DHist(h_recoToLoose_leptonSF_el3,
						 lepPt, abs(lepEta));
		// !! different pt eta xaxis
		sf *= read2DHist(h_recoToLoose_leptonSF_gsf, lepEta, lepPt);
	}

	return sf;
}

float NtupleSFHelper::Get_LeptonSF_tight_vs_loose(float lepPt, float lepEta,
											bool isEle, bool isMu)
{
	float sf = 1.;
	
	if (_analysis == Analyze_2lss1tau) {
		if (isMu) {
			sf = read2DHist(h_looseToTight_leptonSF_mu_2lss,
							lepPt, abs(lepEta));
		}
		else if (isEle) {
			sf = read2DHist(h_looseToTight_leptonSF_el_2lss,
							lepPt, abs(lepEta));
		}
	}
	else if (_analysis == Analyze_3l) {
		if (isMu) {
			sf = read2DHist(h_looseToTight_leptonSF_mu_3l,
							lepPt, abs(lepEta));
		}
		else if (isEle) {
			sf = read2DHist(h_looseToTight_leptonSF_el_3l,
							lepPt, abs(lepEta));
		}
	}

	return sf;
}

float NtupleSFHelper::Get_PUWeight(int nPU)
{
	assert(not _isdata);
	
	int xbin = h_puweight->FindBin(nPU);
	return h_puweight->GetBinContent(xbin);
}

float NtupleSFHelper::Get_TauIDSF(float tauPt, float tauEta, bool isGenMatched)
{
	assert(not _isdata);

	// tau ID efficiency data/MC scale factor
	float tauEff_sf = 1.0;//0.9;
	// tau ID fake rate data/MC scale factor
	float tauFR_sf = tauEta < 1.497 ?
	    readTF(f_fakerate_tau_mvaM_etaL_ratio, tauPt) :
		readTF(f_fakerate_tau_mvaM_etaH_ratio, tauPt);

	return (isGenMatched ? tauEff_sf : tauFR_sf);
}

float NtupleSFHelper::Get_FakeRate(float lepConePt, float lepEta,
							 bool isEle, bool isMuon)
{
	assert(_selection == Control_1lfakeable);
	
	float fakerate = 0;
	
	if (isEle)
		fakerate = read2DHist(h_fakerate_el, lepConePt, lepEta);
	else if (isMuon)
		fakerate = read2DHist(h_fakerate_mu, lepConePt, lepEta);
	
	if (lepConePt < 10.) return 0.;
	
	return fakerate;
}

float NtupleSFHelper::Get_FakeRate(float tauPt, float tauEta)
{
	assert(not _isdata);
	
	float fr_mc = 0;
	float ratio = 0;

	if (abs(tauEta) < 1.479) {
		fr_mc = readTGraph(g_fakerate_tau_mvaM_etaL_mc, tauPt);
		ratio = readTF(f_fakerate_tau_mvaM_etaL_ratio, tauPt);
	}
	else {
		fr_mc = readTGraph(g_fakerate_tau_mvaM_etaH_mc, tauPt);
		ratio = readTF(f_fakerate_tau_mvaM_etaH_ratio, tauPt);
	}
	
	return fr_mc * ratio;
}

float NtupleSFHelper::Get_EleChargeMisIDProb(float elePt, float eleEta,
									   int eleCharge, int tauCharge)
{
	assert(_selection == Control_2los1tau);
	
	// only apply the charge flip rate to the electron that is same sign as tau
	// due to the tau charge requirement in signal region
	if (eleCharge * tauCharge < 0) return 0;

	return read2DHist(h_chargeMisId, elePt, abs(eleEta));
}

float NtupleSFHelper::read2DHist(TH2* h2d, float x, float y)
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

float NtupleSFHelper::evalTGraph(TGraphAsymmErrors* graph, float x)
{
	float x1 = std::max(float(graph->GetXaxis()->GetXmin()+1e-5),
						std::min(float(graph->GetXaxis()->GetXmax()-1e-5), x)
						);
	return graph->Eval(x1);
}

float NtupleSFHelper::readTGraph(TGraphAsymmErrors* graph, float x)
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

float NtupleSFHelper::readTF(TF1* f, float x)
{
	//float x1 = std::max(float(f->GetXaxis()->GetXmin()),
	//					std::min(float(f->GetXaxis()->GetXmax()), x)
	//					);
	return f->Eval(x);
}

#endif
