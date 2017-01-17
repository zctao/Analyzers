#ifndef CU_ttH_EDA_Setups_cc
#define CU_ttH_EDA_Setups_cc

/// Includes
#include "CU_ttH_EDA.h"

void CU_ttH_EDA::Set_up_histograms()
{

	// event count histograms
	h_nProcessed = fs_->make<TH1I>("h_nProcessed","",1,0.5,1.5);
	h_SumGenWeight = fs_->make<TH1D>("h_SumGenWeight","",1,0.5,1.5);
	h_SumGenWeightxPU = fs_->make<TH1D>("h_SumGenWeightxPU","",1,0.5,1.5);
	
	// categories
	TString lep_cat[3] = {"mumu", "ee", "emu"};
	TString btag_cat[2] = {"bloose", "bmedium"};
	
	if (analysis_type == Analyze_2lss1tau) {
		
		for (int il = 0; il < 3; il++) {
			for (int ib = 0; ib <2; ib++) {
				TString h2d_name =
					"h_MVA_ttV_vs_ttbar_"+lep_cat[il]+"_"+btag_cat[ib];
				
				h_MVA_ttV_vs_ttbar[il][ib] =
					fs_->make<TH2D>(h2d_name, ";BDT", 20, -1, 1, 20, -1, 1);
				h_MVA_ttV_vs_ttbar[il][ib]->GetXaxis()->SetTitle("ttbar");
				h_MVA_ttV_vs_ttbar[il][ib]->GetYaxis()->SetTitle("ttV");

				TString hshape_name = 
					"h_MVA_shape_"+lep_cat[il]+"_"+btag_cat[ib];
				
				h_MVA_shape[il][ib] =
					fs_->make<TH1D>(hshape_name, "", 7, 0.5, 7.5);

				// systemtics
				if (!isdata and doSystematics) {

					for (int icsv=0; icsv < 16; ++icsv) {
						//TString h2d_csv_name =
						//	"h_MVA_ttV_vs_ttbar_"+lep_cat[il]+"_"+btag_cat[ib]+"_"
						//	+sysList[icsv];
						//h_MVA_ttV_vs_ttbar_sys[il][ib][icsv] =
						//	fs_->make<TH2D>(h2d_csv_name,";BDT",20,-1,1,20,-1,1);
						//h_MVA_ttV_vs_ttbar_sys[il][ib][icsv]->SetTitle("ttbar");
						//h_MVA_ttV_vs_ttbar_sys[il][ib][icsv]->SetTitle("ttV");

						TString hshape_csv_name =
							"h_MVA_shape_"+lep_cat[il]+"_"+btag_cat[ib]+"_"
							+sysList[icsv];

						h_MVA_shape_csv_sys[il][ib][icsv] =
							fs_->make<TH1D>(hshape_csv_name, "", 7, 0.5, 7.5);
					}
					
					TString thuName[4] = {"xDown","xUp","yDown","yUp"}; //to check
					for (int ithu=0; ithu < 4; ++ithu) {
						TString hshape_thu_name =
							"h_MVA_shape_"+lep_cat[il]+"_"+btag_cat[ib]+"_"
							+thuName[ithu];
						h_MVA_shape_thu_sys[il][ib][ithu] =
							fs_->make<TH1D>(hshape_thu_name, "", 7, 0.5, 7.5);
					}
					
					setup_sysHist = true;
				}
			}
		}

	}

	if (analysis_type == Analyze_3l) {
		h_mTWl = fs_->make<TH1D>("h_mTWl","", 10, 0, 300);
		h_met = fs_->make<TH1D>("h_met", "", 8, 0, 200);
		h_nJets = fs_->make<TH1D>("h_nJets", "", 8, 0, 8);
		h_lep_charge = fs_->make<TH1D>("h_lep_charge", "", 10, -5, 5);
		h_mZ = fs_->make<TH1D>("h_mZ", "", 15, 60, 120);
	}
}

void CU_ttH_EDA::Set_up_FakeRate_Lut()
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

void CU_ttH_EDA::Set_up_TauSF_Lut()
{
	file_fr_tau = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/FR_tau_2016.root").c_str(), "read");
	
	f_fakerate_tau_mvaM_etaL_ratio = (TF1*) file_fr_tau->Get("jetToTauFakeRate/dR03mvaMedium/absEtaLt1_5/fitFunction_data_div_mc_hadTaus_pt");
	f_fakerate_tau_mvaM_etaH_ratio = (TF1*) file_fr_tau->Get("jetToTauFakeRate/dR03mvaMedium/absEta1_5to9_9/fitFunction_data_div_mc_hadTaus_pt");
}

void CU_ttH_EDA::Set_up_ChargeMisID_Lut()
{
	file_eleMisCharge = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/QF_data_el.root").c_str(), "read");
	h_chargeMisId = (TH2F*) file_eleMisCharge->Get("chargeMisId");
}

void CU_ttH_EDA::Set_up_PUWeight_hist()
{
	file_puweight = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/PU_weights/PU_weights_2016_271036_284044.root").c_str(), "read");
	h_puweight = (TH1F*) file_puweight->Get("h_ratio_data_MC");
}

void CU_ttH_EDA::Set_up_LeptonSF_Lut()
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
	if (analysis_type == Analyze_2lss1tau) {
		/// for muon
		file_looseToTight_leptonSF_mu_2lss = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/lepMVAEffSF_m_2lss.root").c_str(),"read");
		h_looseToTight_leptonSF_mu_2lss = (TH2F*)(file_looseToTight_leptonSF_mu_2lss->Get("sf"));
		
		/// for electron
		file_looseToTight_leptonSF_el_2lss = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/lepMVAEffSF_e_2lss.root").c_str(),"read");
		h_looseToTight_leptonSF_el_2lss = (TH2F*)(file_looseToTight_leptonSF_el_2lss->Get("sf"));
	}
	
	if (analysis_type == Analyze_3l) {
		/// for muon
		file_looseToTight_leptonSF_mu_3l = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/lepMVAEffSF_m_3l.root").c_str(),"read");
		h_looseToTight_leptonSF_mu_3l = (TH2F*)(file_looseToTight_leptonSF_mu_3l->Get("sf"));
		
		/// for electron
		file_looseToTight_leptonSF_el_3l = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/leptonSF/lepMVAEffSF_e_3l.root").c_str(),"read");
		h_looseToTight_leptonSF_el_3l = (TH2F*)(file_looseToTight_leptonSF_el_3l->Get("sf"));
	}
}

/// Prepare histograms for trigger/filter counts
int CU_ttH_EDA::Set_up_Run_histograms_triggers()
{
	unsigned int numHLT = trigger_names_no_ver.size();
	h_hlt = fs_->make<TH1D>("h_hlt", ";HLT path", numHLT, 0, numHLT);
	if (!h_hlt)
		return 1;

	TAxis *axis = h_hlt->GetXaxis();
	if (!axis)
		return 1;

	for (unsigned int i = 0; i < numHLT; ++i)
		axis->SetBinLabel(i + 1, trigger_names_no_ver[i].c_str());

	unsigned int numFLT = filter_names_no_ver.size();
	h_flt = fs_->make<TH1D>("h_flt", ";Filter path", numFLT, 0, numFLT);
	if (!h_flt)
		return 1;

	axis = h_flt->GetXaxis();
	if (!axis)
		return 1;

	for (unsigned int i = 0; i < numFLT; ++i)
		axis->SetBinLabel(i + 1, filter_names_no_ver[i].c_str());

	return 0;
}

void CU_ttH_EDA::Set_up_HLT_path_name()
{
	std::vector<std::string> triggerNames = hlt_config.triggerNames();

	Add_HLT_pathName_version(trigger_on_HLT_e, triggerNames);
	Add_HLT_pathName_version(trigger_on_HLT_mu, triggerNames);
	Add_HLT_pathName_version(trigger_on_HLT_ee, triggerNames);
	Add_HLT_pathName_version(trigger_on_HLT_emu, triggerNames);
	Add_HLT_pathName_version(trigger_on_HLT_mumu, triggerNames);
	
}

void CU_ttH_EDA::Add_HLT_pathName_version(std::vector<std::string>& triggers,
										  std::vector<std::string>& triggerNames)
{
	for (std::string & hltPath : triggers) {
		
		bool foundPath = false;

		for (const std::string & pathName : triggerNames) {
			//std::string pathNameNoVer = hlt_config.removeVersion(pathName);
			if (pathName.find(hltPath) != std::string::npos) {
				hltPath = pathName;
				foundPath = true;
				break;
			}
		}

		if (not foundPath) {
			std::cerr <<"WARNING!! Cannot find HLT path "<< hltPath << std::endl;
		}
	}
}

void CU_ttH_EDA::Set_up_trigger_name_vectors()
{
	/// Fill trigger name vectors and counters
	trigger_names = hlt_config.triggerNames();

	trigger_names_no_ver.clear();
	trigger_names_no_ver.push_back("All");
	std::string prefix = "HLT_";
	for (unsigned int i = 0; i < trigger_names.size(); ++i) {
		std::string pathNameNoVer = hlt_config.removeVersion(trigger_names[i]);

		n_trigger_fired[pathNameNoVer] = 0;

		if (trigger_names[i].compare(0, prefix.length(), prefix) == 0)
			trigger_names_no_ver.push_back(pathNameNoVer);
	}

	/// Fill filter name vectors and counters
	filter_names = filter_config.triggerNames();

	filter_names_no_ver.clear();
	filter_names_no_ver.push_back("All");
	for (unsigned int i = 0; i < filter_names.size(); ++i) {
		std::string pathNameNoVer =
			filter_config.removeVersion(filter_names[i]);

		n_filter_fired[pathNameNoVer] = 0;

		filter_names_no_ver.push_back(pathNameNoVer);
	}
}


void CU_ttH_EDA::Set_up_tokens(const edm::ParameterSet &config)
{
	token.event_gen_info =
		consumes<GenEventInfoProduct>(edm::InputTag(std::string("generator")));
	token.event_lhe_info =
		consumes<LHEEventProduct>(edm::InputTag(std::string("externalLHEProducer")));
	token.triggerResults = consumes<edm::TriggerResults>(
		edm::InputTag(std::string("TriggerResults"), std::string(""), hltTag));
	token.filterResults = consumes<edm::TriggerResults>(edm::InputTag(
		std::string("TriggerResults"), std::string(""), filterTag));

	token.vertices = consumes<reco::VertexCollection>(
	    config.getParameter<edm::InputTag>("pv"));
	token.sec_vertices = consumes<reco::VertexCompositePtrCandidateCollection>(
	    config.getParameter<edm::InputTag>("sv"));
	token.PU_info = consumes<std::vector<PileupSummaryInfo>>(
	    config.getParameter<edm::InputTag>("pileup"));
	token.srcRho = consumes<double>(
	    config.getParameter<edm::InputTag>("rho"));
	token.electrons = consumes<pat::ElectronCollection>(
	    config.getParameter<edm::InputTag>("electrons"));
	token.muons = consumes<pat::MuonCollection>(
		config.getParameter<edm::InputTag>("muons"));
	token.taus = consumes<pat::TauCollection>(
	    config.getParameter<edm::InputTag>("taus"));
	token.jets = consumes<pat::JetCollection>(
	    config.getParameter<edm::InputTag>("jets"));
	token.METs = consumes<pat::METCollection>(
	    config.getParameter<edm::InputTag>("mets"));
	token.PF_candidates = consumes<pat::PackedCandidateCollection>(
	    config.getParameter<edm::InputTag>("pfcand"));
	token.BS = consumes<reco::BeamSpot>(
	    config.getParameter<edm::InputTag>("beamspot"));
	//token.top_jets = consumes<boosted::HTTTopJetCollection>(
	//	edm::InputTag("HTTTopJetsPFMatcher", "heptopjets", "p"));
	//token.subfilter_jets = consumes<boosted::SubFilterJetCollection>(
	//	edm::InputTag("CA12JetsCA3FilterjetsPFMatcher", "subfilterjets", "p"));	
	token.MC_particles = consumes<reco::GenParticleCollection>(
	    config.getParameter<edm::InputTag>("prunedgen"));
	token.MC_packed = consumes<pat::PackedGenParticleCollection>(
	    config.getParameter<edm::InputTag>("packedgen"));
	token.genJets = consumes<reco::GenJetCollection>(
	    edm::InputTag(std::string("slimmedGenJets")));
}

void CU_ttH_EDA::Set_up_selection_region(const string & selection_region )
{
	if (selection_region == "signal_2lss1tau") {
		selection_type = Signal_2lss1tau;
		return;
	}

	if (selection_region == "signal_1l2tau") {
		selection_type = Signal_1l2tau;
		return;
	}

	if (selection_region == "signal_3l") {
		selection_type = Signal_3l;
		return;
	}

	if (selection_region == "control_2los1tau") {
		selection_type = Control_2los1tau;
		return;
	}

	if (selection_region == "control_1lfakeable") {
		selection_type = Control_1lfakeable;
		return;
	}

	if (selection_region == "control_WZ") {
		selection_type = Control_WZ;
		return;
	}

	std::cout << "Not valid selection region !" << std::endl;
	assert(0);
	return;
		
}

void CU_ttH_EDA::Set_up_Tree()
{
	eventTree = fs_->make<TTree>("eventTree", "Event tree");
	
	/*
	// If ntuple class is inherited from Root TClass (ToDo)
	//
	//eventTree -> Branch("ntuple_", "CU_ttH_EDA_Ntuple", &ntuple);
	//std::cout << "IsTObject :" << ntuple->IsTObject() <<  std::endl;
	//std::cout << "GetNdata() :" << ntuple->GetNdata() << std::endl;
	//std::cout << "CanSplit() :" << ntuple->CanSplit() << std::endl;
	//ntuple->Dump();
	*/
	evtNtuple.set_up_branches(eventTree);
}


void CU_ttH_EDA::Set_up_BTagCalibration_Readers()
{
	const std::string base =
		std::string(getenv("CMSSW_BASE")) +  "/src/Analyzers/ttH_analyzer/data/";

	BTagCalibration calib_csvv2("csvv2", base + "CSVv2_ichep.csv");
	
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

/*
void CU_ttH_EDA::Set_up_CSV_rootFile()
{
	std::string inputFileHF = "data/csv_rwt_fit_hf_v2_final_2016_09_7test.root";
	std::string inputFileLF = "data/csv_rwt_fit_lf_v2_final_2016_09_7test.root";

	f_CSVwgt_HF = new TFile ((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/" + inputFileHF).c_str());
	f_CSVwgt_LF = new TFile ((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/" + inputFileLF).c_str());

	fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);
}
*/
/*
void CU_ttH_EDA::fillCSVhistos(TFile* fileHF, TFile* fileLF)
{
	for( int iSys=0; iSys<9; iSys++ ){
		for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = NULL;
		for( int iPt=0; iPt<3; iPt++ ){
			for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
		}
	}
	for( int iSys=0; iSys<5; iSys++ ){
		for( int iPt=0; iPt<5; iPt++ ) c_csv_wgt_hf[iSys][iPt] = NULL;
	}

	// CSV reweighting /// only care about the nominal ones
	for( int iSys=0; iSys<9; iSys++ ){
		TString syst_csv_suffix_hf = "final";
		TString syst_csv_suffix_c = "final";
		TString syst_csv_suffix_lf = "final";

		switch( iSys ){
		case 0:
			// this is the nominal case
			break;
		case 1:
			// JESUp
			syst_csv_suffix_hf = "final_JESUp"; syst_csv_suffix_lf = "final_JESUp";
			syst_csv_suffix_c  = "final_cErr1Up";
			break;
		case 2:
			// JESDown
			syst_csv_suffix_hf = "final_JESDown"; syst_csv_suffix_lf = "final_JESDown";
			syst_csv_suffix_c  = "final_cErr1Down";
			break;
		case 3:
			// purity up
			syst_csv_suffix_hf = "final_LFUp"; syst_csv_suffix_lf = "final_HFUp";
			syst_csv_suffix_c  = "final_cErr2Up";
			break;
		case 4:
			// purity down
			syst_csv_suffix_hf = "final_LFDown"; syst_csv_suffix_lf = "final_HFDown";
			syst_csv_suffix_c  = "final_cErr2Down";
			break;
		case 5:
			// stats1 up
			syst_csv_suffix_hf = "final_Stats1Up"; syst_csv_suffix_lf = "final_Stats1Up";
			break;
		case 6:
			// stats1 down
			syst_csv_suffix_hf = "final_Stats1Down"; syst_csv_suffix_lf = "final_Stats1Down";
			break;
		case 7:
			// stats2 up
			syst_csv_suffix_hf = "final_Stats2Up"; syst_csv_suffix_lf = "final_Stats2Up";
			break;
		case 8:
			// stats2 down
			syst_csv_suffix_hf = "final_Stats2Down"; syst_csv_suffix_lf = "final_Stats2Down";
			break;
		}

		for( int iPt=0; iPt<5; iPt++ ) {
			h_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );
		}
		
		if( iSys<5 ){
			for( int iPt=0; iPt<5; iPt++ ) {
				c_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
			}
		}
		
		for( int iPt=0; iPt<4; iPt++ ){
			for( int iEta=0; iEta<3; iEta++ ) {
				h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
			}
		}
	}

	return;
}
*/
#endif
	
