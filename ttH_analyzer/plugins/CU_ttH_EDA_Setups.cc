#ifndef CU_ttH_EDA_Setups_cc
#define CU_ttH_EDA_Setups_cc

/// Includes
#include "CU_ttH_EDA.h"

void CU_ttH_EDA::Close_output_files()
{
	if (analysis_type == Analyze_lepton_jet) {
		fclose(events_e_cut1);
		fclose(events_e_cut2);
		fclose(events_e_cut3);
		fclose(events_e_cut4);
		fclose(events_e_cut5);
		fclose(events_e_cut6);
		fclose(events_e_cut7);

		fclose(events_mu_cut1);
		fclose(events_mu_cut2);
		fclose(events_mu_cut3);
		fclose(events_mu_cut4);
		fclose(events_mu_cut5);
		fclose(events_mu_cut6);
		fclose(events_mu_cut7);
	}

	if (analysis_type == Analyze_dilepton) {
		fclose(events_dimu_cut1);
		fclose(events_dimu_cut2);
		fclose(events_dimu_cut3);
		fclose(events_dimu_cut4);
		fclose(events_dimu_cut5);
		fclose(events_dimu_cut6);
		fclose(events_dimu_cut7);

		fclose(events_diele_cut1);
		fclose(events_diele_cut2);
		fclose(events_diele_cut3);
		fclose(events_diele_cut4);
		fclose(events_diele_cut5);
		fclose(events_diele_cut6);
		fclose(events_diele_cut7);

		fclose(events_elemu_cut1);
		fclose(events_elemu_cut2);
		fclose(events_elemu_cut3);
		fclose(events_elemu_cut4);
		fclose(events_elemu_cut5);
	}
}

void CU_ttH_EDA::Set_up_histograms()
{

	if (analysis_type == Analyze_lepton_jet) {
		h_tth_syncex1_ele =
			fs_->make<TH1D>("h_tth_syncex1_ele", ";cut", 8, 0, 8);
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(2, "Single ele trig");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(3, "==1 electron");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(4, "==0 muons");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(5, ">=4 jets");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(6, ">=2 b-tags");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(7, ">=1 top-tags");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(8, ">=1 Higgs-tags");

		h_tth_syncex1_mu = fs_->make<TH1D>("h_tth_syncex1_mu", ";cut", 8, 0, 8);
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(2, "Single mu trig");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(3, "==1 muon");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(4, "==0 electrons");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(5, ">=4 jets");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(6, ">=2 b-tags");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(7, ">=1 top-tags");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(8, ">=1 Higgs-tags");
	}

	if (analysis_type == Analyze_dilepton) {
		h_tth_syncex1_dimu =
			fs_->make<TH1D>("h_tth_syncex1_dimu", ";cut", 8, 0, 8);
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(2, "Double mu trig");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(3, ">=2 muons");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(4, "Mll > 20");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(5, "Z Veto   ");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(6, ">=2 jets");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(7, "MET > 40");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(8, ">=1 b-tags");

		h_tth_syncex1_diele =
			fs_->make<TH1D>("h_tth_syncex1_diele", ";cut", 8, 0, 8);
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(2, "Double ele trig");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(3, ">=2 electrons");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(4, "Mll > 20");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(5, "Z Veto   ");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(6, ">=2 jets");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(7, "MET > 40");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(8, ">=1 b-tags");

		h_tth_syncex1_elemu =
			fs_->make<TH1D>("h_tth_syncex1_elemu", ";cut", 6, 0, 6);
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(2, "Ele-mu trig");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(3, ">=2 leptons");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(4, "Mll > 20");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(5, ">=2 jets");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(6, ">=1 b-tags");
	}

	if (analysis_type == Analyze_tau_ssleptons) {
		h_MVA_ttV_vs_ttbar =
			fs_->make<TH2D>("h_MVA_ttV_vs_ttbar", ";BDT", 20, -1, 1, 20, -1, 1);
		h_MVA_ttV_vs_ttbar->GetXaxis()->SetTitle("ttbar");
		h_MVA_ttV_vs_ttbar->GetYaxis()->SetTitle("ttV");

		h_MVA_shape =
			fs_->make<TH1D>("h_MVA_shape","", 6, 0.5, 6.5);

		if (!isdata and doSystematics) {
			h_MVA_ttV_vs_ttbar_jesup =
				fs_->make<TH2D>("h_MVA_ttV_vs_ttbar_JESUp", ";BDT", 20, -1, 1, 20, -1, 1);
			h_MVA_ttV_vs_ttbar_jesup->GetXaxis()->SetTitle("ttbar");
			h_MVA_ttV_vs_ttbar_jesup->GetYaxis()->SetTitle("ttV");
			
			h_MVA_shape_jesup =
				fs_->make<TH1D>("h_MVA_shape_JESUp","", 6, 0.5, 6.5);

			h_MVA_ttV_vs_ttbar_jesdown =
				fs_->make<TH2D>("h_MVA_ttV_vs_ttbar_JESDown", ";BDT", 20, -1, 1, 20, -1, 1);
			h_MVA_ttV_vs_ttbar_jesdown->GetXaxis()->SetTitle("ttbar");
			h_MVA_ttV_vs_ttbar_jesdown->GetYaxis()->SetTitle("ttV");
			
			h_MVA_shape_jesdown =
				fs_->make<TH1D>("h_MVA_shape_JESDown","", 6, 0.5, 6.5);
			
			for (int i = 0; i < 16 ; ++i) {
				TString h2d_name = "h_MVA_ttV_vs_ttbar_" + sysList[i];
				h_MVA_ttV_vs_ttbar_sys[i] =
					fs_->make<TH2D>(h2d_name, ";BDT", 20, -1, 1, 20, -1, 1);
				h_MVA_ttV_vs_ttbar_sys[i]->SetTitle("ttbar");
				h_MVA_ttV_vs_ttbar_sys[i]->SetTitle("ttV");

				TString h1d_name = "h_MVA_shape_" + sysList[i];
				h_MVA_shape_sys[i] = fs_->make<TH1D>(h1d_name, "", 6, 0.5, 6.5);
			}

			setup_sysHist = true;
		}

		if (selection_region == "control_1lfakeable" and isdata) {
			TFile* file_fr = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/FR_data_ttH_mva.root").c_str());
			h_fakerate_el = (TH2F*) file_fr->Get("FR_mva075_el_data_comb");
			h_fakerate_mu = (TH2F*) file_fr->Get("FR_mva075_mu_data_comb");
		}
	}

	if (analysis_type == Analyze_ditaus_lepton) {

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

/// Make and open write-out files
void CU_ttH_EDA::Set_up_output_files()
{
	if (analysis_type == Analyze_lepton_jet) {
		events_e_cut1 = fopen("Outputs/CU_events_e_cut1.dat", "w");
		events_e_cut2 = fopen("Outputs/CU_events_e_cut2.dat", "w");
		events_e_cut3 = fopen("Outputs/CU_events_e_cut3.dat", "w");
		events_e_cut4 = fopen("Outputs/CU_events_e_cut4.dat", "w");
		events_e_cut5 = fopen("Outputs/CU_events_e_cut5.dat", "w");
		events_e_cut6 = fopen("Outputs/CU_events_e_cut6.dat", "w");
		events_e_cut7 = fopen("Outputs/CU_events_e_cut7.dat", "w");

		events_mu_cut1 = fopen("Outputs/CU_events_mu_cut1.dat", "w");
		events_mu_cut2 = fopen("Outputs/CU_events_mu_cut2.dat", "w");
		events_mu_cut3 = fopen("Outputs/CU_events_mu_cut3.dat", "w");
		events_mu_cut4 = fopen("Outputs/CU_events_mu_cut4.dat", "w");
		events_mu_cut5 = fopen("Outputs/CU_events_mu_cut5.dat", "w");
		events_mu_cut6 = fopen("Outputs/CU_events_mu_cut6.dat", "w");
		events_mu_cut7 = fopen("Outputs/CU_events_mu_cut7.dat", "w");
	}

	if (analysis_type == Analyze_dilepton) {
		events_dimu_cut1 = fopen("Outputs/CU_events_dimu_cut1.dat", "w");
		events_dimu_cut2 = fopen("Outputs/CU_events_dimu_cut2.dat", "w");
		events_dimu_cut3 = fopen("Outputs/CU_events_dimu_cut3.dat", "w");
		events_dimu_cut4 = fopen("Outputs/CU_events_dimu_cut4.dat", "w");
		events_dimu_cut5 = fopen("Outputs/CU_events_dimu_cut5.dat", "w");
		events_dimu_cut6 = fopen("Outputs/CU_events_dimu_cut6.dat", "w");
		events_dimu_cut7 = fopen("Outputs/CU_events_dimu_cut7.dat", "w");

		events_diele_cut1 = fopen("Outputs/CU_events_diele_cut1.dat", "w");
		events_diele_cut2 = fopen("Outputs/CU_events_diele_cut2.dat", "w");
		events_diele_cut3 = fopen("Outputs/CU_events_diele_cut3.dat", "w");
		events_diele_cut4 = fopen("Outputs/CU_events_diele_cut4.dat", "w");
		events_diele_cut5 = fopen("Outputs/CU_events_diele_cut5.dat", "w");
		events_diele_cut6 = fopen("Outputs/CU_events_diele_cut6.dat", "w");
		events_diele_cut7 = fopen("Outputs/CU_events_diele_cut7.dat", "w");

		events_elemu_cut1 = fopen("Outputs/CU_events_elemu_cut1.dat", "w");
		events_elemu_cut2 = fopen("Outputs/CU_events_elemu_cut2.dat", "w");
		events_elemu_cut3 = fopen("Outputs/CU_events_elemu_cut3.dat", "w");
		events_elemu_cut4 = fopen("Outputs/CU_events_elemu_cut4.dat", "w");
		events_elemu_cut5 = fopen("Outputs/CU_events_elemu_cut5.dat", "w");
	}
}

void CU_ttH_EDA::Set_up_tokens(const edm::ParameterSet &config)
{
	token.event_gen_info =
		consumes<GenEventInfoProduct>(edm::InputTag(std::string("generator")));
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
	tauNtuple.set_up_branches(eventTree);
}

void CU_ttH_EDA::Set_up_BTagCalibration_Readers()
{
	const std::string base =
		std::string(getenv("CMSSW_BASE")) +  "/src/Analyzers/ttH_analyzer/data/";

	BTagCalibration calib_csvv2("csvv2", base + "CSVv2.csv");  // for 76X

	BTagCaliReaders["NA"] = new BTagCalibrationReader(&calib_csvv2,  // calibration instance
									   BTagEntry::OP_RESHAPING, // operating point
									   "iterativefit", // measurement type
									   "central");  // systematics type
	BTagCaliReaders["JESUp"] = 
		new BTagCalibrationReader(&calib_csvv2,	BTagEntry::OP_RESHAPING,
								  "iterativefit", "up_jes");
	BTagCaliReaders["JESDown"] = 
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "down_jes");
	BTagCaliReaders["LFUp"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "up_lf");
	BTagCaliReaders["LFDown"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "down_lf");
	BTagCaliReaders["HFUp"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "up_hf");
	BTagCaliReaders["HFDown"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "down_hf");
	BTagCaliReaders["HFStats1Up"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "up_hfstats1");
	BTagCaliReaders["HFStats1Down"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "down_hfstats1");	
	BTagCaliReaders["HFStats2Up"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "up_hfstats2");
	BTagCaliReaders["HFStats2Down"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "down_hfstats2");	
	BTagCaliReaders["LFStats1Up"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "up_lfstats1");
	BTagCaliReaders["LFStats1Down"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "down_lfstats1");	
	BTagCaliReaders["LFStats2Up"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "up_lfstats2");
	BTagCaliReaders["LFStats2Down"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "down_lfstats2");
	BTagCaliReaders["cErr1Up"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "up_cferr1");
	BTagCaliReaders["cErr1Down"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "down_cferr1");
	BTagCaliReaders["cErr2Up"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "up_cferr2");
	BTagCaliReaders["cErr2Down"] =
		new BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING,
								  "iterativefit", "down_cferr2");	
}

void CU_ttH_EDA::Set_up_CSV_rootFile()
{
	std::string inputFileHF = "data/csv_rwt_fit_hf_76x_2016_02_08.root";
	std::string inputFileLF = "data/csv_rwt_fit_lf_76x_2016_02_08.root";

	TFile* f_CSVwgt_HF = new TFile ((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/" + inputFileHF).c_str());
	TFile* f_CSVwgt_LF = new TFile ((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/" + inputFileLF).c_str());

	fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);
}

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

		for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );
		
		if( iSys<5 ){
			for( int iPt=0; iPt<5; iPt++ ) c_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
		}
		
		for( int iPt=0; iPt<4; iPt++ ){
			for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
		}
	}

	return;
}

#endif
	
