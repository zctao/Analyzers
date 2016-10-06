#ifndef CU_ttH_EDA_Setups_cc
#define CU_ttH_EDA_Setups_cc

/// Includes
#include "CU_ttH_EDA.h"

void CU_ttH_EDA::Set_up_histograms()
{

	// categories
	TString lep_cat[3] = {"mumu", "ee", "emu"};
	TString btag_cat[2] = {"bloose", "bmedium"};
	
	if (analysis_type == Analyze_2lss1tau) {
		h_nProcessed = fs_->make<TH1I>("h_nProcessed","",1,0,1);

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

						h_MVA_shape_sys[il][ib][icsv] =
							fs_->make<TH1D>(hshape_csv_name, "", 7, 0.5, 7.5);
					}

					setup_sysHist = true;
				}
			}
		}


		if (selection_type == Control_1lfakeable and isdata) {
			TFile* file_fr = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/FR_data_ttH_mva.root").c_str());
			
			h_fakerate_el = (TH2F*) file_fr->Get("FR_mva075_el_data_comb");
			h_fakerate_mu = (TH2F*) file_fr->Get("FR_mva075_mu_data_comb");
			h_fakerate_el -> SetDirectory(0);
			h_fakerate_mu -> SetDirectory(0);

			delete file_fr;
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

/*
void CU_ttH_EDA::Set_up_BTagCalibration_Readers()
{
	const std::string base =
		std::string(getenv("CMSSW_BASE")) +  "/src/Analyzers/ttH_analyzer/data/";

	BTagCalibration calib_csvv2("csvv2", base + "CSVv2_ichep.csv");  // 80X
	
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
*/

void CU_ttH_EDA::Set_up_CSV_rootFile()
{
	//std::string inputFileHF = "data/csv_rwt_fit_hf_76x_2016_02_08.root";
	//std::string inputFileLF = "data/csv_rwt_fit_lf_76x_2016_02_08.root";
	std::string inputFileHF = "data/csv_rwt_fit_hf_v2_final_2016_08_5test.root";
	std::string inputFileLF = "data/csv_rwt_fit_lf_v2_final_2016_08_5test.root";

	TFile* f_CSVwgt_HF = new TFile ((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/" + inputFileHF).c_str());
	TFile* f_CSVwgt_LF = new TFile ((std::string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/" + inputFileLF).c_str());

	fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);

	delete f_CSVwgt_HF;
	delete f_CSVwgt_LF;
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

		for( int iPt=0; iPt<5; iPt++ ) {
			h_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );
			h_csv_wgt_hf[iSys][iPt] -> SetDirectory(0);
		}
		
		if( iSys<5 ){
			for( int iPt=0; iPt<5; iPt++ ) {
				c_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
				c_csv_wgt_hf[iSys][iPt] -> SetDirectory(0);
			}
		}
		
		for( int iPt=0; iPt<4; iPt++ ){
			for( int iEta=0; iEta<3; iEta++ ) {
				h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
				h_csv_wgt_lf[iSys][iPt][iEta] -> SetDirectory(0);
			}
		}
	}

	return;
}

#endif
	
