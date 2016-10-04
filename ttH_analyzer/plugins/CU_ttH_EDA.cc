#ifndef CU_ttH_EDA_cc
#define CU_ttH_EDA_cc

// -*- C++ -*-
//
// Package:    ttH-LeptonPlusJets/CU_ttH_EDA
// Class:      CU_ttH_EDA
//
/**\class CU_ttH_EDA CU_ttH_EDA.cc
 ttH-LeptonPlusJets/AnalysisCode/plugins/CU_ttH_EDA.cc

 Description: [one line class summary]

 Implementation:
		 [Notes on implementation]
*/

/// Includes
#include "CU_ttH_EDA.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

/*
 * Function/method section
*/

/// Constructor
CU_ttH_EDA::CU_ttH_EDA(const edm::ParameterSet &iConfig):
	// Analysis type
	config_analysis_type (iConfig.getParameter<string>("analysis_type")),
	// Sync ntuple
	produce_sync_ntuple (iConfig.getParameter<bool>("produce_sync_ntuple")),
	// Systematics
	doSystematics (iConfig.getParameter<bool>("do_systematics")),
	// Sample parameter
	doLumiScale (iConfig.getParameter<bool>("doLumiScale")),
	sampleName (iConfig.getParameter<string>("sampleName")),
	sample_xs (iConfig.getParameter<double>("sample_xs")),
	int_lumi (iConfig.getParameter<double>("int_lumi")),
	// Generic
	verbose_ (iConfig.getParameter<bool>("verbosity")),
	dumpHLT_ (iConfig.getParameter<bool>("print_HLT_event_path")),
	hltcut_off (iConfig.getParameter<bool>("turn_off_HLT_cut")),
	hltTag (iConfig.getParameter<string>("HLT_config_tag")),
	filterTag (iConfig.getParameter<string>("filter_config_tag")),
	// Triggers
	trigger_stats (iConfig.getParameter<bool>("collect_trigger_stats")),
	trigger_on_HLT_e (iConfig.getParameter<std::vector<string>>("HLT_electron_triggers")),
	trigger_on_HLT_mu (iConfig.getParameter<std::vector<string>>("HLT_muon_triggers")),
	trigger_on_HLT_ee (iConfig.getParameter<std::vector<string>>("HLT_electron_electron_triggers")),
	trigger_on_HLT_emu (iConfig.getParameter<std::vector<string>>("HLT_electron_muon_triggers")),
	trigger_on_HLT_mumu (iConfig.getParameter<std::vector<string>>("HLT_muon_muon_triggers")),
	// Cuts
	min_tight_lepton_pT (iConfig.getParameter<double>("min_tight_lepton_pT")),
	min_ele_pT (iConfig.getParameter<double>("min_ele_pT")),
	min_mu_pT (iConfig.getParameter<double>("min_mu_pT")),
	min_tau_pT (iConfig.getParameter<double>("min_tau_pT")),
	min_jet_pT (iConfig.getParameter<double>("min_jet_pT")),
	min_bjet_pT (iConfig.getParameter<double>("min_bjet_pT")),
	max_jet_eta (iConfig.getParameter<double>("max_jet_eta")),
	max_bjet_eta (iConfig.getParameter<double>("max_bjet_eta")),
	min_njets (iConfig.getParameter<int>("min_njets")),
	min_nbtags (iConfig.getParameter<int>("min_nbtags")),
	// JEC
	JECType (iConfig.getParameter<string>("JECType")),
	//jet_corrector (iConfig.getParameter<string>("jet_corrector")),
	isdata (iConfig.getParameter<bool>("using_real_data")),
	selection_region (iConfig.getParameter<string>("selection_region"))
{
	/*
	 * now do whatever initialization is needed
	*/

	/// temporary mock-up parameters
	MAODHelper_era = "2015_74x";
	MAODHelper_sample_nr = 2500;

	Load_configuration_set_type(config_analysis_type);	
	Load_configuration_MAODH(isdata);
	
	// Load_configuration(static_cast<string>("Configs/config_analyzer.yaml"));

	Set_up_selection_region(selection_region);
	
	miniAODhelper.UseCorrectedJets();
	
	Set_up_tokens(iConfig.getParameter<edm::ParameterSet>("input_tags"));
		
	Set_up_histograms();

	Set_up_Tree();

	//Set_up_BTagCalibration_Readers();
	Set_up_CSV_rootFile();

}

/// Destructor
CU_ttH_EDA::~CU_ttH_EDA()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
	//delete BTagCaliReader;
}

// ------------ method called for each event  ------------
void CU_ttH_EDA::analyze(const edm::Event &iEvent,
						 const edm::EventSetup &iSetup)
{
	
	using namespace edm;
	++event_count;
	h_nProcessed->Fill(1);

	/// Create and set up edm:Handles in stack mem.
	edm_Handles handle;
	Set_up_handles(iEvent, handle, token, isdata);
	
	// for signal sample, filter Higgs decay mode
	if (sampleName.Contains("ttH_")) {
		bool rightHiggsDecay = HiggsDecayFilter(*(handle.MC_particles), sampleName);
		if (not rightHiggsDecay)
			return;
	}
	
	/// Declaring local struct for data readout and manipulations
	CU_ttH_EDA_event_vars local;

	/// Triggers have not fired yet. Check_triggers, Check_filters will adjust
	local.pass_single_e = false;
	local.pass_single_mu = false;
	local.pass_double_mu = false;
	local.pass_double_e = false;
	local.pass_elemu = false;
	Update_common_vars(iEvent, local);
	
	//Initialize weights
	local.weight = 1.;
	local.csv_weight = 1.;
	//local.pu_weight = 1.;
	local.gen_weight = 1.;
	local.hlt_sf = 1.;
	local.lepIDEff_sf = 1.;
	
	/// Run checks on event containers via their handles
	Check_triggers(handle.triggerResults, local);
	Check_filters(handle.filterResults);
	
	Check_vertices_set_MAODhelper(handle.vertices);
	// 	Check_beam_spot(BS);	// dumb implementation

	// Setting rho
	auto rho = handle.srcRho;
	miniAODhelper.SetRho(*rho);

	/// Get and set miniAODhelper's jet corrector from the event setup
	//miniAODhelper.SetJetCorrector(
	//	JetCorrector::getJetCorrector(jet_corrector, iSetup));

	edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
	iSetup.get<JetCorrectionsRecord>().get("AK4PF",JetCorParColl);
	const JetCorrectorParameters & JetCorPar = (*JetCorParColl)["Uncertainty"];
	miniAODhelper.SetJetCorrectorUncertainty(JetCorPar);

	if (not isdata)
		local.gen_weight = handle.event_gen_info.product()->weight();

	if (trigger_stats) {
		h_hlt->Fill(0., 1);
		h_flt->Fill(0., 1);
	}

	/*
	  Muons
	*/
	for (const auto& mu : *(handle.muons)){
		if (mu.userFloat("idPreselection") > 0.5 and mu.pt() > min_mu_pT)
			local.mu_preselected.push_back(mu);
	}
	// sort by pT
	local.mu_preselected_sorted = miniAODhelper.GetSortedByPt(local.mu_preselected);
	// loose, fakeable and tight
	// loose = preselected
	for (auto & mu : local.mu_preselected_sorted) {
		miniLepton lepton(mu);
		local.leptons_loose.push_back(lepton);
		
		if (lepton.passFakeableSel()) {
			local.mu_fakeable.push_back(mu);
			local.leptons_fakeable.push_back(lepton);
		
			if (lepton.passTightSel()) {
				local.mu_tight.push_back(mu);
				local.leptons_tight.push_back(lepton);
			}
		}
	}

	/*
	  Electrons
	*/
	for (const auto& ele : *(handle.electrons)) {
		if (ele.userFloat("idPreselection") > 0.5 and ele.pt() > min_ele_pT)
			local.e_preselected.push_back(ele);
	}
	// remove overlap with muons
	local.e_preselected =
		removeOverlapdR(local.e_preselected, local.mu_preselected, 0.05);
	// sort by pT
	local.e_preselected_sorted = miniAODhelper.GetSortedByPt(local.e_preselected);
	// loose, fakeable and tight
	// loose = preselected
	for (auto & ele : local.e_preselected_sorted) {
		miniLepton lepton(ele);
		local.leptons_loose.push_back(lepton);
		
		if (lepton.passFakeableSel()) {
			local.e_fakeable.push_back(ele);
			local.leptons_fakeable.push_back(lepton);
		
			if (lepton.passTightSel()) {
				local.e_tight.push_back(ele);
				local.leptons_tight.push_back(lepton);
			}
		}
	}

	std::sort(local.leptons_fakeable.begin(), local.leptons_fakeable.end(), [] (miniLepton l1, miniLepton l2) { return ptr(l1)->conePt() > ptr(l2)->conePt();});
	std::sort(local.leptons_tight.begin(), local.leptons_tight.end(), [] (miniLepton l1, miniLepton l2) { return ptr(l1)->conePt() > ptr(l2)->conePt();});
	
	/*
	  Taus
	*/
	for (const auto& tau : *(handle.taus)) {
		if (tau.userFloat("idPreselection")>0.5 and tau.pt()>min_tau_pT)
			local.tau_preselected.push_back(tau);  // for cleaning
		if (tau.userFloat("idSelection")>0.5 and tau.pt()>min_tau_pT)
			local.tau_selected.push_back(tau);
	}
	// remove overlap with muons and electrons
	local.tau_selected =
		removeOverlapdR(local.tau_selected, local.e_preselected, 0.4);
	local.tau_selected =
		removeOverlapdR(local.tau_selected, local.mu_preselected, 0.4);
	// sort by pT
	local.tau_selected_sorted = miniAODhelper.GetSortedByPt(local.tau_selected);

	
	// number of selected leptons
	local.n_electrons_loose = static_cast<int>(local.e_preselected.size());
	local.n_electrons_fakeable = static_cast<int>(local.e_fakeable.size());
	local.n_electrons_tight = static_cast<int>(local.e_tight.size());
	local.n_muons_loose = static_cast<int>(local.mu_preselected.size());
	local.n_muons_fakeable = static_cast<int>(local.mu_fakeable.size());
	local.n_muons_tight = static_cast<int>(local.mu_tight.size());
	local.n_taus_pre = static_cast<int>(local.tau_preselected.size());
	local.n_taus = static_cast<int>(local.tau_selected.size());
									
	/*
	  Jets selection
	*/
	local.jets_corrected =
		miniAODhelper.GetCorrectedJets(*(handle.jets),
									   iEvent,iSetup, JECTypes[JECType]);
	local.jets_raw = miniAODhelper.GetSelectedJets(
		local.jets_corrected, min_jet_pT, max_jet_eta, jetID::jetLoose, '-');
	
	// overlap removal by dR
	local.jets_no_mu = removeOverlapdR(local.jets_raw, local.mu_fakeable, 0.4);
	local.jets_no_mu_e = removeOverlapdR(local.jets_no_mu, local.e_fakeable, 0.4);
	local.jets_selected = removeOverlapdR(local.jets_no_mu_e, local.tau_preselected, 0.4);
	
	// sort by pT
	local.jets_selected_sorted =
		miniAODhelper.GetSortedByPt(local.jets_selected);

	/*
	  BTag
	*/
	local.jets_selected_btag_loose = miniAODhelper.GetSelectedJets(
		local.jets_selected_sorted, min_bjet_pT, max_bjet_eta, jetID::jetLoose,
		'L');
	local.jets_selected_btag_medium = miniAODhelper.GetSelectedJets(
	    local.jets_selected_btag_loose, min_bjet_pT, max_bjet_eta,
		jetID::jetLoose, 'M');

	local.n_jets = static_cast<int>(local.jets_selected.size());
	local.n_btags_loose = static_cast<int>(local.jets_selected_btag_loose.size());
	local.n_btags_medium = static_cast<int>(local.jets_selected_btag_medium.size());


	/*
	  MET
	*/
	// TODO: need to propagate JEC and uncertainties to MET
	local.pfMET = handle.METs->front();
	// MHT
	float mht = getMHT(local);
	float met = local.pfMET.pt();
	float metld = 0.00397 * met + 0.00265 * mht;
	local.MHT = mht;
	local.metLD = metld;
	
	/*
	/// Get Corrected MET, !!!not yet used!!!
	// may need to be placed in CU_ttH_EDA_event_vars
	local.MET_corrected =
		handle.METs->front(); // miniAODhelper.GetCorrectedMET( METs.at(0),
							  // pfJets_forMET, iSysType );
	*/
	
	if (analysis_type == Analyze_2lss1tau) {
		
		bool passHLT =
			local.pass_single_e or local.pass_single_mu or local.pass_double_mu or
			local.pass_double_e or local.pass_elemu;
		// Check if all HLTs failed for debug purpose: assert(not passHLT);

		if (hltcut_off) passHLT = true;
					
		// Event selection
		bool pass_event_selection =
			pass_event_sel_2l(local, selection_type)
			and passHLT;

		evtNtuple.initialize();
		
		if (pass_event_selection) {
			assert(local.leptons_fakeable.size() == 2);
				
			// MVA
			MVA_ttbar_vars.Calculate_mvaVars(local.leptons_fakeable,
											 local.tau_selected_sorted,
											 local.jets_selected_sorted,
											 local.pfMET);
			MVA_ttV_vars.Calculate_mvaVars(local.leptons_fakeable,
										   local.tau_selected_sorted,
										   local.jets_selected_sorted,
										   local.pfMET);
				
			double mva_ttar = MVA_ttbar_vars.Get_mvaScore();
			double mva_ttV = MVA_ttV_vars.Get_mvaScore();

			// Write MVA variables to ntuple
			evtNtuple.MVA_2lss_ttbar = mva_ttar;
			evtNtuple.MVA_2lss_ttV = mva_ttV;		
			evtNtuple.MT_met_lep0 = MVA_ttbar_vars.Get_MT_met_lep1();
			evtNtuple.n_jet25_recl = MVA_ttbar_vars.Get_nJet25();
			evtNtuple.mindr_lep0_jet = MVA_ttbar_vars.Get_mindr_lep1_jet();
			evtNtuple.mindr_lep1_jet = MVA_ttbar_vars.Get_mindr_lep2_jet();
			evtNtuple.avg_dr_jet = MVA_ttbar_vars.Get_avg_dr_jet();
			evtNtuple.max_lep_eta = MVA_ttbar_vars.Get_max_lep_eta();
			evtNtuple.lep0_conept = MVA_ttV_vars.Get_lep1_conePt();
			evtNtuple.lep1_conept = MVA_ttV_vars.Get_lep2_conePt();

			auto lep1type = local.leptons_fakeable[0].Type();
			auto lep2type = local.leptons_fakeable[1].Type();
			
			// weights and scale factors
			if (isdata) {
				local.weight = 1.;

				//////////////////////////
				// charge flip background (data driven)
				if (selection_type == Control_2los1tau) {
					float P1_misCharge =
						getEleChargeMisIDProb(local.leptons_fakeable[0],true);
					float P2_misCharge =
						getEleChargeMisIDProb(local.leptons_fakeable[1],true);

					local.weight = P1_misCharge + P2_misCharge;
				}

				//////////////////////////
				// fake lepton background (data driven)
				if (selection_type == Control_1lfakeable) {
					float f1 = getFakeRate(local.leptons_fakeable[0]);
					float f2 = getFakeRate(local.leptons_fakeable[1]);

					if (not local.leptons_fakeable[0].passTightSel()
						and local.leptons_fakeable[1].passTightSel())
						local.weight = f1/(1.-f1);
					else if (local.leptons_fakeable[0].passTightSel() and
							 not local.leptons_fakeable[1].passTightSel())
						local.weight = f2/(1.-f2);
					else if (not local.leptons_fakeable[0].passTightSel() and
							 not local.leptons_fakeable[1].passTightSel())
						local.weight = -f1*f2/((1.-f1)*(1.-f2));
				}
			}
			else {
				/// CSV weight
				//local.csv_weight = getEvtCSVWeight(local.jets_selected, "NA");
				local.csv_weight = getEvtCSVWeight(local.jets_selected, csv_iSys["NA"]);
				/// HLT sf
				if (lep1type == LeptonType::kmu and lep2type == LeptonType::kmu) {
					local.hlt_sf = 1.01;
				}
				else if (lep1type == LeptonType::kele and lep2type == LeptonType::kele) {
					local.hlt_sf = 1.02;
				}
				else {
					local.hlt_sf = 1.02;
				}
				
				/// total event weight
				local.weight =
					local.csv_weight * //local.gen_weight *
					local.hlt_sf * local.lepIDEff_sf;

				//if (local.weight < 0) {
				//	std::cout << "csv_weight : " << local.csv_weight << std::endl;
				//	std::cout << "gen_weight : " << local.gen_weight << std::endl;
				//	std::cout << "hlt_sf : " << local.hlt_sf << std::endl;
				//	std::cout << "lepIDEff_sf : " << local.lepIDEff_sf << std::endl;
				//}
			}

			// determine event category
			int ilep = -1;
			int ibtag = -1;
			
			if (lep1type==LeptonType::kmu and lep2type==LeptonType::kmu)
				ilep = 0;  // mumu
			else if (lep1type==LeptonType::kele and lep2type==LeptonType::kele)
				ilep = 1;  // ee
			else
				ilep = 2;  // emu

			if (local.jets_selected_btag_medium.size()>=2)
				ibtag = 1;  // b-tight
			else
				ibtag = 0;  // b-loose
			
			// 2D hist
			h_MVA_ttV_vs_ttbar[ilep][ibtag]
				->Fill(mva_ttar, mva_ttV, local.weight);
			
			// 1D shape
			int bin = partition2DBDT(mva_ttar, mva_ttV);
			
			h_MVA_shape[ilep][ibtag]
				->Fill(bin, local.weight);
			
			// systematics
			if (!isdata and doSystematics) {

				for (int isys = 0; isys < 16; ++isys) {
					double csv_weight_sys =
						//getEvtCSVWeight(local.jets_selected, sysList[isys]);
						getEvtCSVWeight(local.jets_selected, csv_iSys[sysList[isys]]);
					double evt_weight_sys =
						local.weight / local.csv_weight * csv_weight_sys;

					// 2D
					h_MVA_ttV_vs_ttbar_sys[ilep][ibtag][isys]
						->Fill(mva_ttar, mva_ttV,evt_weight_sys);
					// 1D shape
					h_MVA_shape_sys[ilep][ibtag][isys]
						->Fill(bin, evt_weight_sys);
				}
			}

		}
		
		if (produce_sync_ntuple or pass_event_selection) {
			// no event selection is applied for sync ntuples
			evtNtuple.write_ntuple(local);
			eventTree->Fill();
		}
					
	} // end of analysis_type == Analyze_2lss1tau
	
	if (analysis_type == Analyze_3l) {
		bool passHLT =
			local.pass_single_e or local.pass_single_mu or local.pass_double_mu or
			local.pass_double_e or local.pass_elemu;

		if (hltcut_off) passHLT = true;
		
		// Event selection
		bool pass_event_selection = 
			pass_event_sel_3l(local, selection_type) and passHLT;

		evtNtuple.initialize();

		if (pass_event_selection) {

			// weights and scale factors
			if (isdata) {
				local.weight = 1.;
			}
			else {
				/// CSV weight
				local.csv_weight =
					getEvtCSVWeight(local.jets_selected, csv_iSys["NA"]);
				    //getEvtCSVWeight(local.jets_selected, "NA");

				/// HLT sf: (1) +/- 0.06
				local.hlt_sf = 1.0;

				/// total event weight
				local.weight =
					local.csv_weight * 
					local.hlt_sf * local.lepIDEff_sf;
			}
			
		}

		if (produce_sync_ntuple or pass_event_selection) {
			evtNtuple.write_ntuple(local);
			eventTree->Fill();
		}
		
	}  // end of analysis_type == Analyze_3l
	
}

// ------------ method called once each job just before starting event loop
// ------------
void CU_ttH_EDA::beginJob()
{
	
	TH1::SetDefaultSumw2(true);
	TH2::SetDefaultSumw2(true);

	event_count = 0;
}

// ------------ method called once each job just after ending the event loop
// ------------
void CU_ttH_EDA::endJob() {
	
	if (analysis_type == Analyze_2lss1tau) {
		//std::cout << "Total number of samples: " << event_count << std::endl;
		
		if (not isdata and doLumiScale) {
			/*
			// Rescale histograms for MC
			h_MVA_ttV_vs_ttbar -> Scale(int_lumi * sample_xs / event_count);
			h_MVA_shape -> Scale(int_lumi * sample_xs / event_count);

			
			if (setup_sysHist) {
				for (auto h : h_MVA_ttV_vs_ttbar_sys) {
					h -> Scale(int_lumi * sample_xs / event_count);
				}
				
				for (auto h : h_MVA_shape_sys) {
					h -> Scale(int_lumi * sample_xs / event_count);
				}
			}
			*/
		}

	}

	return;
}

// ------------ method called when starting to processes a run  ------------
void CU_ttH_EDA::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup)
{
	/// Update HLTConfigProvider(s) for the new run
	bool hlt_config_changed = true; // init() updates this one
	if (!hlt_config.init(iRun, iSetup, hltTag, hlt_config_changed)) {
		std::cerr << "Warning, didn't find trigger process HLT,\t" << hltTag
				  << std::endl;

		return;
	}

	if (hlt_config_changed)
		std::cout << "New " << hltTag << " config has been loaded.\n";

	bool filter_config_changed = true; // init() updates this one
	if (!filter_config.init(iRun, iSetup, filterTag, filter_config_changed)) {
		std::cerr << "Warning, didn't find filter process HLT,\t" << filterTag
				  << std::endl;

		return;
	}

	if (filter_config_changed)
		std::cout << "New " << filterTag << " config has been loaded.\n";

	/// Set up filter and trigger name vectors and maps
	if (trigger_stats) {
		Set_up_trigger_name_vectors();

		if (Set_up_Run_histograms_triggers() != 0) {
			std::cerr << "Setting up histograms for trigger/filter counts has "
					  << "failed\n";

			return;
		}
	}
	
}

// ------------ method called when ending the processing of a run  ------------
void CU_ttH_EDA::endRun(const edm::Run &, const edm::EventSetup &)
{

	if (trigger_stats)
		End_Run_hist_fill_triggers();

	// report on triggers fired
	if (trigger_stats && dumpHLT_) {
		std::cout
			<< "***************************************************************"
			<< std::endl;
		std::cout << "  Summary for HLT: Total number of events = "
				  << event_count << std::endl;
		for (std::map<std::string, unsigned long>::const_iterator iter =
				 n_trigger_fired.begin();
			 iter != n_trigger_fired.end(); ++iter) {
			std::string name = iter->first;
			double eff = 100 * double(iter->second) / double(event_count);
			printf("\t %s,\t %lu,\t %.1f \n", name.c_str(), iter->second, eff);
		}
		std::cout
			<< "***************************************************************"
			<< std::endl;
		std::cout << "  Summary for Filters: Total number of events = "
				  << event_count << std::endl;
		for (std::map<std::string, unsigned long>::const_iterator iter =
				 n_filter_fired.begin();
			 iter != n_filter_fired.end(); ++iter) {
			std::string name = iter->first;
			double eff = 100 * double(iter->second) / double(event_count);
			printf("\t %s,\t %lu,\t %.1f \n", name.c_str(), iter->second, eff);
		}
		std::cout
			<< "***************************************************************"
			<< std::endl;
	}

	std::cout
		<< "***************************************************************"
		<< std::endl;
	std::cout << "  Total number of events = " << event_count << std::endl;
	std::cout
		<< "***************************************************************"
		<< std::endl;
}

// ------------ method called when starting to processes a luminosity block
// ------------
/*
void CU_ttH_EDA::beginLuminosityBlock(const edm::LuminosityBlock &,
const edm::EventSetup &)
{
}
*/

// ------------ method called when ending the processing of a luminosity block
// ------------
/*
void CU_ttH_EDA::endLuminosityBlock(const edm::LuminosityBlock & lumi, const edm::EventSetup & iConfig)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void CU_ttH_EDA::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	// The following says we do not know what parameters are allowed so do no
	// validation
	// Please change this to state exactly what you do use, even if it is no
	// parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

// define this as a CMSSW plugin
DEFINE_FWK_MODULE(CU_ttH_EDA);

#endif
