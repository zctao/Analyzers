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
	// Turn off event selection
	turn_off_event_sel (iConfig.getParameter<bool>("turn_off_event_sel")),
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
	// JEC
	JECType (iConfig.getParameter<string>("JECType")),
	//jet_corrector (iConfig.getParameter<string>("jet_corrector")),
	isdata (iConfig.getParameter<bool>("using_real_data")),
	selection_region (iConfig.getParameter<string>("selection_region")),
	// debug flag
	debug (iConfig.getParameter<bool>("debug_mode"))
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

	if (not isdata) {
		Set_up_LeptonSF_Lut();
	}

}

/// Destructor
CU_ttH_EDA::~CU_ttH_EDA()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
	//delete BTagCaliReader;

	if (not isdata) {
		delete file_recoToLoose_leptonSF_mu1_b;
		delete file_recoToLoose_leptonSF_mu1_e;
		delete file_recoToLoose_leptonSF_mu2;
		delete file_recoToLoose_leptonSF_mu3;
		
		delete file_recoToLoose_leptonSF_el;
		delete file_recoToLoose_leptonSF_gsf;
		
		if (analysis_type == Analyze_2lss1tau) {
			delete file_looseToTight_leptonSF_mu_2lss;
			delete file_looseToTight_leptonSF_el_2lss;
		}
		if (analysis_type == Analyze_3l) {
			delete file_looseToTight_leptonSF_mu_3l;
			delete file_looseToTight_leptonSF_el_3l;
		}
	}
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
	local.mc_weight = 1.;
	local.mc_weight_scale_muF0p5 = 1.;
	local.mc_weight_scale_muF2 = 1.;
	local.mc_weight_scale_muR0p5 = 1.;
	local.mc_weight_scale_muR2 = 1.;
	local.hlt_sf = 1.;
	local.lepIDEff_sf = 1.;
	local.pu_weight = 1.;

	local.npuTrue = -1.;
	local.npuInTime = -1.;
	
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

	// JEC
	edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
	iSetup.get<JetCorrectionsRecord>().get("AK4PF",JetCorParColl);
	const JetCorrectorParameters & JetCorPar = (*JetCorParColl)["Uncertainty"];
	miniAODhelper.SetJetCorrectorUncertainty(JetCorPar);

	
	if (not isdata) {
		// pileup
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		
		for (PVI = handle.PU_info->begin(); PVI != handle.PU_info->end();
			 ++PVI) {
			int BX = PVI->getBunchCrossing();
			if (BX == 0) { // "0" is the in-time crossing, negative values are the early crossings, positive are late
				local.npuTrue = PVI -> getTrueNumInteractions();
				local.npuInTime = PVI -> getPU_NumInteractions();
				break;
			}
		}
		
		// MC weights
		double genWeight = handle.event_gen_info.product()->weight();
		
		local.mc_weight = genWeight / abs(genWeight);

		genWeightSum += genWeight / abs(genWeight);

		if (handle.event_lhe_info.isValid()) {
			if (handle.event_lhe_info->weights().size() > 6) {
				local.mc_weight_scale_muF0p5 = // muF = 0.5 | muR = 1
					local.mc_weight * (handle.event_lhe_info->weights()[2].wgt)/
					(handle.event_lhe_info->originalXWGTUP()); 
				local.mc_weight_scale_muF2 = // muF = 2 | muR = 1
					local.mc_weight * (handle.event_lhe_info->weights()[1].wgt)/
					(handle.event_lhe_info->originalXWGTUP()); 
				local.mc_weight_scale_muR0p5 = // muF = 1 | muR = 0.5
					local.mc_weight * (handle.event_lhe_info->weights()[6].wgt)/
					(handle.event_lhe_info->originalXWGTUP()); 
				local.mc_weight_scale_muR2 = // muF = 1 | muR = 2
					local.mc_weight * (handle.event_lhe_info->weights()[3].wgt)/
					(handle.event_lhe_info->originalXWGTUP());
			}
		}
	}

	if (trigger_stats) {
		h_hlt->Fill(0., 1);
		h_flt->Fill(0., 1);
	}

	/*
	  Muons
	*/
	for (const auto& mu : *(handle.muons)){
		if (mu.userFloat("idPreselection") > 0.5 and
			mu.pt() > 5. and abs(mu.eta()) < 2.4)
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
		if (ele.userFloat("idPreselection") > 0.5 and
			ele.pt() > 7. and abs(ele.eta()) < 2.5)
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
		if (tau.userFloat("idPreselection")>0.5 and
			tau.pt() > 20. and abs(tau.eta()) < 2.3)
			local.tau_preselected.push_back(tau);  // for cleaning
	}

	// remove overlap with muons and electrons
	local.tau_preselected =
		removeOverlapdR(local.tau_preselected, local.e_preselected, 0.3);
	local.tau_preselected =
		removeOverlapdR(local.tau_preselected, local.mu_preselected, 0.3);
	
	// sort by pT
	local.tau_preselected_sorted =
		miniAODhelper.GetSortedByPt(local.tau_preselected);
	
	for (const auto& tau : local.tau_preselected_sorted) {
		if (tau.userFloat("idSelection")>0.5)
			local.tau_selected.push_back(tau);
	}
	
	// number of selected leptons
	local.n_electrons_loose = static_cast<int>(local.e_preselected.size());
	local.n_electrons_fakeable = static_cast<int>(local.e_fakeable.size());
	local.n_electrons_tight = static_cast<int>(local.e_tight.size());
	local.n_muons_loose = static_cast<int>(local.mu_preselected.size());
	local.n_muons_fakeable = static_cast<int>(local.mu_fakeable.size());
	local.n_muons_tight = static_cast<int>(local.mu_tight.size());
	local.n_taus_pre = static_cast<int>(local.tau_preselected.size());
	local.n_taus = static_cast<int>(local.tau_selected.size());

	if (debug) {
		std::cout << "************************************" << std::endl;
		std::cout << "run:ls:event : " << local.run_nr <<":"<<local.lumisection_nr<<":"<<local.event_nr<<std::endl;
		std::cout << "lepton selection : " << std::endl;
		std::cout << "n_muon_fakeable : " << local.n_muons_fakeable << std::endl;
		std::cout << "n_muon_tight : " << local.n_muons_tight << std::endl;
		std::cout << "fakeable muons : " << std::endl;
		for (auto & mu : local.mu_fakeable) {
			std::cout<<"pt eta phi : " << mu.pt()<<" "<<mu.eta()<<" "<<mu.phi()<<std::endl;
		}
		std::cout << "n_electrons_fakeable : " << local.n_electrons_fakeable << std::endl;	
		std::cout << "n_electrons_tight : " << local.n_electrons_tight << std::endl;
		std::cout << "fakeable electrons : " << std::endl;
		for (auto & e : local.e_fakeable) {
			std::cout<<"pt eta phi : " << e.pt()<<" "<<e.eta()<<" "<<e.phi()<<std::endl;
		}
		std::cout << "n_taus_pre : " << local.n_taus_pre << std::endl;
		std::cout << "n_taus : " << local.n_taus << std::endl;
		std::cout << "preselected taus : " << std::endl;
		for (auto & tau : local.tau_preselected) {
			std::cout <<"pt eta phi : "<< tau.pt()<<" "<<tau.eta()<<" "<<tau.phi()<<std::endl;
		}
	}
	
	/*
	  Jets selection
	*/
	if (isdata) { // no smearing or systematics estimation with data sample
		local.jets_raw = *(handle.jets);
	}
	else {
		local.jets_raw =
			miniAODhelper.GetCorrectedJets(*(handle.jets), iEvent, iSetup,
									   handle.genJets, JECTypes[JECType]);
	}
	
	// overlap removal by dR
	std::vector<pat::Jet> jets_no_mu = removeOverlapdR(local.jets_raw, local.mu_fakeable, 0.4);
	std::vector<pat::Jet> jets_no_mu_e = removeOverlapdR(jets_no_mu, local.e_fakeable, 0.4);
	std::vector<pat::Jet> jets_cleaned = removeOverlapdR(jets_no_mu_e, local.tau_preselected, 0.4);

	// selected jets
	local.jets_selected = miniAODhelper.GetSelectedJets(
	    jets_cleaned, 25., 2.4, jetID::jetLoose, '-');

	if (debug) {
		std::cout << "all slimmedJets : " << std::endl;
		for (const auto & j : *(handle.jets)) {
			std::cout << "pt eta phi E : " << j.pt()<<" "<<j.eta()<<" "<<j.phi()<<" "<<j.energy()<<std::endl;
		}
		std::cout << "after JEC : " << std::endl;
		for (const auto & j : local.jets_raw) {
			std::cout << "pt eta phi E : " << j.pt()<<" "<<j.eta()<<" "<<j.phi()<<" "<<j.energy()<<std::endl;
		}
		std::cout << "after overlap removal : " << std::endl;
		for (const auto & j : jets_cleaned) {
			std::cout << "pt eta phi E CSV : " << j.pt()<<" "<<j.eta()<<" "<<j.phi()<<" "<<j.energy()<<" "<<j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<<std::endl;
		}
		std::cout << "selected jets : " << std::endl;
		for (const auto & j : local.jets_selected) {
			std::cout << "pt eta phi E CSV : " << j.pt()<<" "<<j.eta()<<" "<<j.phi()<<" "<<j.energy()<<" "<<j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<<std::endl;
		}
	}
	
	// sort by pT
	local.jets_selected_sorted =
		miniAODhelper.GetSortedByPt(local.jets_selected);

	/*
	  BTag
	*/
	local.jets_selected_btag_loose = miniAODhelper.GetSelectedJets(
		local.jets_selected_sorted, 20., 2.5, jetID::jetLoose,
		'L');
	local.jets_selected_btag_medium = miniAODhelper.GetSelectedJets(
	    local.jets_selected_btag_loose,20., 2.5,
		jetID::jetLoose, 'M');

	local.n_jets = static_cast<int>(local.jets_selected.size());
	local.n_btags_loose = static_cast<int>(local.jets_selected_btag_loose.size());
	local.n_btags_medium = static_cast<int>(local.jets_selected_btag_medium.size());

	if (debug) {
		std::cout << "n_jets : " << local.n_jets << std::endl;
		std::cout << "n_btags_loose : " << local.n_btags_loose << std::endl;
		std::cout << "n_btags_medium : " << local.n_btags_medium << std::endl;
	}
	
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

	if (debug) {
		std::cout << "MET : " << std::endl;
		std::cout << "pt : " << met << std::endl;
		std::cout << "mht : " << mht << std::endl;
	}
	
	
	if (analysis_type == Analyze_2lss1tau) {
		
		// event category index, value set by pass_event_sel_2l()
		int ilep = -1;
		int ibtag = -1;
		
		// Event selection	
		bool pass_event_selection =
			pass_event_sel_2l(local, selection_type, ilep, ibtag);
		
		// match HLT path by category
		bool matchHLT = false;
		if (ilep == 0) // mumu
			matchHLT = local.pass_single_mu or local.pass_double_mu;
		else if (ilep == 1) // ee
			matchHLT = local.pass_single_e or local.pass_double_e;
		else if (ilep == 2) // emu
			matchHLT = local.pass_single_e or local.pass_single_mu or
				local.pass_elemu;

		pass_event_selection = pass_event_selection and (matchHLT or hltcut_off);
		
		// Check if all HLTs failed for debug purpose: assert(not matchHLT);
		
		double mva_ttbar = -9999.;
		double mva_ttV = -9999.;
		
		if (turn_off_event_sel or pass_event_selection) {
			
			// Calculate event-level MVA variables
			MVA_ttbar_vars.Calculate_mvaVars(local.leptons_fakeable,
											 local.tau_selected,
											 local.jets_selected_sorted,
											 local.pfMET);
			MVA_ttV_vars.Calculate_mvaVars(local.leptons_fakeable,
										   local.tau_selected,
										   local.jets_selected_sorted,
										   local.pfMET);
			
			mva_ttbar = MVA_ttbar_vars.Get_mvaScore();
			mva_ttV = MVA_ttV_vars.Get_mvaScore();
		
		}	

		if (pass_event_selection) {
			assert(local.leptons_fakeable.size() >= 2);
			
			// Weights and scale factors
			//////////////////////////
			// signal region
			if (selection_type == Signal_2lss1tau) {
				if (isdata) {
					local.weight = 1.;
				}
				else {
					/// CSV weight
					//local.csv_weight = getEvtCSVWeight(local.jets_selected, "NA");
					local.csv_weight = getEvtCSVWeight(local.jets_selected, csv_iSys["NA"]);
					/// HLT sf
					local.hlt_sf = getLepHLTSF(ilep);

					// LeptonSF
					local.lepIDEff_sf *= getLeptonSF(local.leptons_fakeable[0]);
					local.lepIDEff_sf *= getLeptonSF(local.leptons_fakeable[1]);
					
					/// total event weight
					local.weight =
						local.csv_weight * local.mc_weight *
						local.hlt_sf * local.lepIDEff_sf;
				}
				
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

			//////////////////////////
			// charge flip background (data driven)
			if (selection_type == Control_2los1tau) {
				float P1_misCharge =
					getEleChargeMisIDProb(local.leptons_fakeable[0],true);
				float P2_misCharge =
					getEleChargeMisIDProb(local.leptons_fakeable[1],true);
				
				local.weight = P1_misCharge + P2_misCharge;
			}

			// Fill histograms
			// 2D hist
			h_MVA_ttV_vs_ttbar[ilep][ibtag]
				->Fill(mva_ttbar, mva_ttV, local.weight);
			
			// 1D shape
			int bin = partition2DBDT(mva_ttbar, mva_ttV);
			
			h_MVA_shape[ilep][ibtag]
				->Fill(bin, local.weight);
			
			// systematics
			if (!isdata and doSystematics) {

				// CSV reweight
				for (int isys = 0; isys < 16; ++isys) {
					double csv_weight_sys =
						//getEvtCSVWeight(local.jets_selected, sysList[isys]);
						getEvtCSVWeight(local.jets_selected, csv_iSys[sysList[isys]]);
					double evt_weight_sys = // TODO: protect control region on MC
						local.weight / local.csv_weight * csv_weight_sys;

					// 2D
					//h_MVA_ttV_vs_ttbar_sys[ilep][ibtag][isys]
					//	->Fill(mva_ttbar, mva_ttV,evt_weight_sys);
					// 1D shape
					h_MVA_shape_csv_sys[ilep][ibtag][isys]
						->Fill(bin, evt_weight_sys);
				}

				// theoretical uncertainty from renormalization and
				// factorization scale
				h_MVA_shape_thu_sys[ilep][ibtag][0]->
					Fill(bin,local.mc_weight_scale_muF0p5);
				h_MVA_shape_thu_sys[ilep][ibtag][1]->
					Fill(bin,local.mc_weight_scale_muF2);
				h_MVA_shape_thu_sys[ilep][ibtag][2]->
					Fill(bin,local.mc_weight_scale_muR0p5);
				h_MVA_shape_thu_sys[ilep][ibtag][3]->
					Fill(bin,local.mc_weight_scale_muR2);
			}
			
		}

		// Write Ntuple
		if (turn_off_event_sel or pass_event_selection) {

			evtNtuple.initialize();
			
			evtNtuple.write_ntuple(local);
			
			// Write MVA variables to ntuple
			evtNtuple.MVA_2lss_ttbar = mva_ttbar;
			evtNtuple.MVA_2lss_ttV = mva_ttV;		
			evtNtuple.MT_met_lep0 = MVA_ttbar_vars.Get_MT_met_lep1();
			evtNtuple.n_jet25_recl = MVA_ttbar_vars.Get_nJet25();
			evtNtuple.mindr_lep0_jet = MVA_ttbar_vars.Get_mindr_lep1_jet();
			evtNtuple.mindr_lep1_jet = MVA_ttbar_vars.Get_mindr_lep2_jet();
			evtNtuple.avg_dr_jet = MVA_ttbar_vars.Get_avg_dr_jet();
			evtNtuple.max_lep_eta = MVA_ttbar_vars.Get_max_lep_eta();
			evtNtuple.lep0_conept = MVA_ttV_vars.Get_lep1_conePt();
			evtNtuple.lep1_conept = MVA_ttV_vars.Get_lep2_conePt();
			
			eventTree->Fill();
		}
					
	} // end of analysis_type == Analyze_2lss1tau
	
	if (analysis_type == Analyze_3l) {
		bool matchHLT = true;  // FIXME

		if (hltcut_off) matchHLT = true;
		
		// Event selection
		bool pass_event_selection = 
			pass_event_sel_3l(local, selection_type) and matchHLT;

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

		if (turn_off_event_sel or pass_event_selection) {
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
	genWeightSum = 0.;
}

// ------------ method called once each job just after ending the event loop
// ------------
void CU_ttH_EDA::endJob() {

	// fill GenWeightSum histogram
	h_SumGenWeight -> SetBinContent(1,genWeightSum);
	
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

	if (hlt_config_changed) {
		std::cout << "New " << hltTag << " config has been loaded.\n";
		
		if (not hltcut_off) Set_up_HLT_path_name();
	}

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
