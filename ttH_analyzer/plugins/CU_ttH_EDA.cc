#ifndef CU_ttH_EDA_cc
#define CU_ttH_EDA_cc


// -*- C++ -*-
//
// Package:    ttH-LeptonPlusJets/CU_ttH_EDA
// Class:      CU_ttH_EDA
// 
/**\class CU_ttH_EDA CU_ttH_EDA.cc ttH-LeptonPlusJets/AnalysisCode/plugins/CU_ttH_EDA.cc

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
CU_ttH_EDA::CU_ttH_EDA(const edm::ParameterSet& iConfig)
{
	/*
	 * FIXME: an interface to a config reader is needed!!!
	 * 
	 * now do whatever initialization is needed
	*/
	verbose_ = false;
	dumpHLT_ = false;
	
	hltTag = "HLT";
	filterTag = "PAT";
	
	Set_up_tokens();
	Set_up_histograms();
	Set_up_output_files();
	
	
	/// temporary mock-up parameters
	std::string era = "2012_53x";
	MAODHelper_sample_nr = 2500;
	
	total_xs = 831.76;
	sample_n = 25446993;
	int_lumi = 10000;
	weight_sample = int_lumi * total_xs / sample_n;
	
	
	miniAODhelper.SetUp(era,
		MAODHelper_sample_nr,
		analysisType::LJ,	// LJ, DIL, TauLJ, TauDIL
		false	// is data
	);
	
	
	/// Setting up cuts
	min_tight_lepton_pT = 20;
}



/// Destructor
CU_ttH_EDA::~CU_ttH_EDA()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
	Close_output_files();
}



// ------------ method called for each event  ------------
void CU_ttH_EDA::analyze(const edm::Event &iEvent,
	const edm::EventSetup &iSetup)
{
	using namespace edm;
	++event_count;
	
	/// Declaring local struct for data readout and manipulations
	CU_ttH_EDA_event_vars local;
	
	
	/// Triggers have not fired yet. Check_triggers, Check_filters will adjust
	local.pass_single_e = false;
	local.pass_single_mu = false;
        local.pass_double_mu = false;
        local.pass_double_e = false;
        local.pass_elemu = false;
	Update_common_vars(iEvent, local);
	
	
	/// Create and set up edm:Handles in stack mem.
	edm_Handles handle;
	Set_up_handles(iEvent, handle, token);
	
	
	/// Run checks on event containers via their handles
	Check_triggers(handle.triggerResults, local);
	Check_filters(handle.filterResults);
	Check_vertices_set_MAODhelper(handle.vertices);
// 	Check_beam_spot(BS);	// dumb implementation
	
	
	/// Get and set miniAODhelper's jet corrector from the event setup
	miniAODhelper.SetJetCorrector(
		JetCorrector::getJetCorrector("ak4PFchsL1L2L3", iSetup));
	
	
	// 	weight_gen = event_gen_info.product()->weight();
	local.weight = weight_sample * (handle.event_gen_info.product()->weight());
	
	h_hlt->Fill(0., 1);
	h_flt->Fill(0., 1);
	
	
	
	/// Lepton selection
	local.e_selected =
		miniAODhelper.GetSelectedElectrons(*(handle.electrons),
			min_tight_lepton_pT, electronID::electronPhys14M);
	local.mu_selected =
		miniAODhelper.GetSelectedMuons(*(handle.muons),
			min_tight_lepton_pT, muonID::muonTight);
	
	local.n_electrons = static_cast<int>(local.e_selected.size());
	local.n_muons = static_cast<int>(local.mu_selected.size());
	
	
	/// Sort leptons by pT
	local.mu_selected_sorted = miniAODhelper.GetSortedByPt(local.mu_selected);
	local.e_selected_sorted = miniAODhelper.GetSortedByPt(local.e_selected);
	
	
	
	/// Jet selection
	local.jets_raw = miniAODhelper.GetUncorrectedJets(handle.jets);
	local.jets_no_mu = miniAODhelper.RemoveOverlaps(local.mu_selected,
		local.jets_raw);
	local.jets_no_mu_e = miniAODhelper.RemoveOverlaps(local.e_selected,
		local.jets_no_mu);
	local.jets_corrected = miniAODhelper.GetCorrectedJets(local.jets_no_mu_e,
		iEvent, iSetup);
	local.jets_selected = miniAODhelper.GetSelectedJets(local.jets_corrected,
		25., 2.4, jetID::jetLoose, '-');
	local.jets_selected_tag = miniAODhelper.GetSelectedJets(local.jets_corrected,
		25., 2.4, jetID::jetLoose, 'M');
	
	local.n_jets = static_cast<int>(local.jets_selected.size());
	local.n_btags = static_cast<int>(local.jets_selected_tag.size());
	
	
	/// Sort jets by pT
	local.jets_selected_sorted = miniAODhelper.GetSortedByPt(
		local.jets_selected);
	local.jets_selected_tag_sorted = miniAODhelper.GetSortedByPt(
		local.jets_selected_tag);
	
	
	
	/// Top and Higgs tagging using collections through handles. adjusts local.<tag>
	Top_tagger(handle.top_jets, local);
	Higgs_tagger(handle.subfilter_jets, local);	// has some weird memory features on Top_tagger
	
	
	/// Get Corrected MET, !!!not yet used!!!
	//may need to be placed in CU_ttH_EDA_event_vars
	local.MET_corrected = handle.METs->front();	//miniAODhelper.GetCorrectedMET( METs.at(0), pfJets_forMET, iSysType );
	
	
	/// Check tags, fill hists, print events
	Check_Fill_Print_ej(local);
	Check_Fill_Print_muj(local);
        Check_Fill_Print_dimuj(local);
        Check_Fill_Print_dielej(local);
        Check_Fill_Print_elemuj(local);
}



// ------------ method called once each job just before starting event loop  ------------
void CU_ttH_EDA::beginJob()
{
	TH1::SetDefaultSumw2(true);
	
	event_count = 0;
}



// ------------ method called once each job just after ending the event loop  ------------
void CU_ttH_EDA::endJob() 
{
	return;
}



// ------------ method called when starting to processes a run  ------------
void CU_ttH_EDA::beginRun(edm::Run const &iRun,
	edm::EventSetup const &iSetup)
{
	/// Update HLTConfigProvider(s) for the new run
	bool hlt_config_changed = true;	// init() updates this one
	if (!hlt_config.init(iRun, iSetup, hltTag, hlt_config_changed)) {
		std::cerr << "Warning, didn't find trigger process HLT,\t" <<
			hltTag << std::endl;
		
		return;
	}
	
	if (hlt_config_changed)
		std::cout << "New " << hltTag << " config has been loaded.\n";
	
	
	bool filter_config_changed = true;	// init() updates this one
	if (!filter_config.init(iRun, iSetup, filterTag, filter_config_changed)) {
		std::cerr << "Warning, didn't find filter process HLT,\t" <<
			filterTag << std::endl;
		
		return;
	}
	
	if (filter_config_changed)
		std::cout << "New " << filterTag << " config has been loaded.\n";
	
	
	/// Set up filter and trigger name vectors and maps
	Set_up_name_vectors();
	
	
	if (Set_up_Run_histograms() != 0) {
		std::cerr << "Setting up histograms for trigger/filter counts has " <<
			"failed\n";
		
		return;
	}
}



// ------------ method called when ending the processing of a run  ------------

void CU_ttH_EDA::endRun(edm::Run const&, edm::EventSetup const&)
{
	// report results of sync exercises
	std::cout << "***************************************************************" << std::endl;
	std::cout << "\t Synchronization for mu" << std::endl;
	std::cout << "Selection \t Number of events\n";
	for (int i = 0; i < h_tth_syncex1_mu->GetNbinsX(); ++i)
		printf("%s\t %.0f\n",
			h_tth_syncex1_mu->GetXaxis()->GetBinLabel(i + 1),
			h_tth_syncex1_mu->GetBinContent(i + 1));
	
	std::cout << "***************************************************************" << std::endl;
	std::cout << "\t Synchronization for e" << std::endl;
	std::cout << "Selection \t Number of events\n";
	for (int i = 0; i < h_tth_syncex1_ele->GetNbinsX(); ++i)
		printf("%s\t %.0f\n",
			h_tth_syncex1_ele->GetXaxis()->GetBinLabel(i + 1),
			h_tth_syncex1_ele->GetBinContent(i + 1));
	
	std::cout << "***************************************************************" << std::endl;
	std::cout << "\t Synchronization for di-mu" << std::endl;
	std::cout << "Selection \t Number of events\n";
	for (int i = 0; i < h_tth_syncex1_dimu->GetNbinsX(); ++i)
		printf("%s\t %.0f\n",
			h_tth_syncex1_dimu->GetXaxis()->GetBinLabel(i + 1),
			h_tth_syncex1_dimu->GetBinContent(i + 1));
	
	std::cout << "***************************************************************" << std::endl;
	
	std::cout << "\t Synchronization for di-ele" << std::endl;
	std::cout << "Selection \t Number of events\n";
	for (int i = 0; i < h_tth_syncex1_diele->GetNbinsX(); ++i)
		printf("%s\t %.0f\n",
			h_tth_syncex1_diele->GetXaxis()->GetBinLabel(i + 1),
			h_tth_syncex1_diele->GetBinContent(i + 1));
	
	std::cout << "***************************************************************" << std::endl;
		
	std::cout << "\t Synchronization for ele-mu" << std::endl;
	std::cout << "Selection \t Number of events\n";
	for (int i = 0; i < h_tth_syncex1_elemu->GetNbinsX(); ++i)
		printf("%s\t %.0f\n",
			h_tth_syncex1_elemu->GetXaxis()->GetBinLabel(i + 1),
			h_tth_syncex1_elemu->GetBinContent(i + 1));
	
	std::cout << "***************************************************************" << std::endl;

        
        End_Run_hist_fill();
	
	
	
	// report on triggers fired
	if (dumpHLT_) {
		std::cout << "***************************************************************" << std::endl;
		std::cout << "  Summary for HLT: Total number of events = " << event_count << std::endl;
		for (std::map<std::string, unsigned long>::const_iterator iter = n_trigger_fired.begin();
			iter != n_trigger_fired.end(); ++iter) {
			std::string name = iter->first;
			double eff = 100 * double(iter->second) / double(event_count);
			printf("\t %s,\t %lu,\t %.1f \n", name.c_str(), iter->second, eff);
		}
		std::cout << "***************************************************************" << std::endl;
		std::cout << "  Summary for Filters: Total number of events = " << event_count << std::endl;
		for (std::map<std::string, unsigned long>::const_iterator iter = n_filter_fired.begin();
			iter != n_filter_fired.end(); ++iter) {
			std::string name = iter->first;
			double eff = 100 * double(iter->second) / double(event_count);
			printf("\t %s,\t %lu,\t %.1f \n", name.c_str(), iter->second, eff);
		}
		std::cout << "***************************************************************" << std::endl;
	}

	std::cout << "***************************************************************" << std::endl;
	std::cout << "  Total number of events = " << event_count << std::endl;
	std::cout << "***************************************************************" << std::endl;
}



// ------------ method called when starting to processes a luminosity block  ------------
/*
void CU_ttH_EDA::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/



// ------------ method called when ending the processing of a luminosity block  ------------
/*
void CU_ttH_EDA::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void CU_ttH_EDA::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}



//define this as a CMSSW plugin
DEFINE_FWK_MODULE(CU_ttH_EDA);



#endif
