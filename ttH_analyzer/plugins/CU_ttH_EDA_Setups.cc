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
	h_GenWeightProcessed = fs_->make<TH1D>("h_GenWeightProcessed","",3,-1.5,1.5);
	h_GenWeightxPUProcessed = fs_->make<TH1D>("h_GenWeightxPUProcessed","",1,0.5,1.5);
	
	// categories
	TString lep_cat[3] = {"mumu", "ee", "emu"};
	TString btag_cat[2] = {"bloose", "bmedium"};
	
	if (analysis_type == Analyze_2lss1tau) {

	}

	if (analysis_type == Analyze_3l) {
		/*
		h_mTWl = fs_->make<TH1D>("h_mTWl","", 10, 0, 300);
		h_met = fs_->make<TH1D>("h_met", "", 8, 0, 200);
		h_nJets = fs_->make<TH1D>("h_nJets", "", 8, 0, 8);
		h_lep_charge = fs_->make<TH1D>("h_lep_charge", "", 10, -5, 5);
		h_mZ = fs_->make<TH1D>("h_mZ", "", 15, 60, 120);
		*/
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

#endif
	
