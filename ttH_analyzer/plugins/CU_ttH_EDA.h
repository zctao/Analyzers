#ifndef CU_ttH_EDA_h
#define CU_ttH_EDA_h

/// Core libraries
#include <memory>
#include <cstdio>	// printf, fprintf
#include <stdexcept> // standard exceptions
#include <cmath>

/// CMSSW user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1Trigger/interface/L1HFRingsFwd.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"

/// ROOT includes
#include "TH1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TParameter.h"
#include "TLorentzVector.h"

/// Higgs and top tagger
#include "MiniAOD/BoostedObjects/interface/HTTTopJet.h"
#include "MiniAOD/BoostedObjects/interface/SubFilterJet.h"
#include "BoostedTTH/BoostedAnalyzer/interface/BoostedUtils.hpp"
#include "MiniAOD/MiniAODHelper/interface/TopTagger.h"
#include "MiniAOD/MiniAODHelper/interface/HiggsTagger.h"

/// structs for holding multiple edm::Handle and EDGetTokenT
#include "CU_ttH_EDA_Handles.h"

/// Configuration reader
#include "yaml-cpp/yaml.h"

/*
 *
 * enum for analysis type.
 * Purpose: allows faster in-function comparisons
 *
 */

enum analysis_types {
	Analyze_lepton_jet,
	Analyze_dilepton,
	Analyze_taus_lepton_jet,
	Analyze_taus_dilepton
};

/*
 *
 * struct for per-event variables used in analyze(...)
 *
 */

struct CU_ttH_EDA_event_vars {
	double weight; // total event weight (there are partial weights)

	/// Common, run parameters
	int run_nr;
	int event_nr;
	int lumisection_nr;

	/// Number of tags per event
	int n_electrons;
	int n_muons;
	int n_noniso_taus;
	int n_loose_taus;
	int n_medium_taus;
	int n_tight_taus;
	int n_jets;
	//int n_btags;
	int n_btags_loose;
	int n_btags_medium;
	int n_ttags;
	int n_Htags;

	/// Passing-trigger flags
	bool pass_single_e;
	bool pass_single_mu;
	bool pass_double_mu;
	bool pass_double_e;
	bool pass_elemu;
	
	/// Particle container vectors
	std::vector<pat::Electron> e_selected;
	std::vector<pat::Electron> e_selected_sorted;
	std::vector<pat::Muon> mu_selected;
	std::vector<pat::Muon> mu_selected_sorted;

	std::vector<pat::Tau> noniso_tau_selected;
	std::vector<pat::Tau> noniso_tau_selected_sorted;
	std::vector<pat::Tau> loose_tau_selected;
	std::vector<pat::Tau> loose_tau_selected_sorted;
	std::vector<pat::Tau> medium_tau_selected;
	std::vector<pat::Tau> medium_tau_selected_sorted;
	std::vector<pat::Tau> tight_tau_selected;
	std::vector<pat::Tau> tight_tau_selected_sorted;
	
	std::vector<pat::Jet> jets_raw;
	std::vector<pat::Jet> jets_no_mu;
	std::vector<pat::Jet> jets_no_mu_e;
	std::vector<pat::Jet> jets_corrected;
	std::vector<pat::Jet> jets_selected;
	std::vector<pat::Jet> jets_selected_sorted;
	//std::vector<pat::Jet> jets_selected_tag;
	//std::vector<pat::Jet> jets_selected_tag_sorted;
	std::vector<pat::Jet> jets_selected_tag_loose;
	std::vector<pat::Jet> jets_selected_tag_loose_sorted;
	std::vector<pat::Jet> jets_selected_tag_medium;
	std::vector<pat::Jet> jets_selected_tag_medium_sorted;

	double weight_gen;
	
	/// Other quantities
	pat::MET MET_corrected;
	
	double dimuon_mass;
	double dielectron_mass;
	double dilepton_mass;
	double ditau_mass;
};

struct CU_ttH_EDA_gen_vars {
	/// Generated particle container vectors
	std::vector<reco::GenParticle> x; // mediators
	std::vector<reco::GenParticle> tops;
	reco::CandidateCollection x_daughters; // or edm::OwnVector<reco::Candidate>
	reco::CandidateCollection top_daughters;
	reco::CandidateCollection w_daughters;

	std::vector<int> tau_class;

	double ditop_mass;
	double ditau_mass;
};

 /*
 *
 * Analyzer class
 *
 */

class CU_ttH_EDA : public edm::EDAnalyzer
{
  public:
	explicit CU_ttH_EDA(const edm::ParameterSet &);
	~CU_ttH_EDA();

	static void fillDescriptions(edm::ConfigurationDescriptions &);

  private:
	/*
	 * Function/method section
	*/

	/// Standard EDAnalyzer functions
	void analyze(const edm::Event &, const edm::EventSetup &) override;

	void beginJob() override;
	void endJob() override;

	void beginRun(const edm::Run &, const edm::EventSetup &) override;
	void endRun(const edm::Run &, const edm::EventSetup &) override;

	// 	virtual void beginLuminosityBlock(edm::LuminosityBlock const&,
	// 		edm::EventSetup const&) override;
	// 	virtual void endLuminosityBlock(edm::LuminosityBlock const&,
	// 		edm::EventSetup const&) override;

	/// One-time-run functions
	void Close_output_files();		 // at ~CU_ttH_EDA()
	void Load_configuration(string); // at CU_ttH_EDA(), runs _MAODH()
	void Load_configuration_set_type(const string &); // sets analysis_type
	void Load_configuration_MAODH(bool); // runs miniAODhelper.SetUp
	void Set_up_histograms();			 // at CU_ttH_EDA()
	void Set_up_output_files();			 // at CU_ttH_EDA()
	void Set_up_tokens();				 // at CU_ttH_EDA()
	
	void Setup_Tree();

	int Set_up_Run_histograms_triggers(); // at beginRun(), after
										  // Set_up_name_vectors()

	void Set_up_trigger_name_vectors(); // at beginRun()

	int End_Run_hist_fill_triggers(); // fill histograms at endRun()

	/// Per-event functions
	void Update_common_vars(const edm::Event &, CU_ttH_EDA_event_vars &);

	/// Object checks. Returns 1 in case of an error
	int Check_beam_spot(edm::Handle<reco::BeamSpot>);
	int Check_triggers(edm::Handle<edm::TriggerResults>,
					   CU_ttH_EDA_event_vars &); // adjusts event variables
	int Check_filters(edm::Handle<edm::TriggerResults>);
	int Check_vertices_set_MAODhelper(edm::Handle<reco::VertexCollection>, reco::Vertex &);

	// trigger iterator, part of Check_triggers()
	bool Check_triggers_iterator(const vector<string> &,
								 edm::Handle<edm::TriggerResults>);

	/// Taggers. Returns 1 in case of an error
	int Higgs_tagger(Handle<boosted::SubFilterJetCollection>,
					 CU_ttH_EDA_event_vars &); // FIXME: uses b-tag medium WP
	int Top_tagger(Handle<boosted::HTTTopJetCollection>,
				   CU_ttH_EDA_event_vars &);
	TopTagger toptagger;

	/// Other functions
	void Check_Fill_Print_ej(CU_ttH_EDA_event_vars &);
	void Check_Fill_Print_muj(CU_ttH_EDA_event_vars &);
	void Check_Fill_Print_dimuj(CU_ttH_EDA_event_vars &);
	void Check_Fill_Print_dielej(CU_ttH_EDA_event_vars &);
	void Check_Fill_Print_elemuj(CU_ttH_EDA_event_vars &);
	//void Check_Fill_Print_dileptauh(CU_ttH_EDA_event_vars &);
	//void Check_Fill_Print_eleditauh(CU_ttH_EDA_event_vars &);
	//void Check_Fill_Print_muditauh(CU_ttH_EDA_event_vars &);
	//bool pass_cut(CU_ttH_EDA_event_vars &, string);
	//bool pass_multi_cuts(CU_ttH_EDA_event_vars &, std::vector<string>, bool, TH1D*, int);
	void Make_Ntuple(CU_ttH_EDA_gen_vars &, CU_ttH_EDA_event_vars &, TTree *);
	void Fill_Tau_Eff_Hist(CU_ttH_EDA_gen_vars &, CU_ttH_EDA_event_vars &);

	template <typename T1, typename T2>
		std::vector<T1>
		removeOverlap(const std::vector<T1>& v1, const std::vector<T2>& v2, double dR = 0.02);
	
	/// Gen information functions
	const reco::GenParticle* getGenTau(const pat::Tau &);
	void Get_GenInfo(Handle<reco::GenParticleCollection>,
					 Handle<pat::PackedGenParticleCollection>,
					 CU_ttH_EDA_gen_vars &);
	void printDecayChain(const reco::Candidate &, int &, int,
						 bool);
	const reco::Candidate* get_last_in_decay_chain(const reco::Candidate* );
	void get_stable_daughters(const reco::Candidate&,
							std::vector<const reco::Candidate *>& );
	int tau_classifier(std::vector<const reco::Candidate*>& );
	
	//
	template <class lepton>
	int Print_event_in_file1(FILE *, lepton &, std::vector<pat::Jet> &,
							 CU_ttH_EDA_event_vars &);
	template <class lep1, class lep2>
	// template<class lepton2>
	int Print_event_in_file1_dilepton(FILE *, lep1 &, lep2 &, double,
									  std::vector<pat::Jet> &,
									  CU_ttH_EDA_event_vars &);

	/*
	* Variable section
	*/

	analysis_types analysis_type;
	std::string config_analysis_type;
	
	/// debug flags
	bool verbose_;
	bool dumpHLT_;
	std::string hltTag;
	std::string filterTag;
	
	edm_Tokens token; // common tokens for all events

	/// Triggers, paths: configs filled/updated via run
	HLTConfigProvider hlt_config;
	HLTConfigProvider filter_config;

	/// Triggers, paths.
	// Used for trigger statistics, filled via run (trigger_stats = true)
	bool trigger_stats;
	
	// counters (trigger_stats = true)
	std::map<std::string, unsigned long> n_trigger_fired; // HLT
	std::map<std::string, unsigned long> n_filter_fired;

	std::vector<std::string> trigger_names;
	std::vector<std::string> trigger_names_no_ver;

	std::vector<std::string> filter_names;
	std::vector<std::string> filter_names_no_ver;

	// triggers of interest. Provided by config
	std::vector<std::string> trigger_on_HLT_e;	// single electron trigger
	std::vector<std::string> trigger_on_HLT_mu;   // single muon trigger
	std::vector<std::string> trigger_on_HLT_ee;   // dielectron tigger
	std::vector<std::string> trigger_on_HLT_emu;  // electron+muon trigger
	std::vector<std::string> trigger_on_HLT_mumu; // dimuon trigger

	/// Output file is opened/closed through CMS py config
	edm::Service<TFileService> fs_;

	/// Common sample parameters
	unsigned long event_count; // running event counter

	double total_xs;	  // total cross section
	double int_lumi;	  // integrated luminosity
	double sample_n;	  // total nr of events. Should be long if compatible
	double weight_sample; // int lumi * xs / sample_n
	double weight_gen;

	/// primary vertex
	reco::Vertex pv;
	
	/// Cuts
	float min_ldg_lepton_pT;
	float min_subldg_lepton_pT;
	float min_tau_pT;
	float min_jet_pT;
	float min_bjet_pT;
	float max_jet_eta;
	float max_bjet_eta;
	int min_njets;
	int min_nbtags;

	std::vector<string> cuts;

	std::string jet_corrector;

	/// Selection helper
	MiniAODHelper miniAODhelper;
	bool isdata;
	char MAODHelper_b_tag_strength;

	int MAODHelper_sample_nr; // past insample_, in-development var. for
							  // MAODHelper?
	std::string MAODHelper_era;

	/// Histograms
	TH1D *h_tth_syncex1_ele;
	TH1D *h_tth_syncex1_mu;
	TH1D *h_tth_syncex1_dimu;
	TH1D *h_tth_syncex1_diele;
	TH1D *h_tth_syncex1_elemu;

	TH1D *h_tth_syncex_dileptauh;
	TH1D *h_tth_syncex_eleditauh;
	TH1D *h_tth_syncex_muditauh;

	TH1D *h_hlt;
	TH1D *h_flt;

	// 	TH1D* h_electron_selection;
	// 	TH1D* h_muon_selection;

	// Jet multiplicity
	TH1D *h_njets;
	TH1D *h_nbtags;

	// Tau
	TH1D *h_ntauID;
	// Tau selection efficiency
	TH1D *h_num_genHadTau;
	TH1D *h_genHadTau_pt;
	TH1D *h_genHadTau_eta;
	TH1D *h_genHadTau_phi;
	//TH1D *h_num_FakeTau;
	TH1D *h_num_selectedTau_noniso;
	TH1D *h_selectedTau_noniso_genpt;
	TH1D *h_selectedTau_noniso_geneta;
	TH1D *h_selectedTau_noniso_genphi;
	TH1D *h_num_selectedTau_loose;
	TH1D *h_selectedTau_loose_genpt;
	TH1D *h_selectedTau_loose_geneta;
	TH1D *h_selectedTau_loose_genphi;
	TH1D *h_num_selectedTau_medium;
	TH1D *h_selectedTau_medium_genpt;
	TH1D *h_selectedTau_medium_geneta;
	TH1D *h_selectedTau_medium_genphi;
	TH1D *h_num_selectedTau_tight;
	TH1D *h_selectedTau_tight_genpt;
	TH1D *h_selectedTau_tight_geneta;
	TH1D *h_selectedTau_tight_genphi;


	/// Write-out files
	FILE *events_e_cut1, *events_e_cut2, *events_e_cut3, *events_e_cut4,
		*events_e_cut5, *events_e_cut6, *events_e_cut7;

	FILE *events_mu_cut1, *events_mu_cut2, *events_mu_cut3, *events_mu_cut4,
		*events_mu_cut5, *events_mu_cut6, *events_mu_cut7;

	FILE *events_dimu_cut1, *events_dimu_cut2, *events_dimu_cut3,
		*events_dimu_cut4, *events_dimu_cut5, *events_dimu_cut6,
		*events_dimu_cut7;

	FILE *events_diele_cut1, *events_diele_cut2, *events_diele_cut3,
		*events_diele_cut4, *events_diele_cut5, *events_diele_cut6,
		*events_diele_cut7;

	FILE *events_elemu_cut1, *events_elemu_cut2, *events_elemu_cut3,
		*events_elemu_cut4, *events_elemu_cut5;

	/// Technical variables

	// 	int return_code;

	/// Legacy
	// 	int NjetMin;
	// 	int NjetMax;
	// 	int NjetBins;
	// 	int NtagMin;
	// 	int NtagMax;
	// 	int NtagBins;
	// 	int NpuMin;
	// 	int NpuMax;
	// 	int NpuBins;
	//
	// 	int numEvents_;
	// 	int numHbb;
	// 	int numHdecayToOneParticle;
	// 	int numHdecayToOneGluon;
	// 	int numHdecayToOnePhoton;
	// 	int numHdecayToOneZBoson;
	// 	int numHdecayToOneWBoson;
	//
	// 	double mySample_xSec_;
	// 	double mySample_nGen_;
	// 	double intLumi_;
	//
	// 	double ptmax;
	// 	int NptBins;
	//
	// 	double jetptmax;
	// 	int NjetptBins;

	/// tree & branches
	TTree *eventTree;
	
	// All sorted by pt
	int n_electrons;
	int n_muons;
	int n_loose_taus;
	int n_medium_taus;
	int n_tight_taus;
	int n_jets;
	int n_btags;

	float pv_x;
	float pv_y;
	float pv_z;
	
	//std::vector<TLorentzVector> electrons; // sorted by pt
	std::vector<float> e_pt;
	std::vector<float> e_eta;
	std::vector<float> e_phi;
	std::vector<float> e_mass;
	std::vector<int> e_charges;
	std::vector<float> e_vz;
	std::vector<float> e_vx;
	std::vector<float> e_vy;
	std::vector<float> e_vtx_dz;
	std::vector<float> e_vtx_dxy;
	std::vector<bool> e_isGsfCtfScPixChargeConsistent;
	
	//std::vector<TLorentzVector> muons; // sorted by pt
	std::vector<float> mu_pt;
	std::vector<float> mu_eta;
	std::vector<float> mu_phi;
	std::vector<float> mu_mass;
	std::vector<int> mu_charges;
	std::vector<float> mu_vz;
	std::vector<float> mu_vx;
	std::vector<float> mu_vy;
	std::vector<float> mu_vtx_dz;
	std::vector<float> mu_vtx_dxy;
	std::vector<float> mu_relTrkPtError;
	
	//std::vector<TLorentzVector> loose_taus; // sorted by pt
	std::vector<float> loose_tau_pt;
	std::vector<float> loose_tau_eta;
	std::vector<float> loose_tau_phi;
	std::vector<float> loose_tau_mass;
	//std::vector<TLorentzVector> medium_taus; // sorted by pt
	std::vector<float> medium_tau_pt;
	std::vector<float> medium_tau_eta;
	std::vector<float> medium_tau_phi;
	std::vector<float> medium_tau_mass;
	//std::vector<TLorentzVector> tight_taus; // sorted by pt
	std::vector<float> tight_tau_pt;
	std::vector<float> tight_tau_eta;
	std::vector<float> tight_tau_phi;
	std::vector<float> tight_tau_mass;
	std::vector<int> Ltau_charges;
	std::vector<int> Mtau_charges;
	std::vector<int> Ttau_charges;
	std::vector<float> Ltau_vz;
	std::vector<float> Mtau_vz;
	std::vector<float> Ttau_vz;
	std::vector<float> Ltau_vx;
	std::vector<float> Mtau_vx;
	std::vector<float> Ttau_vx;
	std::vector<float> Ltau_vy;
	std::vector<float> Mtau_vy;
	std::vector<float> Ttau_vy;
	std::vector<float> Ltau_vtx_dz;
	std::vector<float> Mtau_vtx_dz;
	std::vector<float> Ttau_vtx_dz;
	std::vector<float> Ltau_vtx_dxy;
	std::vector<float> Mtau_vtx_dxy;
	std::vector<float> Ttau_vtx_dxy;
	std::vector<int> Ltau_decaymode;
	std::vector<int> Mtau_decaymode;
	std::vector<int> Ttau_decaymode;
	std::vector<size_t> Ltau_ntrk;
	std::vector<size_t> Mtau_ntrk;
	std::vector<size_t> Ttau_ntrk;
	// dxy_PCA
	// dxy_Sig
	// flightLength
	// secondary vertex
	// leading track
	
	
	//std::vector<TLorentzVector> jets; // jets_selected_sorted
	std::vector<float> jet_pt;
	std::vector<float> jet_eta;
	std::vector<float> jet_phi;
	std::vector<float> jet_mass;
	std::vector<float> jet_charges;
	std::vector<float> jet_vz;
	std::vector<float> jet_vx;
	std::vector<float> jet_vy;
	//std::vector<TLorentzVector> bjets; // jets_selected_tag_sorted
	std::vector<float> bjet_pt;
	std::vector<float> bjet_eta;
	std::vector<float> bjet_phi;
	std::vector<float> bjet_mass;
	std::vector<float> bjet_charges;
	std::vector<float> bjet_vz;
	std::vector<float> bjet_vx;
	std::vector<float> bjet_vy;

	// MET
	float MET_x;
	float MET_y;
	
	// MC truth
	std::vector<int> gen_x_pdgId;
	std::vector<int> gen_x_status;
	//std::vector<TLorentzVector> gen_x;
	std::vector<float> gen_x_pt;
	std::vector<float> gen_x_eta;
	std::vector<float> gen_x_phi;
	std::vector<float> gen_x_mass;
	std::vector<float> gen_x_vx;
	std::vector<float> gen_x_vy;
	std::vector<float> gen_x_vz;
	
	std::vector<int> gen_top_pdgId;
	std::vector<int> gen_top_status;
	//std::vector<TLorentzVector> gen_top;
	std::vector<float> gen_top_pt;
	std::vector<float> gen_top_eta;
	std::vector<float> gen_top_phi;
	std::vector<float> gen_top_mass;
	std::vector<float> gen_top_vx;
	std::vector<float> gen_top_vy;
	std::vector<float> gen_top_vz;
	
	std::vector<int> gen_xDaug_pdgId;
	std::vector<int> gen_xDaug_status;
	//std::vector<TLorentzVector> gen_xDaug;
	std::vector<float> gen_xDaug_pt;
	std::vector<float> gen_xDaug_eta;
	std::vector<float> gen_xDaug_phi;
	std::vector<float> gen_xDaug_mass;
	std::vector<float> gen_xDaug_vx;
	std::vector<float> gen_xDaug_vy;
	std::vector<float> gen_xDaug_vz;

	std::vector<int> gen_tau_class;

	std::vector<int> gen_topDaug_pdgId;
	std::vector<int> gen_topDaug_status;
	//std::vector<TLorentzVector> gen_topDaug;
	std::vector<float> gen_topDaug_pt;
	std::vector<float> gen_topDaug_eta;
	std::vector<float> gen_topDaug_phi;
	std::vector<float> gen_topDaug_mass;
	std::vector<float> gen_topDaug_vx;
	std::vector<float> gen_topDaug_vy;
	std::vector<float> gen_topDaug_vz;
	
	std::vector<int> gen_wDaug_pdgId;
	std::vector<int> gen_wDaug_status;
	//std::vector<TLorentzVector> gen_wDaug;
	std::vector<float> gen_wDaug_pt;
	std::vector<float> gen_wDaug_eta;
	std::vector<float> gen_wDaug_phi;
	std::vector<float> gen_wDaug_mass;
	std::vector<float> gen_wDaug_vx;
	std::vector<float> gen_wDaug_vy;
	std::vector<float> gen_wDaug_vz;
};



template <typename T1, typename T2>
std::vector<T1>
CU_ttH_EDA::removeOverlap(const std::vector<T1> &v1, const std::vector<T2> &v2, double dR)
{
	std::vector<T1> res;
	for (const auto& o1: v1) {
		bool keep = true;
		for (const auto& o2: v2)
			if (miniAODhelper.DeltaR(&o1, &o2) < dR)
				keep = false;
		if (keep)
			res.push_back(o1);
	}
	return res;
}

#endif
