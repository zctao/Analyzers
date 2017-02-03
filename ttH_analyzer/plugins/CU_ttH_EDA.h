#ifndef CU_ttH_EDA_h
#define CU_ttH_EDA_h

/// Core libraries
#include <memory>
#include <cstdio>	// printf, fprintf
#include <stdexcept> // standard exceptions
#include <map>

/// CMSSW user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
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
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

/// BTag Calibration
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

//#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"

/// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TAxis.h"
#include "TFile.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"

/// Higgs and top tagger
//#include "MiniAOD/BoostedObjects/interface/HTTTopJet.h"
//#include "MiniAOD/BoostedObjects/interface/SubFilterJet.h"
//#include "BoostedTTH/BoostedAnalyzer/interface/BoostedUtils.hpp"
#include "MiniAOD/MiniAODHelper/interface/TopTagger.h"
#include "MiniAOD/MiniAODHelper/interface/HiggsTagger.h"

/// structs for holding multiple edm::Handle and EDGetTokenT
#include "CU_ttH_EDA_Handles.h"

#include "Analyzers/ttH_analyzer/interface/miniLepton.h"
#include "Analyzers/ttH_analyzer/interface/CU_ttH_EDA_Ntuple.h"
#include "Analyzers/ttH_analyzer/interface/kinMVA_ttbar_2lss.h"
#include "Analyzers/ttH_analyzer/interface/kinMVA_ttV_2lss.h"

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
	Analyze_1l2tau,
	Analyze_2lss1tau,
	Analyze_3l
};

// enum for event selection region
enum Selection_types {
	Signal_2lss1tau,
	Signal_1l2tau,
	Signal_3l,
	Control_2los1tau,
	Control_1lfakeable,
	Control_WZ
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

	// virtual void beginLuminosityBlock(const edm::LuminosityBlock &,
	//								  const edm::EventSetup &) override;
	// virtual void endLuminosityBlock(const edm::LuminosityBlock &,
	//								const edm::EventSetup &) override;

	/// One-time-run functions
	void Load_configuration(string); // at CU_ttH_EDA(), runs _MAODH()
	void Load_configuration_set_type(const string &); // sets analysis_type
	void Load_configuration_MAODH(bool); // runs miniAODhelper.SetUp
	void Set_up_histograms();			 // at CU_ttH_EDA()
	void Set_up_tokens(const edm::ParameterSet &);
	void Set_up_Tree();
	void Set_up_selection_region(const string &);
	void Set_up_BTagCalibration_Readers();
	void Set_up_CSV_rootFile();
	void fillCSVhistos(TFile*, TFile*);
	void Set_up_LeptonSF_Lut();
	void Set_up_FakeRate_Lut();
	void Set_up_TauSF_Lut();
	void Set_up_ChargeMisID_Lut();
	void Set_up_PUWeight_hist();

	int Set_up_Run_histograms_triggers(); // at beginRun(), after
										  // Set_up_name_vectors()

	void Set_up_trigger_name_vectors(); // at beginRun()
	void Set_up_HLT_path_name();
	void Add_HLT_pathName_version(std::vector<std::string>&,
								  std::vector<std::string>&);

	int End_Run_hist_fill_triggers(); // fill histograms at endRun()

	/// Per-event functions
	void Update_common_vars(const edm::Event &, CU_ttH_EDA_event_vars &);

	/// Object checks. Returns 1 in case of an error
	int Check_beam_spot(edm::Handle<reco::BeamSpot>);
	int Check_triggers(edm::Handle<edm::TriggerResults>,
					   CU_ttH_EDA_event_vars &); // adjusts event variables
	int Check_filters(edm::Handle<edm::TriggerResults>);
	int Check_vertices_set_MAODhelper(edm::Handle<reco::VertexCollection>);

	// trigger iterator, part of Check_triggers()
	bool Check_triggers_iterator(const vector<string> &,
								 edm::Handle<edm::TriggerResults>);

	std::vector<pat::Jet> GetCorrectedJets(const std::vector<pat::Jet>&, JetCorrectionUncertainty*, const std::string);
	std::vector<pat::Tau> GetCorrectedTaus(const std::vector<pat::Tau>&, float, const std::string);
	
	void printDecayChain(const reco::Candidate &p, int &index, int mother_index,
						 bool details);
	int HiggsDaughterPdgId(const std::vector<reco::GenParticle>&);
	bool HiggsDecayFilter(const std::vector<reco::GenParticle>&, const TString&);

	template <typename T1, typename T2>
		std::vector<T1>
		removeOverlapdR(const std::vector<T1>& v1, const std::vector<T2>& v2, double dR = 0.02);

	float getMHT(CU_ttH_EDA_event_vars &);

	// MC truth matching
	template <typename T>
		int MatchGenParticle_Type(const T&, const std::vector<reco::GenParticle>&);
	const reco::GenParticle* getMatchedGenParticle(const pat::Electron&, const std::vector<reco::GenParticle>&);
	const reco::GenParticle* getMatchedGenParticle(const pat::Muon&, const std::vector<reco::GenParticle>&);
	const reco::GenParticle* getMatchedGenParticle(const pat::Tau&, const std::vector<reco::GenParticle>&);
	
	// event selection
	bool pass_event_sel_2l(CU_ttH_EDA_event_vars &, Selection_types, int&, int&);
	bool pass_event_sel_3l(CU_ttH_EDA_event_vars &, Selection_types);

	// MVA
	int partition2DBDT(double, double);

	// csv reweighting
	//double getEvtCSVWeight(std::vector<pat::Jet> &, int); // use root file
	double getEvtCSVWeight(std::vector<pat::Jet> &, const std::string &); // use csv file
	double getJetCSVWeight(pat::Jet &, std::string /*pass by copy*/);

	float read2DHist(TH2*, float, float);
	float readTGraph(TGraphAsymmErrors*, float);
	float evalTGraph(TGraphAsymmErrors*, float);
	float readTF(TF1*, float);
	
	// Electron Charge misId
	float getEleChargeMisIDProb(const miniLepton&, int);
	
	// Lepton fake rate
	float getFakeRate(const miniLepton&);
	float getFakeRate(const pat::Tau&);

	// HLT scale factor
	float getLepHLTSF(int);
	
	// lepton scale factor
	float getLeptonSF(const miniLepton&);
	float getLeptonSF_loose(const miniLepton&);
	float getLeptonSF_tight_vs_loose(const miniLepton&);

	// PU weight
	float getPUWeight(int);
	
	/*
	* Variable section
	*/

	// Analysis type
	analysis_types analysis_type;
	std::string config_analysis_type;
	
	// flag for sync ntuple
	bool turn_off_event_sel;

	bool doSystematics;

	/// Common sample parameters
	bool doLumiScale;
	unsigned long event_count; // running event counter
	TString sampleName; 
	double sample_xs;	  // cross section	
	double int_lumi;	  // integrated luminosity
	double sample_n;	  // total nr of events. Should be long if compatible
	double genWeightSum;
	double genWeightxPUSum;
	//double weight_sample; // int lumi * xs / sample_n
	
	/// debug flags
	bool verbose_;
	bool dumpHLT_;
	
	edm_Tokens token; // common tokens for all events

	/// Triggers, paths: configs filled/updated via run
	bool hltcut_off;
	HLTConfigProvider hlt_config;
	HLTConfigProvider filter_config;

	/// Triggers, paths.
	// Used for trigger statistics, filled via run (trigger_stats = true)
	std::string hltTag;
	std::string filterTag;

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
	
	/// Selection helper
	MiniAODHelper miniAODhelper;

	// tauES
	std::string TESType;  // "tauESUp" or "tauESDown" or  "NA"
	
	//JEC
	std::string JECType;
	// sysType::sysType defined in MiniAODHelper.h
	std::map<std::string,sysType::sysType> JECTypes  
		= {{"NA", sysType::NA},
		   {"JERUp", sysType::JERup},{"JERDown", sysType::JERdown},
		   {"JESUp", sysType::JESup},{"JESDown", sysType::JESdown}};
	bool doJERsmear;

	// CSV WP
	double csv_loose_wp;
	double csv_medium_wp;
	double csv_tight_wp;
	
	bool isdata;
	char MAODHelper_b_tag_strength;
	int MAODHelper_sample_nr; // past insample_, in-development var. for
							  // MAODHelper?
	std::string MAODHelper_era;
	
	// event selection region
	std::string selection_region;
    Selection_types selection_type;

	// debug flag
	bool debug;
	
	/// Histograms
	TH1D *h_hlt;
	TH1D *h_flt;

	TH1I *h_nProcessed;
	TH1D *h_SumGenWeight;
	TH1D *h_SumGenWeightxPU;
	TH1D *h_GenWeightProcessed;  // cross check
	TH1D *h_GenWeightxPUProcessed;  // cross check
	
	//TH2D *h_MVA_ttV_vs_ttbar[3][2];
	//TH1D *h_MVA_shape[3][2];
	//TH2D *h_MVA_ttV_vs_ttbar_sys[3][2][16];
	//TH1D *h_MVA_shape_csv_sys[3][2][16];
	//TH1D *h_MVA_shape_thu_sys[3][2][4];
	//bool setup_sysHist = false;

	// for WZ control region
	//TH1D *h_mTWl;  // mT of the lepton from W
	//TH1D *h_met;
	//TH1D *h_nJets;
	//TH1D *h_lep_charge;
	//TH1D *h_mZ;

	// Fake lepton rate lookup histograms
	TFile *file_fr_lep;
	TH2F *h_fakerate_el;
	TH2F *h_fakerate_mu;
	// jet to tau fake rate
	TFile* file_fr_tau;
	TGraphAsymmErrors *g_fakerate_tau_mvaM_etaL_mc;
	TGraphAsymmErrors *g_fakerate_tau_mvaM_etaH_mc;
	//TGraphAsymmErrors *g_fakerate_tau_mvaT_etaL_mc;
	//TGraphAsymmErrors *g_fakerate_tau_mvaT_etaH_mc;
	TF1 *f_fakerate_tau_mvaM_etaL_ratio;
	TF1 *f_fakerate_tau_mvaM_etaH_ratio;
	//TF1 *f_fakerate_tau_mvaT_etaL_ratio;
	//TF1 *f_fakerate_tau_mvaT_etaH_ratio;

	// Electron Charge MisID lookup histogram
	TFile *file_eleMisCharge;
	TH2F *h_chargeMisId;
	
	// Lepton ID scale factor lookup tables
	// Input files
	TFile* file_recoToLoose_leptonSF_mu1_b;
	TFile* file_recoToLoose_leptonSF_mu1_e;
	TFile* file_recoToLoose_leptonSF_mu2;
	TFile* file_recoToLoose_leptonSF_mu3;
	TFile* file_recoToLoose_leptonSF_el;
	TFile* file_recoToLoose_leptonSF_gsf;
	TFile* file_looseToTight_leptonSF_mu_2lss;
	TFile* file_looseToTight_leptonSF_el_2lss;
	TFile* file_looseToTight_leptonSF_mu_3l;
	TFile* file_looseToTight_leptonSF_el_3l;
	
	TGraphAsymmErrors *h_recoToLoose_leptonSF_mu1_b;
	TGraphAsymmErrors *h_recoToLoose_leptonSF_mu1_e;
	TH2F *h_recoToLoose_leptonSF_mu2;
	TGraphAsymmErrors *h_recoToLoose_leptonSF_mu3;
	TH2F *h_recoToLoose_leptonSF_el1;
	TH2F *h_recoToLoose_leptonSF_el2;
	TH2F *h_recoToLoose_leptonSF_el3;
	TH2F *h_recoToLoose_leptonSF_gsf;
	TH2F *h_looseToTight_leptonSF_mu_2lss;
	TH2F *h_looseToTight_leptonSF_el_2lss;
	TH2F *h_looseToTight_leptonSF_mu_3l;
	TH2F *h_looseToTight_leptonSF_el_3l;

	// PU weight
	TFile* file_puweight;
	TH1F *h_puweight;
	
	// tree and ntuple
	TTree *eventTree;
	CU_ttH_EDA_Ntuple evtNtuple;

	// MVA
	kinMVA_ttbar_2lss MVA_ttbar_vars;
	kinMVA_ttV_2lss MVA_ttV_vars;
	
	std::string sysList[16] =
		{"LFUp","LFDown","HFUp","HFDown",
		 "HFStats1Up","HFStats1Down","HFStats2Up","HFStats2Down",
		 "LFStats1Up","LFStats1Down","LFStats2Up","LFStats2Down",
		 "cErr1Up","cErr1Down","cErr2Up","cErr2Down"};

	//std::map<std::string, BTagCalibrationReader*> BTagCaliReaders;
	BTagCalibrationReader* BTagCaliReader;
	/*
	// or read SF from root file
	TFile* f_CSVwgt_HF;
	TFile* f_CSVwgt_LF;
	TH1D* h_csv_wgt_hf[9][5];
	TH1D* c_csv_wgt_hf[9][5];
	TH1D* h_csv_wgt_lf[9][4][3];
	
	// CSV reweight iSys map
	std::map<std::string, int> csv_iSys
		= {{"JESUp",7}, {"JESDown",8}, {"LFUp",9}, {"LFDown",10}, {"HFUp",11},
		   {"HFDown",12}, {"HFStats1Up",13}, {"HFStats1Down",14},
		   {"HFStats2Up",15}, {"HFStats2Down",16}, {"LFStats1Up",17},
		   {"LFStats1Down",18}, {"LFStats2Up",19}, {"LFStats2Down",20},
		   {"cErr1Up",21},{"cErr1Down",22},{"cErr2Up",23},{"cErr2Down",24},
		   {"NA",0}};
	*/
	
};

template <typename T1, typename T2>
std::vector<T1>
CU_ttH_EDA::removeOverlapdR(const std::vector<T1> &v1, const std::vector<T2> &v2, double dR)
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
