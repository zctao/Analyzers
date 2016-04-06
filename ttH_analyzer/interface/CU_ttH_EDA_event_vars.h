#ifndef CU_ttH_EDA_event_vars_h
#define CU_ttH_EDA_event_vars_h

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

#include <vector>

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
	int n_taus;
	int n_jets;
	int n_btags;
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
	std::vector<pat::Tau> tau_selected;
	std::vector<pat::Tau> tau_selected_sorted;
	
	std::vector<pat::Jet> jets_raw;
	std::vector<pat::Jet> jets_no_mu;
	std::vector<pat::Jet> jets_no_mu_e;
	std::vector<pat::Jet> jets_corrected;
	std::vector<pat::Jet> jets_selected;
	std::vector<pat::Jet> jets_selected_sorted;
	std::vector<pat::Jet> jets_selected_tag;
	std::vector<pat::Jet> jets_selected_tag_sorted;

	/// Other quantities
	pat::MET pfMET;
	pat::MET MET_corrected;
	double dimuon_mass;
	double dielectron_mass;
	double dilepton_mass;

	double MHT;
	double metLD;
	
	// Gen Particles
	std::vector<reco::GenParticle> genHiggs;
	std::vector<reco::GenParticle> genTops;
	reco::CandidateCollection genHiggs_daughters; // or edm::OwnVector<reco::Candidate>
	reco::CandidateCollection genTop_daughters;
	reco::CandidateCollection genW_daughters;
};

#endif
