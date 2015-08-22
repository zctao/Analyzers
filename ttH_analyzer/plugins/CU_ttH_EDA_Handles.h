#ifndef CU_ttH_EDA_Handles_h
#define CU_ttH_EDA_Handles_h

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

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
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1Trigger/interface/L1HFRingsFwd.h"

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

/// Higgs and top tagger
#include "MiniAOD/BoostedObjects/interface/HTTTopJet.h"
#include "MiniAOD/BoostedObjects/interface/SubFilterJet.h"
#include "BoostedTTH/BoostedAnalyzer/interface/BoostedUtils.hpp"

/*
 *
 * edm::Handle and edm::EDGetTokenT container structs.
 *
 */

using namespace edm;

struct edm_Handles
{
	Handle<GenEventInfoProduct> event_gen_info;
	Handle<edm::TriggerResults> triggerResults;
	Handle<edm::TriggerResults> filterResults;

	Handle<reco::VertexCollection> vertices;
	Handle<reco::VertexCompositePtrCandidateCollection> sec_vertices;
	Handle<std::vector<PileupSummaryInfo>> PU_info;

	Handle<pat::ElectronCollection> electrons;
	Handle<pat::MuonCollection> muons;

	Handle<pat::JetCollection> jets;
	Handle<pat::METCollection> METs;

	Handle<reco::GenParticleCollection> MC_particles;
	Handle<pat::PackedCandidateCollection> PF_candidates;

	Handle<reco::BeamSpot> BS;

	Handle<boosted::HTTTopJetCollection> top_jets;
	Handle<boosted::SubFilterJetCollection> subfilter_jets;
};

struct edm_Tokens
{
	EDGetTokenT<GenEventInfoProduct> event_gen_info;
	EDGetTokenT<edm::TriggerResults> triggerResults;
	EDGetTokenT<edm::TriggerResults> filterResults;

	EDGetTokenT<reco::VertexCollection> vertices;
	EDGetTokenT<reco::VertexCompositePtrCandidateCollection> sec_vertices;
	EDGetTokenT<std::vector<PileupSummaryInfo>> PU_info;

	EDGetTokenT<pat::ElectronCollection> electrons;
	EDGetTokenT<pat::MuonCollection> muons;

	EDGetTokenT<pat::JetCollection> jets;
	EDGetTokenT<pat::METCollection> METs;

	EDGetTokenT<reco::GenParticleCollection> MC_particles;
	EDGetTokenT<pat::PackedCandidateCollection> PF_candidates;

	EDGetTokenT<reco::BeamSpot> BS;

	EDGetTokenT<boosted::HTTTopJetCollection> top_jets;
	EDGetTokenT<boosted::SubFilterJetCollection> subfilter_jets;
};

/// Set up handles with getByToken from edm::Event
void Set_up_handles(const Event &, edm_Handles &, edm_Tokens &);

#endif
