#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

using namespace std;

class MCHiggsDecayModeFilter : public edm::EDFilter {
public:
	explicit MCHiggsDecayModeFilter(const edm::ParameterSet &);
	~MCHiggsDecayModeFilter();

private:
	virtual void beginJob();
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual void endJob();

	// ----------member data ---------------------------
	edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken_;
	string DecayMode;
	bool Enable;

	//map<string, int> decayid = {{"hww", 24}, {"hzz", 23}, {"htt", 15}, {"hbb", 5}};
};


MCHiggsDecayModeFilter::MCHiggsDecayModeFilter(const edm::ParameterSet & iConfig) :
	DecayMode (iConfig.getParameter<string>("DecayMode")),
	Enable (iConfig.getParameter<bool>("Enable"))
{
	prunedGenToken_ = consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));
}

MCHiggsDecayModeFilter::~MCHiggsDecayModeFilter()
{
}

bool
MCHiggsDecayModeFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	if (not Enable)	return true;
	
	// Get GenParticles
	edm::Handle<reco::GenParticleCollection> GenParticles;
	iEvent.getByToken(prunedGenToken_, GenParticles);

	for (size_t i = 0; i < GenParticles->size(); ++i) {
		const reco::GenParticle *p = &(*GenParticles)[i];

		// Get Higgs
		if (p->pdgId() != 25) continue;

		int ndaugs = p->numberOfDaughters();
		if (ndaugs != 2) continue;
		const reco::Candidate *d1 = p->daughter(0);
		const reco::Candidate *d2 = p->daughter(1);
		
		if ( abs(d1->pdgId()) != abs(d2->pdgId()) ) continue;

		int daug_id = abs(d1->pdgId());

		int required_id = 0;	
		if (DecayMode == "ttH_htt")
			required_id = 15;
		else if (DecayMode == "ttH_hww")
			required_id = 24;
		else if (DecayMode == "ttH_hzz")
			required_id = 23;

		return daug_id == required_id;
		
	}

	return false;
}

// ------------ method called once each job just before starting event loop  ------------
void
MCHiggsDecayModeFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MCHiggsDecayModeFilter::endJob()
{
}

// define this as a plug-in
DEFINE_FWK_MODULE(MCHiggsDecayModeFilter);
