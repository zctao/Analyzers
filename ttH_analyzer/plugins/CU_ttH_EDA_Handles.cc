#ifndef CU_ttH_EDA_Handles_CC
#define CU_ttH_EDA_Handles_CC

#include "CU_ttH_EDA_Handles.h"

/// Set up handles with getByToken from edm::Event
void Set_up_handles(const Event &iEvent, edm_Handles &handle, edm_Tokens &token,
					bool isdata)
{		
	iEvent.getByToken(token.triggerResults, handle.triggerResults);
	iEvent.getByToken(token.filterResults, handle.filterResults);

	iEvent.getByToken(token.vertices, handle.vertices);
	iEvent.getByToken(token.sec_vertices, handle.sec_vertices);
	iEvent.getByToken(token.PU_info, handle.PU_info);
	iEvent.getByToken(token.srcRho, handle.srcRho);

	iEvent.getByToken(token.electrons, handle.electrons);
	iEvent.getByToken(token.muons, handle.muons);
	iEvent.getByToken(token.taus, handle.taus);

	iEvent.getByToken(token.jets, handle.jets);
	iEvent.getByToken(token.METs, handle.METs);

	iEvent.getByToken(token.PF_candidates, handle.PF_candidates);

	iEvent.getByToken(token.BS, handle.BS);

	//iEvent.getByToken(token.top_jets, handle.top_jets);
	//iEvent.getByToken(token.subfilter_jets, handle.subfilter_jets);

	if (!isdata) {
		iEvent.getByToken(token.event_gen_info, handle.event_gen_info);
		iEvent.getByToken(token.MC_particles, handle.MC_particles);
		iEvent.getByToken(token.MC_packed, handle.MC_packed);
		iEvent.getByToken(token.genJets, handle.genJets);
	}
}

#endif
