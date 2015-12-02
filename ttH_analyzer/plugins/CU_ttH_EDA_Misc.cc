#ifndef CU_ttH_EDA_Misc_cc
#define CU_ttH_EDA_Misc_cc

/// Includes
#include "CU_ttH_EDA.h"

int CU_ttH_EDA::End_Run_hist_fill_triggers()
{
	TAxis *axis = h_hlt->GetXaxis();
	if (!axis)
		return 1;

	for (std::map<std::string, unsigned long>::const_iterator iter =
			 n_trigger_fired.begin();
		 iter != n_trigger_fired.end(); ++iter) {

		int bin_num = axis->FindBin((iter->first).c_str());
		h_hlt->Fill(bin_num - 1, (iter->second));
	}

	axis = h_flt->GetXaxis();
	if (!axis)
		return 1;

	for (std::map<std::string, unsigned long>::const_iterator iter =
			 n_filter_fired.begin();
		 iter != n_filter_fired.end(); ++iter) {

		int bin_num = axis->FindBin((iter->first).c_str());
		h_flt->Fill(bin_num - 1, (iter->second));
	}

	return 0;
}

void CU_ttH_EDA::Update_common_vars(const edm::Event &iEvent,
									CU_ttH_EDA_event_vars &local)
{
	local.run_nr = iEvent.id().run();
	local.event_nr = iEvent.id().event();
	local.lumisection_nr = iEvent.id().luminosityBlock();
}

int CU_ttH_EDA::Check_beam_spot(edm::Handle<reco::BeamSpot> BS)
{
	if (!BS.isValid())
		return 1;

	// 	math::XYZPoint beamSpotPosition;
	// 	beamSpotPosition.SetCoordinates(0, 0, 0);
	// //		double BSx = 0, BSy = 0, BSz = 0;
	//
	// 	double BSx = BS->x0();
	// 	double BSy = BS->y0();
	// 	double BSz = BS->z0();
	// 	beamSpotPosition = BS->position();
	//
	//		if (verbose_)
	//			printf("\t BeamSpot: x = %.2f,\t y = %.2f,\t z = %.2f \n",
	//				BSx, BSy, BSz );
	if (verbose_)
		printf("\t BeamSpot: x = %.2f,\t y = %.2f,\t z = %.2f \n", BS->x0(),
			   BS->y0(), BS->z0());

	return 0;
}

int CU_ttH_EDA::Check_triggers(edm::Handle<edm::TriggerResults> triggerResults,
							   CU_ttH_EDA_event_vars &local)
{
	if (!triggerResults.isValid()) {
		std::cerr << "Trigger results not valid for tag " << hltTag
				  << std::endl;
		return 1;
	}

	/// scan: trigger_on_HLT_<type>
	local.pass_single_e =
		Check_triggers_iterator(trigger_on_HLT_e, triggerResults);
	local.pass_single_mu =
		Check_triggers_iterator(trigger_on_HLT_mu, triggerResults);
	local.pass_double_e =
		Check_triggers_iterator(trigger_on_HLT_ee, triggerResults);
	local.pass_elemu =
		Check_triggers_iterator(trigger_on_HLT_emu, triggerResults);
	local.pass_double_mu =
		Check_triggers_iterator(trigger_on_HLT_mumu, triggerResults);

	if (trigger_stats) {
		for (std::vector<std::string>::const_iterator trigger =
				 trigger_names.begin();
			 trigger != trigger_names.end(); ++trigger) {
			std::string pathName = *trigger;
			unsigned int hltIndex = hlt_config.triggerIndex(pathName);

			if (hltIndex >= triggerResults->size())
				continue;

			bool trigger_accept = triggerResults->accept(hltIndex);
			int prescale = -1; // hlt_config.prescaleValue(iEvent, iSetup,
							   // pathName);

			if (verbose_ && dumpHLT_)
				std::cout << " =====>  HLT: path name = " << pathName
						  << ",\t prescale = " << prescale
						  << ",\t pass = " << trigger_accept << std::endl;

			std::string pathNameNoVer = hlt_config.removeVersion(pathName);

			if (trigger_accept)
				++n_trigger_fired[pathNameNoVer];
		}
	}

	return 0;
}

bool CU_ttH_EDA::Check_triggers_iterator(
	const vector<string> &triggers,
	edm::Handle<edm::TriggerResults> triggerResults)
{
	for (std::vector<std::string>::const_iterator trigger = triggers.begin();
		 trigger != triggers.end(); ++trigger) {
		unsigned int hltIndex = hlt_config.triggerIndex(*trigger);

		if (hltIndex >= triggerResults->size())
			continue;

		if (triggerResults->accept(hltIndex))
			return true;
	}

	return false;
}

int CU_ttH_EDA::Check_filters(edm::Handle<edm::TriggerResults> filterResults)
{
	if (!filterResults.isValid()) {
		std::cerr << "Trigger results not valid for tag " << filterTag
				  << std::endl;
		return 1;
	}

	if (trigger_stats) {
		for (std::vector<std::string>::const_iterator trigger =
				 filter_names.begin();
			 trigger != filter_names.end(); ++trigger) {
			std::string pathName = *trigger;
			unsigned int hltIndex = filter_config.triggerIndex(pathName);

			if (hltIndex >= filterResults->size())
				continue;

			bool filter_accept = filterResults->accept(hltIndex);
			int prescale = -1; // filter_config.prescaleValue(iEvent, iSetup,
							   // pathName);

			if (verbose_ && dumpHLT_)
				std::cout << " =====>  Filter: path name = " << pathName
						  << ",\t prescale = " << prescale
						  << ",\t pass = " << filter_accept << std::endl;

			std::string pathNameNoVer = filter_config.removeVersion(pathName);

			if (filter_accept)
				++n_filter_fired[pathNameNoVer];
		}
	}

	return 0;
}

int CU_ttH_EDA::Check_vertices_set_MAODhelper(
			    edm::Handle<reco::VertexCollection> vertices,
			    reco::Vertex& vertex)
{
	/// Primary vertex handling
	if (!vertices.isValid())
		return 1;

	//reco::Vertex vertex;
	int n_PVs = 0;

	for (reco::VertexCollection::const_iterator vtx = vertices->begin();
		 vtx != vertices->end(); ++vtx) {

		if (vtx->isFake() || vtx->ndof() < 4.0 || abs(vtx->z()) > 24.0 ||
			abs(vtx->position().Rho()) > 2.0)
			continue;

		if (n_PVs == 0)
			vertex = *vtx;

		++n_PVs;
	}

	if (verbose_)
		printf("\t Event PV: x = %.3f,\t y = %.3f,\t z = %.3f \n", vertex.x(),
			   vertex.y(), vertex.z());

	if (n_PVs > 0)
		miniAODhelper.SetVertex(
			vertex); // FIXME?: overload miniAODhelper::SetVertex(reco::Vertex&)
	
	return 0;
}

/// Taggers
int CU_ttH_EDA::Higgs_tagger(
	Handle<boosted::SubFilterJetCollection> subfilter_jets,
	CU_ttH_EDA_event_vars &local)
{
	local.n_Htags = 0;

	if (!subfilter_jets.isValid())
		return 1;

	boosted::SubFilterJetCollection subfilterjets =
		BoostedUtils::GetSortedByPt(*subfilter_jets);

	for (boosted::SubFilterJetCollection::iterator higgsJet =
			 subfilterjets.begin();
		 higgsJet != subfilterjets.end(); ++higgsJet) {
		// pt and eta requirements for top jet
		if (higgsJet->fatjet.pt() <= 250. || abs(higgsJet->fatjet.eta()) >= 1.8)
			continue;

		int numBtagFiltJets = 0;
		std::vector<pat::Jet> filterjets = higgsJet->filterjets;
		int numFiltJets = filterjets.size();
		for (int ijet = 0; ijet < numFiltJets; ++ijet) {
			if (verbose_) {
				printf("\t\t filt jet %2d:\t pT = %.1f,\t eta = %.2f,\t phi = "
					   "%.2f,\t CSVv2 = %+5.3f,\t CSVv1 = %+5.3f \n",
					   ijet, filterjets[ijet].pt(), filterjets[ijet].eta(),
					   filterjets[ijet].phi(),
					   filterjets[ijet].bDiscriminator(
						   "combinedInclusiveSecondaryVertexV2BJetTags"),
					   filterjets[ijet].bDiscriminator(
						   "combinedSecondaryVertexBJetTags"));
			}

			if (filterjets[ijet].pt() <= 20. ||
				abs(filterjets[ijet].eta()) >= 2.5)
				continue;

			// b-tag medium WP
			if (filterjets[ijet].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags") < 0.814)
				continue;

			++numBtagFiltJets;
		}

		if (verbose_) {
			printf("\t Higgs jet %2d:\t pT = %.1f,\t eta = %.2f,\t phi = "
				   "%.2f,\t numFiltJets = %2d,\t numBtagFiltJets = %2d\n",
				   int(higgsJet - subfilterjets.begin()), higgsJet->fatjet.pt(),
				   higgsJet->fatjet.eta(), higgsJet->fatjet.phi(), numFiltJets,
				   numBtagFiltJets);
		}

		if (numBtagFiltJets >= 2)
			++local.n_Htags;
	}

	return 0;
}

int CU_ttH_EDA::Top_tagger(Handle<boosted::HTTTopJetCollection> top_jets,
						   CU_ttH_EDA_event_vars &local)
{
	local.n_ttags = 0;

	if (!top_jets.isValid())
		return 1;

	boosted::HTTTopJetCollection heptopjets =
		BoostedUtils::GetSortedByPt(*top_jets);

	for (boosted::HTTTopJetCollection::iterator topJet = heptopjets.begin();
		 topJet != heptopjets.end(); ++topJet) {
		// pt and eta requirements on top jet
		if (topJet->fatjet.pt() <= 250. || abs(topJet->fatjet.eta()) >= 1.8)
			continue;

		// pt and eta requirements on subjets
		if (topJet->nonW.pt() <= 20 || abs(topJet->nonW.eta()) >= 2.5 ||
			topJet->W1.pt() <= 20 || abs(topJet->W1.eta()) >= 2.5 ||
			topJet->W2.pt() <= 20 || abs(topJet->W2.eta()) >= 2.5)
			continue;

		// must be top-tagged
		if (toptagger.GetTopTaggerOutput(*topJet)<=-1) 
			continue;

		++local.n_ttags;
	}

	return 0;
}


///
void CU_ttH_EDA::Fill_Tau_Eff_Hist(CU_ttH_EDA_gen_vars &gen,
					   CU_ttH_EDA_event_vars &local)
{
	// Note: |eta| <= 2.3  pT > 20 GeV cut on genTau
	int nGenTau = 0;
	for (size_t i = 0; i < gen.x_daughters.size(); ++i) {
		if ( abs(gen.x_daughters[i].pdgId()) != 15)  // check if it is tau
			continue;
		if (gen.tau_class[i] != 1)  // check if it is hadronic
			continue;
		if (gen.x_daughters[i].pt() < 20
			or abs(gen.x_daughters[i].eta()) > 2.3)  // pT and eta cuts
			continue;
		
		++nGenTau;
		h_genHadTau_pt  -> Fill(gen.x_daughters[i].pt());
		h_genHadTau_eta -> Fill(gen.x_daughters[i].eta());
		h_genHadTau_phi -> Fill(gen.x_daughters[i].phi());
	}
	h_num_genHadTau -> Fill(nGenTau);

	int nTau_noniso = 0;
	for ( auto & itau : local.noniso_tau_selected ) {
		const reco::GenParticle* genTau = getGenTau(itau);
		if (!genTau) continue;  // check if there is associated genTau
		if ( genTau->pt() < 20 or abs(genTau->eta()) > 2.3 )
			continue;

		++nTau_noniso;
		h_selectedTau_noniso_genpt  -> Fill(genTau->pt());
		h_selectedTau_noniso_geneta -> Fill(genTau->eta());
		h_selectedTau_noniso_genphi -> Fill(genTau->phi());
	}
	h_num_selectedTau_noniso -> Fill(nTau_noniso);

	int nTau_loose = 0;
	for ( auto & itau : local.loose_tau_selected ) {
		const reco::GenParticle* genTau = getGenTau(itau);
		if (!genTau) continue;  // check if there is associated genTau
		if ( genTau->pt() < 20 or abs(genTau->eta()) > 2.3 )
			continue;

		++nTau_loose;
		h_selectedTau_loose_genpt  -> Fill(genTau->pt());
		h_selectedTau_loose_geneta -> Fill(genTau->eta());
		h_selectedTau_loose_genphi -> Fill(genTau->phi());
	}
	h_num_selectedTau_loose -> Fill(nTau_loose);
	
	int nTau_medium = 0;
	for ( auto & itau : local.medium_tau_selected ) {
		const reco::GenParticle* genTau = getGenTau(itau);
		if (!genTau) continue;  // check if there is associated genTau
		if ( genTau->pt() < 20 or abs(genTau->eta()) > 2.3 )
			continue;

		++nTau_medium;
		h_selectedTau_medium_genpt  -> Fill(genTau->pt());
		h_selectedTau_medium_geneta -> Fill(genTau->eta());
		h_selectedTau_medium_genphi -> Fill(genTau->phi());
	}
	h_num_selectedTau_medium -> Fill(nTau_medium);
	
	int nTau_tight = 0;
	for ( auto & itau : local.tight_tau_selected ) {
		const reco::GenParticle* genTau = getGenTau(itau);
		if (!genTau) continue;  // check if there is associated genTau
		if ( genTau->pt() < 20 or abs(genTau->eta()) > 2.3 )
			continue;

		++nTau_tight;
		h_selectedTau_tight_genpt  -> Fill(genTau->pt());
		h_selectedTau_tight_geneta -> Fill(genTau->eta());
		h_selectedTau_tight_genphi -> Fill(genTau->phi());
	}
	h_num_selectedTau_tight -> Fill(nTau_tight);
}

const reco::GenParticle* CU_ttH_EDA::getGenTau(const pat::Tau& patTau) {
	
	std::vector<reco::GenParticleRef> associatedGenParticles = patTau.genParticleRefs();

	//From Pat::Tau tutorial	
	for (std::vector<reco::GenParticleRef>::const_iterator igen = associatedGenParticles.begin();
		 igen != associatedGenParticles.end(); ++igen) {
		if ( igen->isAvailable() ) {
			const reco::GenParticleRef& genParticle = (*igen);
			if ( abs(genParticle->pdgId()) == 15 )
				return genParticle.get();;
		}
	}

	return 0;
}

///
void CU_ttH_EDA::Get_GenInfo(Handle<reco::GenParticleCollection> pruned,
							 Handle<pat::PackedGenParticleCollection> packed,
							 CU_ttH_EDA_gen_vars &gen)
{
	for (size_t i = 0; i < pruned->size(); ++i) {

		const reco::GenParticle *p = &(*pruned)[i];
		int ndaugs = p->numberOfDaughters();

		if (ndaugs == 0)
			continue;
		
		const reco::Candidate *d = p->daughter(0);
		if (p->pdgId() == d->pdgId())
			continue;

		// -- top --
		if (abs(p->pdgId()) == 6) {

			gen.tops.push_back(*p);

			// -- top daughters --
			for (size_t j = 0; j < (unsigned)ndaugs; ++j) {
				const reco::Candidate *immed_tdaug = p->daughter(j);
				const reco::Candidate *tdaug = get_last_in_decay_chain(immed_tdaug);
				// TO BE CHECKED: keep first or last b in the chain?		
				gen.top_daughters.push_back(*tdaug);

				// -- w daughters --
				if (abs(tdaug->pdgId()) == 24) {
					
					int ndaugs_w = tdaug->numberOfDaughters();
					
					for (size_t k = 0; k < (unsigned)ndaugs_w; ++k) {
						const reco::Candidate *wdaug = tdaug->daughter(k);
						gen.w_daughters.push_back(*wdaug);
					} // end of w daughter loop
				}

			} // end of top daughter loop

		}
		// -- mediator --
		else if (p->pdgId() == 25 /*|| p->pdgId()== 21 || p->pdgId() == 22 || p->pdgId() == 23*/) {

			// ---------------------------------------------------------------------
			// print decay chain
			// bool printDetails = false;
			// if (p->pdgId() == 21 || p->pdgId() == 23) {
			//	int tmp_index = 0;
			//	cout <<"index"<< "\t"<< "id"<< "\t"<< "stat"<< "\t"<< "mother"
			//		 << "\t"<< "nDaug";
			//	if (printDetails)
			//		cout << "\t"<< "pt"<< "\t"<< "eta"<< "\t"<< "phi" << endl;
			//	else
			//		cout << endl;
			//
			//	cout << tmp_index <<"\t"<< p->pdgId() << "\t" << p->status()
			//		 <<"\t"<< "n/a"<<"\t"<< p->numberOfDaughters();
			//	if (printDetails)
			//		cout << "\t" << p->pt() << "\t" << p->eta() << "\t"
			//			 << p->phi() << endl;
			//	else
			//		cout << endl;
			//
			//	printDecayChain(*p, tmp_index, 0, printDetails);
			//}
			// ---------------------------------------------------------------------

			gen.x.push_back(*p);

			// -- mediator daughters --
			for (size_t j = 0; j < (unsigned)ndaugs; ++j) {
				const reco::Candidate *immed_xdaug = p->daughter(j);
				const reco::Candidate *xdaug = get_last_in_decay_chain(immed_xdaug);
				gen.x_daughters.push_back(*xdaug);

				// -- tau and tau daughers --
				if (abs(xdaug->pdgId()) == 15) {
					//int ndaugs_tau = xdaug->numberOfDaughters();
					std::vector<const reco::Candidate*> stabledaughters;
					get_stable_daughters(*xdaug,stabledaughters);
					gen.tau_class.push_back( tau_classifier(stabledaughters) );
				}
				else
					gen.tau_class.push_back(-1);  // not applicable
				
			} // end of mediator daughter loop

		} else
			continue;

	} // end of pruned loop
}

void CU_ttH_EDA::printDecayChain(const reco::Candidate &p, int &index,
								 int mother_index, bool details)
{
	int ndaug = p.numberOfDaughters();

	for (int j = 0; j < ndaug; ++j) {
		const reco::Candidate *daug = p.daughter(j);
		index++;
		cout << index << "\t" << daug->pdgId() << "\t" << daug->status() << "\t"
			 << mother_index << "\t" << daug->numberOfDaughters();
		if (details) {
			cout << "\t" << daug->pt() << "\t" << daug->eta() << "\t"
				 << daug->phi() << endl;
		} else
			cout << endl;

		if (daug->status() != 1)
			printDecayChain(*daug, index, index, details);
	}
}

const reco::Candidate* CU_ttH_EDA::get_last_in_decay_chain(const reco::Candidate* p)
{
	int id = p->pdgId();
	int ndaug = p->numberOfDaughters();
	bool decay_to_itself = false;
	int isame = -99;

	for (int j = 0; j < ndaug; ++j) {
 		const reco::Candidate* daug = p->daughter(j);
		if (daug -> pdgId() == id) {
			decay_to_itself = true;
			isame = j;
		}
	}

	if (!decay_to_itself) {
		const reco::Candidate* last_p = p;
		return last_p;
	}
	else {	
		const reco::Candidate* next_p = p->daughter(isame);
		return get_last_in_decay_chain(next_p);
	}
}

void CU_ttH_EDA::get_stable_daughters(const reco::Candidate& p,
									 std::vector<const reco::Candidate*>& stabledaughters)
{
	int ndaug = p.numberOfDaughters();

	for(int j = 0; j < ndaug; ++j) {
		const reco::Candidate* daug = p.daughter(j);
		if (daug->status() == 1) {
			stabledaughters.push_back(daug);
		}
		else
			get_stable_daughters(*daug,stabledaughters);
	}
}

int CU_ttH_EDA::tau_classifier(std::vector<const reco::Candidate*>& stabledaughters)
{
	bool found_lepton = false;
	for (auto & sdaug: stabledaughters) {
		if( abs(sdaug->pdgId())==11 or abs(sdaug->pdgId())==13 )
			found_lepton = true;
	}

	if (found_lepton)
		return 0;  // leptonic
	else
		return 1;  // hadronic
}

//
void CU_ttH_EDA::Make_Ntuple(CU_ttH_EDA_gen_vars &gen, CU_ttH_EDA_event_vars &local, TTree *eventTree)
{
	// clear variables
		
	n_electrons = -99;
	n_muons = -99;
	n_loose_taus = -99;
	n_loose_taus = -99;
	n_tight_taus = -99;
	n_jets = -99;
	n_btags = -99;

	pv_x = -99.9;
	pv_y = -99.9;
	pv_z = -99.9;
	
	//electrons.clear();
	e_pt.clear();
	e_eta.clear();
	e_phi.clear();
	e_mass.clear();
	e_charges.clear();
	e_vtx_dz.clear();
	e_vtx_dxy.clear();
	e_vz.clear();
	e_vx.clear();
	e_vy.clear();

	//muons.clear();
	mu_pt.clear();
	mu_eta.clear();
	mu_phi.clear();
	mu_mass.clear();
	mu_charges.clear();
	mu_vtx_dz.clear();
	mu_vtx_dxy.clear();
	mu_vz.clear();
	mu_vx.clear();
	mu_vy.clear();

	//loose_taus.clear();
	loose_tau_pt.clear();
	loose_tau_eta.clear();
	loose_tau_phi.clear();
	loose_tau_mass.clear();
	Ltau_charges.clear();
	Ltau_vtx_dz.clear();
	Ltau_vtx_dxy.clear();
	Ltau_vx.clear();
	Ltau_vy.clear();
	Ltau_vz.clear();
	Ltau_decaymode.clear();
	Ltau_ntrk.clear();
	//medium_taus.clear();
	medium_tau_pt.clear();
	medium_tau_eta.clear();
	medium_tau_phi.clear();
	medium_tau_mass.clear();
	Mtau_charges.clear();
	Mtau_vtx_dz.clear();
	Mtau_vtx_dxy.clear();
	Mtau_vx.clear();
	Mtau_vy.clear();
	Mtau_vz.clear();
	Mtau_decaymode.clear();
	Mtau_ntrk.clear();
	//tight_taus.clear();
	tight_tau_pt.clear();
	tight_tau_eta.clear();
	tight_tau_phi.clear();
	tight_tau_mass.clear();
	Ttau_charges.clear();
	Ttau_vtx_dz.clear();
	Ttau_vtx_dxy.clear();
	Ttau_vx.clear();
	Ttau_vy.clear();
	Ttau_vz.clear();
	Ttau_decaymode.clear();
	Ttau_ntrk.clear();

	//jets.clear();
	jet_pt.clear();
	jet_eta.clear();
	jet_phi.clear();
	jet_mass.clear();
	jet_charges.clear();
	jet_vz.clear();
	jet_vx.clear();
	jet_vy.clear();
	
	//bjets.clear();
	bjet_pt.clear();
	bjet_eta.clear();
	bjet_phi.clear();
	bjet_mass.clear();
	bjet_charges.clear();
	bjet_vz.clear();
	bjet_vx.clear();
	bjet_vy.clear();

	gen_x_pdgId.clear();
	gen_x_status.clear();
	//gen_x.clear();
	gen_x_pt.clear();
	gen_x_eta.clear();
	gen_x_phi.clear();
	gen_x_mass.clear();
	gen_x_vx.clear();
	gen_x_vy.clear();
	gen_x_vz.clear();
	
	gen_top_pdgId.clear();
	gen_top_status.clear();
	//gen_top.clear();
	gen_top_pt.clear();
	gen_top_eta.clear();
	gen_top_phi.clear();
	gen_top_mass.clear();
	gen_top_vx.clear();
	gen_top_vy.clear();
	gen_top_vz.clear();
	
	gen_xDaug_pdgId.clear();
	gen_xDaug_status.clear();
	//gen_xDaug.clear();
	gen_xDaug_pt.clear();
	gen_xDaug_eta.clear();
	gen_xDaug_phi.clear();
	gen_xDaug_mass.clear();
	gen_xDaug_vx.clear();
	gen_xDaug_vy.clear();
	gen_xDaug_vz.clear();
	
	gen_tau_class.clear();

	gen_topDaug_pdgId.clear();
	gen_topDaug_status.clear();
	//gen_topDaug.clear();
	gen_topDaug_pt.clear();
	gen_topDaug_eta.clear();
	gen_topDaug_phi.clear();
	gen_topDaug_vx.clear();
	gen_topDaug_vy.clear();
	gen_topDaug_vz.clear();
	
	gen_wDaug_pdgId.clear();
	gen_wDaug_status.clear();
	//gen_wDaug.clear();
	gen_wDaug_pt.clear();
	gen_wDaug_eta.clear();
	gen_wDaug_phi.clear();
	gen_wDaug_mass.clear();
	gen_wDaug_vx.clear();
	gen_wDaug_vy.clear();
	gen_wDaug_vz.clear();

	// Number of tags per event
	n_electrons = local.n_electrons;
	n_muons = local.n_muons;
	n_loose_taus = local.n_loose_taus;
	n_medium_taus = local.n_medium_taus;
	n_tight_taus = local.n_tight_taus;
	n_jets = local.n_jets;
	n_btags = local.n_btags;

	pv_x = pv.x();
	pv_y = pv.y();
	pv_z = pv.z();
	
	// electrons
	for (auto & ele : local.e_selected_sorted) {
		e_pt.push_back(ele.pt());
		e_eta.push_back(ele.eta());
		e_phi.push_back(ele.phi());
		e_mass.push_back(ele.mass());
		e_charges.push_back(ele.charge());
		if ( ele.gsfTrack().isAvailable() ) {
			e_vtx_dz.push_back( ele.gsfTrack()->dz(pv.position()) );
			e_vtx_dxy.push_back( ele.gsfTrack()->dxy(pv.position()) );
		}
		//e_vx.push_back(ele.vx());
		//e_vy.push_back(ele.vy());
		//e_vz.push_back(ele.vz());
		e_isGsfCtfScPixChargeConsistent.push_back(ele.isGsfCtfScPixChargeConsistent());
	}

	// muons
	for (auto & mu : local.mu_selected_sorted) {
		mu_pt.push_back(mu.pt());
		mu_eta.push_back(mu.eta());
		mu_phi.push_back(mu.phi());
		mu_mass.push_back(mu.mass());
		mu_charges.push_back(mu.charge());
		if ( mu.muonBestTrack().isAvailable() ) {
			// innerTrack? GlobalTrack?
			mu_vtx_dz.push_back( mu.muonBestTrack()->dz(pv.position()) );
			mu_vtx_dxy.push_back( mu.muonBestTrack()->dxy(pv.position()) );
			mu_relTrkPtError.push_back( mu.muonBestTrack()->ptError() /
										  mu.muonBestTrack()->pt());
		}
		//mu_vx.push_back(mu.vx());
		//mu_vy.push_back(mu.vy());
		//mu_vz.push_back(mu.vz());
	}

	// taus
	for (auto & tau : local.loose_tau_selected_sorted) {
	    loose_tau_pt.push_back(tau.pt());
		loose_tau_eta.push_back(tau.eta());
		loose_tau_phi.push_back(tau.phi());
		loose_tau_mass.push_back(tau.mass());
		Ltau_charges.push_back(tau.charge());
		Ltau_decaymode.push_back(tau.decayMode());
		//Ltau_ntrk.push_back(tau.signalTracks().size());
		// tau vertex
		//Ltau_vx.push_back(tau.vx());
		//Ltau_vy.push_back(tau.vy());
		//Ltau_vz.push_back(tau.vz());
		//Ltau_vr.push_back( sqrt(tau.vx()*tau.vx()+tau.vy()*tau.vy()) );
	}

	for (auto & tau : local.medium_tau_selected_sorted) {
	    medium_tau_pt.push_back(tau.pt());
		medium_tau_eta.push_back(tau.eta());
		medium_tau_phi.push_back(tau.phi());
		medium_tau_mass.push_back(tau.mass());
		Mtau_charges.push_back(tau.charge());
		Mtau_decaymode.push_back(tau.decayMode());
		//Mtau_ntrk.push_back(tau.signalTracks().size());
		//Mtau_vx.push_back(tau.vx());
		//Mtau_vy.push_back(tau.vy());
		//Mtau_vz.push_back(tau.vz());
		//Mtau_vr.push_back( sqrt(tau.vx()*tau.vx()+tau.vy()*tau.vy()) );
	}

	for (auto & tau : local.tight_tau_selected_sorted) {
		tight_tau_pt.push_back(tau.pt());
		tight_tau_eta.push_back(tau.eta());
		tight_tau_phi.push_back(tau.phi());
		tight_tau_mass.push_back(tau.mass());
		Ttau_charges.push_back(tau.charge());
		Ttau_decaymode.push_back(tau.decayMode());
		//Ttau_ntrk.push_back(tau.signalTracks().size());
		//Ttau_vx.push_back(tau.vx());
		//Ttau_vy.push_back(tau.vy());
		//Ttau_vz.push_back(tau.vz());
		//Ttau_vr.push_back( sqrt(tau.vx()*tau.vx()+tau.vy()*tau.vy()) );
	}

	// jets
	for (auto & jet : local.jets_selected_sorted) {
		jet_pt.push_back(jet.pt());
		jet_eta.push_back(jet.eta());
		jet_phi.push_back(jet.phi());
		jet_mass.push_back(jet.mass());
		jet_charges.push_back(jet.jetCharge());
		// jet vertex
		//jet_vx.push_back(jet.vx());
		//jet_vy.push_back(jet.vy());
		//jet_vz.push_back(jet.vz());
		//jet_vr.push_back( sqrt(jet.vx()*jet.vx()+jet.vy()*jet.vy()) );
	}

	// b-jets
	for (auto & bjet : local.jets_selected_tag_sorted) {
		bjet_pt.push_back(bjet.pt());
		bjet_eta.push_back(bjet.eta());
		bjet_phi.push_back(bjet.phi());
		bjet_mass.push_back(bjet.mass());
		// b vertex
		//bjet_vx.push_back(bjet.vx());
		//bjet_vy.push_back(bjet.vy());
		//bjet_vz.push_back(bjet.vz());
		//bjet_vr.push_back( sqrt(bjet.vx()*bjet.vx()+bjet.vy()*bjet.vy()) );
	}

	// MET
	MET_x = local.MET_corrected.px();
	MET_y = local.MET_corrected.py();
	
	// GenParticle Information
	// mediators
	for (auto & mediator : gen.x) {
		gen_x_pdgId.push_back(mediator.pdgId());
		gen_x_status.push_back(mediator.status());
		gen_x_pt.push_back(mediator.pt());
		gen_x_eta.push_back(mediator.eta());
		gen_x_phi.push_back(mediator.phi());
		gen_x_mass.push_back(mediator.mass());
		//gen_x_vr.push_back(sqrt(mediator.vx()*mediator.vx()
		//						+mediator.vy()*mediator.vy()));
		gen_x_vx.push_back(mediator.vx());
		gen_x_vy.push_back(mediator.vy());
		gen_x_vz.push_back(mediator.vz());
	}

	// tops
	for (auto & top : gen.tops) {
		gen_top_pdgId.push_back(top.pdgId());
		gen_top_status.push_back(top.status());
		gen_top_pt.push_back(top.pt());
		gen_top_eta.push_back(top.eta());
		gen_top_phi.push_back(top.phi());
		gen_top_mass.push_back(top.mass());
		//gen_top_vr.push_back(sqrt(top.vx()*top.vx()+top.vy()*top.vy()));
		gen_top_vx.push_back(top.vx());
		gen_top_vy.push_back(top.vy());
		gen_top_vz.push_back(top.vz());
	}

	// mediator daughters
	for (auto & xdaug : gen.x_daughters) {
		gen_xDaug_pdgId.push_back(xdaug.pdgId());
		gen_xDaug_status.push_back(xdaug.status());
		gen_xDaug_pt.push_back(xdaug.pt());
		gen_xDaug_eta.push_back(xdaug.eta());
		gen_xDaug_phi.push_back(xdaug.phi());
		gen_xDaug_mass.push_back(xdaug.mass());
		//gen_xDaug_vr.push_back(sqrt(xdaug.vx()*xdaug.vx()+xdaug.vy()*xdaug.vy()));
		gen_xDaug_vx.push_back(xdaug.vx());
		gen_xDaug_vy.push_back(xdaug.vy());
		gen_xDaug_vz.push_back(xdaug.vz());
	}

	gen_tau_class = gen.tau_class;

	// top daughters
	for (auto & tdaug : gen.top_daughters) {
		gen_topDaug_pdgId.push_back(tdaug.pdgId());
		gen_topDaug_status.push_back(tdaug.status());
		gen_topDaug_pt.push_back(tdaug.pt());
		gen_topDaug_eta.push_back(tdaug.eta());
		gen_topDaug_phi.push_back(tdaug.phi());
		gen_topDaug_mass.push_back(tdaug.mass());
		//gen_topDaug_vr.push_back(sqrt(tdaug.vx()*tdaug.vx()+tdaug.vy()*tdaug.vy()));
		gen_topDaug_vx.push_back(tdaug.vx());
		gen_topDaug_vy.push_back(tdaug.vy());
		gen_topDaug_vz.push_back(tdaug.vz());
	}

	// w daughters
	for (auto & wdaug : gen.w_daughters) {
		gen_wDaug_pdgId.push_back(wdaug.pdgId());
		gen_wDaug_status.push_back(wdaug.status());
		gen_wDaug_pt.push_back(wdaug.pt());
		gen_wDaug_eta.push_back(wdaug.eta());
		gen_wDaug_phi.push_back(wdaug.phi());
		gen_wDaug_mass.push_back(wdaug.mass());
		//gen_wDaug_vr.push_back(sqrt(wdaug.vx()*wdaug.vx()+wdaug.vy()*wdaug.vy()));
		gen_wDaug_vx.push_back(wdaug.vx());
		gen_wDaug_vy.push_back(wdaug.vy());
		gen_wDaug_vz.push_back(wdaug.vz());
	}

	eventTree->Fill();
}


/// Other functions
void CU_ttH_EDA::Check_Fill_Print_muj(CU_ttH_EDA_event_vars &local)
{
	int fill_itr = 0;

	h_tth_syncex1_mu->Fill(0.5 + fill_itr++); // fills 0.5 first
	if (local.pass_single_mu) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_mu_cut1, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_muons == 1) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_mu_cut2, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_electrons == 0) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_mu_cut3, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_jets >= 4) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_mu_cut4, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_btags >= 2) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_mu_cut5, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_ttags >= 1) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_mu_cut6, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_Htags >= 1) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr);
		Print_event_in_file1(events_mu_cut7, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;
}

void CU_ttH_EDA::Check_Fill_Print_ej(CU_ttH_EDA_event_vars &local)
{
	int fill_itr = 0;

	h_tth_syncex1_ele->Fill(0.5 + fill_itr++); // fills 0.5 first

	if (local.pass_single_e) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_e_cut1, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_electrons == 1) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_e_cut2, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_muons == 0) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_e_cut3, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_jets >= 4) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_e_cut4, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_btags >= 2) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_e_cut5, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_ttags >= 1) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_e_cut6, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_Htags >= 1) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr);
		Print_event_in_file1(events_e_cut7, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;
}

// template<class lepton>
void CU_ttH_EDA::Check_Fill_Print_dimuj(CU_ttH_EDA_event_vars &local)
{
	std::vector<pat::Muon> muon1, muon2;
	if (local.n_muons >= 2) {
		muon1.push_back(local.mu_selected_sorted[0]);
		muon2.push_back(local.mu_selected_sorted[1]);

		TLorentzVector mu1, mu2;
		mu1.SetPtEtaPhiM(local.mu_selected_sorted[0].pt(),
						 local.mu_selected_sorted[0].eta(),
						 local.mu_selected_sorted[0].phi(),
						 local.mu_selected_sorted[0].mass());
		mu2.SetPtEtaPhiM(local.mu_selected_sorted[1].pt(),
						 local.mu_selected_sorted[1].eta(),
						 local.mu_selected_sorted[1].phi(),
						 local.mu_selected_sorted[1].mass());
		local.dimuon_mass = (mu1 + mu2).M();
	} else if (local.n_muons == 1) {
		muon1.push_back(local.mu_selected_sorted[0]);
	}

	int fill_itr = 0;

	h_tth_syncex1_dimu->Fill(0.5 + fill_itr++); // fills 0.5 first

	if (local.pass_double_mu) {
		h_tth_syncex1_dimu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_dimu_cut1, muon1, muon2,
									  local.dimuon_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_muons >= 2) {
		if (local.mu_selected_sorted[0].charge() !=
			local.mu_selected_sorted[1].charge()) {
			h_tth_syncex1_dimu->Fill(0.5 + fill_itr++);
			Print_event_in_file1_dilepton(events_dimu_cut2, muon1, muon2,
										  local.dimuon_mass,
										  local.jets_selected_sorted, local);
		} else {
			std::cout << "Found event with two same-charge muons" << std::endl;
		}
		if (local.n_muons > 2) {
			std::cout << "Found event with more than 2 muons" << std::endl;
		}
	} else
		return;

	if (local.dimuon_mass >= 20) {
		h_tth_syncex1_dimu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_dimu_cut3, muon1, muon2,
									  local.dimuon_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.dimuon_mass <= 76 or local.dimuon_mass >= 106) {
		h_tth_syncex1_dimu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_dimu_cut4, muon1, muon2,
									  local.dimuon_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_jets >= 2) {
		h_tth_syncex1_dimu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_dimu_cut5, muon1, muon2,
									  local.dimuon_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.MET_corrected.pt() > 40) {
		h_tth_syncex1_dimu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_dimu_cut6, muon1, muon2,
									  local.dimuon_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_btags >= 1) {
		h_tth_syncex1_dimu->Fill(0.5 + fill_itr);
		Print_event_in_file1_dilepton(events_dimu_cut7, muon1, muon2,
									  local.dimuon_mass,
									  local.jets_selected_sorted, local);
	} else
		return;
}

void CU_ttH_EDA::Check_Fill_Print_dielej(CU_ttH_EDA_event_vars &local)
{
	std::vector<pat::Electron> electron1, electron2;
	if (local.n_electrons >= 2) {
		electron1.push_back(local.e_selected_sorted[0]);
		electron2.push_back(local.e_selected_sorted[1]);

		TLorentzVector ele1, ele2;
		ele1.SetPtEtaPhiM(local.e_selected_sorted[0].pt(),
						  local.e_selected_sorted[0].eta(),
						  local.e_selected_sorted[0].phi(),
						  local.e_selected_sorted[0].mass());
		ele2.SetPtEtaPhiM(local.e_selected_sorted[1].pt(),
						  local.e_selected_sorted[1].eta(),
						  local.e_selected_sorted[1].phi(),
						  local.e_selected_sorted[1].mass());
		local.dielectron_mass = (ele1 + ele2).M();
	} else if (local.n_electrons == 1) {
		electron1.push_back(local.e_selected_sorted[0]);
	}

	int fill_itr = 0;

	h_tth_syncex1_diele->Fill(0.5 + fill_itr++); // fills 0.5 first

	if (local.pass_double_e) {
		h_tth_syncex1_diele->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_diele_cut1, electron1, electron2,
									  local.dielectron_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_electrons >= 2) {
		if (local.e_selected_sorted[0].charge() !=
			local.e_selected_sorted[1].charge()) {
			h_tth_syncex1_diele->Fill(0.5 + fill_itr++);
			Print_event_in_file1_dilepton(events_diele_cut2, electron1,
										  electron2, local.dielectron_mass,
										  local.jets_selected_sorted, local);
		} else {
			std::cout << "Found event with two same-charge electrons"
					  << std::endl;
		}
		if (local.n_electrons > 2) {
			std::cout << "Found event with more than 2 electrons" << std::endl;
		}
	} else
		return;

	if (local.dielectron_mass >= 20) {
		h_tth_syncex1_diele->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_diele_cut3, electron1, electron2,
									  local.dielectron_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.dielectron_mass <= 76 or local.dielectron_mass >= 106) {
		h_tth_syncex1_diele->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_diele_cut4, electron1, electron2,
									  local.dielectron_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_jets >= 2) {
		h_tth_syncex1_diele->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_diele_cut5, electron1, electron2,
									  local.dielectron_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.MET_corrected.pt() > 40) {
		h_tth_syncex1_diele->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_diele_cut6, electron1, electron2,
									  local.dielectron_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_btags >= 1) {
		h_tth_syncex1_diele->Fill(0.5 + fill_itr);
		Print_event_in_file1_dilepton(events_diele_cut7, electron1, electron2,
									  local.dielectron_mass,
									  local.jets_selected_sorted, local);
	} else
		return;
}

void CU_ttH_EDA::Check_Fill_Print_elemuj(CU_ttH_EDA_event_vars &local)
{
	std::vector<pat::Electron> electron1;
	std::vector<pat::Muon> muon1;
	if (local.n_electrons >= 1) {
		electron1.push_back(local.e_selected_sorted[0]);
	}
	if (local.n_muons >= 1) {
		muon1.push_back(local.mu_selected_sorted[0]);
	}
	if ((local.n_electrons >= 1) and (local.n_muons >= 1)) {
		// electron1.push_back( local.e_selected_sorted[0] );
		// muon1.push_back( local.mu_selected_sorted[0] );

		TLorentzVector ele1, mu1;
		ele1.SetPtEtaPhiM(local.e_selected_sorted[0].pt(),
						  local.e_selected_sorted[0].eta(),
						  local.e_selected_sorted[0].phi(),
						  local.e_selected_sorted[0].mass());
		mu1.SetPtEtaPhiM(local.mu_selected_sorted[0].pt(),
						 local.mu_selected_sorted[0].eta(),
						 local.mu_selected_sorted[0].phi(),
						 local.mu_selected_sorted[0].mass());
		local.dilepton_mass = (ele1 + mu1).M();
	}

	int fill_itr = 0;

	h_tth_syncex1_elemu->Fill(0.5 + fill_itr++); // fills 0.5 first

	if (local.pass_elemu) {
		h_tth_syncex1_elemu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_elemu_cut1, muon1, electron1,
									  local.dilepton_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if ((local.n_electrons >= 1) and (local.n_muons >= 1)) {
		if (local.e_selected_sorted[0].charge() !=
			local.mu_selected_sorted[0].charge()) {
			h_tth_syncex1_elemu->Fill(0.5 + fill_itr++);
			Print_event_in_file1_dilepton(events_elemu_cut2, muon1, electron1,
										  local.dilepton_mass,
										  local.jets_selected_sorted, local);
		} else {
			std::cout << "Found event with two same-charge leptons"
					  << std::endl;
		}
		if (local.n_electrons > 2) {
			std::cout << "Found event with more than 2 leptons" << std::endl;
		}
	} else
		return;

	if (local.dilepton_mass >= 20) {
		h_tth_syncex1_elemu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_elemu_cut3, muon1, electron1,
									  local.dilepton_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_jets >= 2) {
		h_tth_syncex1_elemu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_elemu_cut4, muon1, electron1,
									  local.dilepton_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_btags >= 1) {
		h_tth_syncex1_elemu->Fill(0.5 + fill_itr);
		Print_event_in_file1_dilepton(events_elemu_cut5, muon1, electron1,
									  local.dilepton_mass,
									  local.jets_selected_sorted, local);
	} else
		return;
}

template <class lepton>
int CU_ttH_EDA::Print_event_in_file1(FILE *file, lepton &lpt,
									 std::vector<pat::Jet> &jets,
									 CU_ttH_EDA_event_vars &local)
{
	// print generic event info
	fprintf(file, "%6d %8d %10d   ", local.run_nr, local.lumisection_nr,
			local.event_nr);

	// print lepton info
	if (lpt.size() != 0)
		fprintf(file, "%6.2f %+4.2f %+4.2f   ", lpt[0].pt(), lpt[0].eta(),
				lpt[0].phi());
	else
		fprintf(file, "%6.2f %+4.2f %+4.2f   ", -99., -99., -99.);

	// print jet(s) info
	if (jets.size() >= 4) {
		fprintf(file,
				"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
				jets[0].pt(), jets[1].pt(), jets[2].pt(), jets[3].pt(),
				jets[0].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"),
				jets[1].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"),
				jets[2].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"),
				jets[3].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"));
	} else {
		switch (jets.size()) {
		case 3:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					jets[0].pt(), jets[1].pt(), jets[2].pt(), -99.,
					jets[0].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					jets[1].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					jets[2].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					-99.);
			break;

		case 2:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					jets[0].pt(), jets[1].pt(), -99., -99.,
					jets[0].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					jets[1].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					-99., -99.);
			break;

		case 1:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					jets[0].pt(), -99., -99., -99.,
					jets[0].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					-99., -99., -99.);
			break;

		default:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					-99., -99., -99., -99., -99., -99., -99., -99.);
		}
	}

	// print number of tags
	fprintf(file, "%2d  %2d   %2d  %2d\n", local.n_jets, local.n_btags,
			local.n_ttags, local.n_Htags);

	return ferror(file);
}

template <class lep1, class lep2>
// template<class lepton2>
int CU_ttH_EDA::Print_event_in_file1_dilepton(FILE *file, lep1 &lepton1,
											  lep2 &lepton2,
											  double dilepton_mass,
											  std::vector<pat::Jet> &jets,
											  CU_ttH_EDA_event_vars &local)
{
	// std::vector<lepton> leptons;
	// switch (dilepton_type) {
	//        case 'dimuon':
	//        leptons

	// print generic event info
	fprintf(file, "%6d %8d %10d   ", local.run_nr, local.lumisection_nr,
			local.event_nr);

	// print lepton info
	if (lepton1.size() == 1) {
		fprintf(file, "%6.2f %+4.2f %+4.2f   ", lepton1[0].pt(),
				lepton1[0].eta(), lepton1[0].phi());
	} else {
		fprintf(file, "%6.2f %+4.2f %+4.2f   ", -99., -99., -99.);
	}

	if (lepton2.size() == 1) {
		fprintf(file, "%6.2f %+4.2f %+4.2f   ", lepton2[0].pt(),
				lepton2[0].eta(), lepton2[0].phi());
	} else {
		fprintf(file, "%6.2f %+4.2f %+4.2f   ", -99., -99., -99.);
	}

	// print MET and dilepton mass
	fprintf(file, "%6.2f  %6.2f   ", dilepton_mass, local.MET_corrected.pt());

	// print jet(s) info
	if (jets.size() >= 4) {
		fprintf(file,
				"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
				jets[0].pt(), jets[1].pt(), jets[2].pt(), jets[3].pt(),
				jets[0].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"),
				jets[1].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"),
				jets[2].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"),
				jets[3].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"));
	} else {
		switch (jets.size()) {
		case 3:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					jets[0].pt(), jets[1].pt(), jets[2].pt(), -99.,
					jets[0].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					jets[1].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					jets[2].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					-99.);
			break;

		case 2:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					jets[0].pt(), jets[1].pt(), -99., -99.,
					jets[0].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					jets[1].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					-99., -99.);
			break;

		case 1:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					jets[0].pt(), -99., -99., -99.,
					jets[0].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					-99., -99., -99.);
			break;

		default:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					-99., -99., -99., -99., -99., -99., -99., -99.);
		}
	}

	// print number of tags
	fprintf(file, "%2d  %2d\n", local.n_jets, local.n_btags);

	return ferror(file);
}

///

#endif
