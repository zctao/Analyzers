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
	edm::Handle<reco::VertexCollection> vertices)
{
	/// Primary vertex handling
	if (!vertices.isValid())
		return 1;

	reco::Vertex vertex;
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

void CU_ttH_EDA::Check_Fill_Print_dimutauh(CU_ttH_EDA_event_vars &local)
{
	int fill_itr = 0;

	h_tth_syncex_dimutauh->Fill(0.5 + fill_itr++); // fill 0.5 first

	if (local.pass_double_mu)
		h_tth_syncex_dimutauh->Fill(0.5 + fill_itr++);
	else
		return;

	bool SameSignLeptons = false;
	if (local.n_muons >= 2) {
		
		for (int i=0; i<local.n_muons; ++i) {
			for (int j=i+1; j<local.n_muons; ++j) {
				if (local.mu_selected_sorted[i].charge() ==
					local.mu_selected_sorted[j].charge())
					SameSignLeptons = true;
			}
		}

		if (SameSignLeptons)
			h_tth_syncex_dimutauh->Fill(0.5 + fill_itr++);
		else
			return;	
	}	
	else
		return;

	if (local.n_taus >= 1)
		h_tth_syncex_dimutauh->Fill(0.5 + fill_itr++);
	else
		return;
	
	if (local.n_jets >= min_njets)
		h_tth_syncex_dimutauh->Fill(0.5 + fill_itr++);
	else
		return;

	if (local.n_btags >= min_nbtags)
		h_tth_syncex_dimutauh->Fill(0.5 + fill_itr++);
	else
		return;	
}

void CU_ttH_EDA::Check_Fill_Print_dieletauh(CU_ttH_EDA_event_vars &local)
{
	int fill_itr = 0;

	h_tth_syncex_dieletauh->Fill(0.5 + fill_itr++); // fill 0.5 first

	if (local.pass_double_e)
		h_tth_syncex_dieletauh->Fill(0.5 + fill_itr++);
	else
		return;

	bool SameSignLeptons = false;
	if (local.n_electrons >= 2) {
		
		for (int i=0; i<local.n_electrons; ++i) {
			for (int j=i+1; j<local.n_electrons; ++j) {
				if (local.e_selected_sorted[i].charge() ==
					local.e_selected_sorted[j].charge())
					SameSignLeptons = true;
			}
		}

		if (SameSignLeptons)
			h_tth_syncex_dieletauh->Fill(0.5 + fill_itr++);
		else
			return;	
	}	
	else
		return;

	if (local.n_taus >= 1)
		h_tth_syncex_dieletauh->Fill(0.5 + fill_itr++);
	else
		return;
	
	if (local.n_jets >= min_njets)
		h_tth_syncex_dieletauh->Fill(0.5 + fill_itr++);
	else
		return;

	if (local.n_btags >= min_nbtags)
		h_tth_syncex_dieletauh->Fill(0.5 + fill_itr++);
	else
		return;	
}

void CU_ttH_EDA::Check_Fill_Print_elemutauh(CU_ttH_EDA_event_vars &local)
{
	int fill_itr = 0;

	h_tth_syncex_elemutauh->Fill(0.5 + fill_itr++); // fill 0.5 first

	if (local.pass_elemu)
		h_tth_syncex_elemutauh->Fill(0.5 + fill_itr++);
	else
		return;

	bool SameSignLeptons = false;
	if (local.n_electrons >= 1 and local.n_muons >= 1) {

		for (int i=0; i<local.n_electrons; ++i) {
			for (int j=0; j<local.n_muons; ++j) {
				if (local.e_selected_sorted[i].charge() ==
					local.mu_selected_sorted[j].charge())
					SameSignLeptons = true;
			}
		}

		if (SameSignLeptons)
			h_tth_syncex_elemutauh->Fill(0.5 + fill_itr++);
		else
			return;
	}
	else
		return;

	if (local.n_taus >= 1)
		h_tth_syncex_elemutauh->Fill(0.5 + fill_itr++);
	else
		return;
	
	if (local.n_jets >= min_njets)
		h_tth_syncex_elemutauh->Fill(0.5 + fill_itr++);
	else
		return;

	if (local.n_btags >= min_nbtags)
		h_tth_syncex_elemutauh->Fill(0.5 + fill_itr++);
	else
		return;
}

void CU_ttH_EDA::Check_Fill_Print_eleditauh(CU_ttH_EDA_event_vars &local) {
	int fill_itr = 0;

	h_tth_syncex_eleditauh->Fill(0.5 + fill_itr++); // fill 0.5 first

	if (local.pass_single_e)
		h_tth_syncex_eleditauh->Fill(0.5 + fill_itr++);
	else
		return;

	if (local.n_electrons >= 1)
		h_tth_syncex_eleditauh->Fill(0.5 + fill_itr++);
	else
		return;

	if (local.n_taus >= 2)
		h_tth_syncex_eleditauh->Fill(0.5 + fill_itr++);
	else
		return;

	double mTT = -99;
	TLorentzVector tau1,tau2;
	tau1.SetPtEtaPhiM(local.tau_selected_sorted[0].pt(),
					  local.tau_selected_sorted[0].eta(),
					  local.tau_selected_sorted[0].phi(),
					  local.tau_selected_sorted[0].mass());
	tau2.SetPtEtaPhiM(local.tau_selected_sorted[1].pt(),
					  local.tau_selected_sorted[1].eta(),
					  local.tau_selected_sorted[1].phi(),
					  local.tau_selected_sorted[1].mass());
	mTT = (tau1+tau2).M();
	if (mTT > 120 and mTT < 130)
		h_tth_syncex_eleditauh->Fill(0.5 + fill_itr++);
	else
		return;

	if (local.n_jets >= min_njets)
		h_tth_syncex_eleditauh->Fill(0.5 + fill_itr++);
	else
		return;

	if (local.n_btags >= min_nbtags)
		h_tth_syncex_eleditauh->Fill(0.5 + fill_itr++);
	else
		return;
}

void CU_ttH_EDA::Check_Fill_Print_muditauh(CU_ttH_EDA_event_vars &local) {
	int fill_itr = 0;

	h_tth_syncex_muditauh->Fill(0.5 + fill_itr++); // fill 0.5 first

	if (local.pass_single_mu)
		h_tth_syncex_muditauh->Fill(0.5 + fill_itr++);
	else
		return;

	// need to check if necessary
	if (local.n_muons >= 1)
		h_tth_syncex_muditauh->Fill(0.5 + fill_itr++);
	else
		return;

	if (local.n_taus >= 2)
		h_tth_syncex_muditauh->Fill(0.5 + fill_itr++);
	else
		return;

	double mTT = -99;
	TLorentzVector tau1,tau2;
	tau1.SetPtEtaPhiM(local.tau_selected_sorted[0].pt(),
					  local.tau_selected_sorted[0].eta(),
					  local.tau_selected_sorted[0].phi(),
					  local.tau_selected_sorted[0].mass());
	tau2.SetPtEtaPhiM(local.tau_selected_sorted[1].pt(),
					  local.tau_selected_sorted[1].eta(),
					  local.tau_selected_sorted[1].phi(),
					  local.tau_selected_sorted[1].mass());
	mTT = (tau1+tau2).M();
	if (mTT > 120 and mTT < 130)
		h_tth_syncex_muditauh->Fill(0.5 + fill_itr++);
	else
		return;

	if (local.n_jets >= min_njets)
		h_tth_syncex_muditauh->Fill(0.5 + fill_itr++);
	else
		return;

	if (local.n_btags >= min_nbtags)
		h_tth_syncex_muditauh->Fill(0.5 + fill_itr++);
	else
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
				const reco::Candidate *tdaug = p->daughter(j);
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
		else if (p->pdgId() == 25 ||
				 p->pdgId()== 21 || p->pdgId() == 22 || p->pdgId() == 23) {

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
				const reco::Candidate *xdaug = p->daughter(j);
				gen.x_daughters.push_back(*xdaug);
			} // end of mediator daughter loop

		} else
			continue;

	} // end of pruned loop
}

void CU_ttH_EDA::Write_to_Tree(CU_ttH_EDA_gen_vars &gen, TTree *eventTree)
{
	// clear variables
	x_pdgId.clear();
	x_status.clear();
	x_pt.clear();
	x_eta.clear();
	x_phi.clear();
	x_mass.clear();

	top_pdgId.clear();
	top_status.clear();
	top_pt.clear();
	top_eta.clear();
	top_phi.clear();
	top_mass.clear();

	xDaug_pdgId.clear();
	xDaug_status.clear();
	xDaug_pt.clear();
	xDaug_eta.clear();
	xDaug_phi.clear();
	xDaug_mass.clear();

	topDaug_pdgId.clear();
	topDaug_status.clear();
	topDaug_pt.clear();
	topDaug_eta.clear();
	topDaug_phi.clear();
	topDaug_mass.clear();

	wDaug_pdgId.clear();
	wDaug_status.clear();
	wDaug_pt.clear();
	wDaug_eta.clear();
	wDaug_phi.clear();
	wDaug_mass.clear();

	// mediators
	for (size_t i = 0; i < gen.x.size(); ++i) {
		x_pdgId.push_back(gen.x[i].pdgId());
		x_status.push_back(gen.x[i].status());
		x_pt.push_back(gen.x[i].pt());
		x_eta.push_back(gen.x[i].eta());
		x_phi.push_back(gen.x[i].phi());
		x_mass.push_back(gen.x[i].mass());
	}

	// tops
	for (size_t i = 0; i < gen.tops.size(); ++i) {
		top_pdgId.push_back(gen.tops[i].pdgId());
		top_status.push_back(gen.tops[i].status());
		top_pt.push_back(gen.tops[i].pt());
		top_eta.push_back(gen.tops[i].eta());
		top_phi.push_back(gen.tops[i].phi());
		top_mass.push_back(gen.tops[i].mass());
	}

	// mediator daughters
	for (size_t i = 0; i < gen.x_daughters.size(); ++i) {
		xDaug_pdgId.push_back(gen.x_daughters[i].pdgId());
		xDaug_status.push_back(gen.x_daughters[i].status());
		xDaug_pt.push_back(gen.x_daughters[i].pt());
		xDaug_eta.push_back(gen.x_daughters[i].eta());
		xDaug_phi.push_back(gen.x_daughters[i].phi());
		xDaug_mass.push_back(gen.x_daughters[i].mass());
	}

	// top daughters
	for (size_t i = 0; i < gen.top_daughters.size(); ++i) {
		topDaug_pdgId.push_back(gen.top_daughters[i].pdgId());
		topDaug_status.push_back(gen.top_daughters[i].status());
		topDaug_pt.push_back(gen.top_daughters[i].pt());
		topDaug_eta.push_back(gen.top_daughters[i].eta());
		topDaug_phi.push_back(gen.top_daughters[i].phi());
		topDaug_mass.push_back(gen.top_daughters[i].mass());
	}

	// w daughters
	for (size_t i = 0; i < gen.w_daughters.size(); ++i) {
		wDaug_pdgId.push_back(gen.w_daughters[i].pdgId());
		wDaug_status.push_back(gen.w_daughters[i].status());
		wDaug_pt.push_back(gen.w_daughters[i].pt());
		wDaug_eta.push_back(gen.w_daughters[i].eta());
		wDaug_phi.push_back(gen.w_daughters[i].phi());
		wDaug_mass.push_back(gen.w_daughters[i].mass());
	}

	eventTree->Fill();
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

#endif
