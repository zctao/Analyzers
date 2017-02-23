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

unsigned int CU_ttH_EDA::getTriggerBits(edm::Handle<edm::TriggerResults> triggerResults, std::vector<std::string>& names, HLTConfigProvider& config)
{
	unsigned int bits = 0;
	unsigned int iname = 0;

	for (const auto & n : names) {
		unsigned int index = config.triggerIndex(n);

		if (index >= triggerResults->size()) {
			std::cerr << "Failed to find " << n << std::endl;
			continue;
		}

		if (triggerResults->accept(index))
			bits |= 1<< iname;

		++iname;
	}

	return bits;
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


/// Other functions

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

float CU_ttH_EDA::getMHT(CU_ttH_EDA_event_vars &local)
{
	float MHT_x = 0;
	float MHT_y = 0;

	for (auto & mu : local.mu_fakeable) {
		MHT_x -= mu.px();
		MHT_y -= mu.py();
	}

	for (auto & ele : local.e_fakeable) {
		MHT_x -= ele.px();
		MHT_y -= ele.py();
	}

	for (auto & tau : local.tau_selected) {
		MHT_x -= tau.px();
		MHT_y -= tau.py();
	}

	for (auto & jet : local.jets_selected_sorted) {
		MHT_x -= jet.px();
		MHT_y -= jet.py();
	}

	return sqrt(MHT_x * MHT_x + MHT_y * MHT_y);
}

bool CU_ttH_EDA::pass_event_sel_2l(CU_ttH_EDA_event_vars &local,
								   Selection_types selection_region,
								   int &ilep, int& ibtag)
{

	if (debug) {
		std::cout << std::endl;
		std::cout << "start event selection:" << std::endl;
	}
	
	//////////////////////////
    /// Lepton number
	// at least 2 fakeable leptons
	// no more than 2 tight leptons
	if (!(local.leptons_fakeable.size() >= 2 and local.leptons_tight.size()<=2)) {
		if (debug) {
			std::cout << "FAIL lepton number requirement" << std::endl;
		}
		return false;
	}
	// at least 1 selected tau
	// i.e. pass "byMediumIsolationMVArun2v1DBdR03oldDMwLT"
	if (local.tau_selected.size() < 1) {
		if (debug) {
			std::cout << "FAIL tau number requirement" << std::endl;
		}
		return false;
	}

	//////////////////////////
	/// Lepton pt
	float minpt_ldg = 25.;
	float minpt_subldg = 10;
	if (local.leptons_fakeable[1].Type() == LeptonType::kele)
		minpt_subldg = 15.;

	bool passLeptonPt =
		local.leptons_fakeable[0].pt() > minpt_ldg and
		local.leptons_fakeable[1].pt() > minpt_subldg;

	if (not passLeptonPt) {
		if (debug) {
			std::cout << "FAIL lepton pT cut" << std::endl;
			std::cout << "conePt : " << local.leptons_fakeable[0].conePt();
			std::cout << " " << local.leptons_fakeable[1].conePt() << std::endl;
		}
		return false;
	}

	//////////////////////////
	/// veto two loose leptons with invariant mass < 12 GeV
	bool passPairMassVeto = true;
	
	for (auto it = local.leptons_loose.begin(); it != local.leptons_loose.end()-1; ++it) {
		for (auto it2 = it+1; it2 != local.leptons_loose.end(); ++it2) {
			if ((it->p4() + it2->p4()).M() < 12.) {
				passPairMassVeto = false;
				break;
			}
		}
		if (not passPairMassVeto) break;
	}
	if (not passPairMassVeto) {
		if (debug) {
			std::cout << "FAIL any pair of loose leptons has invariant mass >= 12 GeV" << std::endl;
		}
		return false;
	}
	
	//////////////////////////
	/// Tight charge
	bool passTightCharge =
		local.leptons_fakeable[0].tightCharge() and
		local.leptons_fakeable[1].tightCharge();
	
	if (not passTightCharge) {
		if (debug) {
			std::cout << "FAIL tight charge" << std::endl;
			std::cout << "tight charge 1 : " << local.leptons_fakeable[0].tightCharge() << std::endl;
			if (local.leptons_fakeable[0].Type() == LeptonType::kele) {
				std::cout << "e" << std::endl;
				std::cout << "isGsfCtfScPixChargeConsistent : " << local.e_fakeable[0].isGsfCtfScPixChargeConsistent() << std::endl;
				std::cout << "isGsfScPixChargeConsistent : " << local.e_fakeable[0].isGsfScPixChargeConsistent() << std::endl;
			}
			
			std::cout << "tight charge 2 : " << local.leptons_fakeable[1].tightCharge() << std::endl;
			if (local.leptons_fakeable[1].Type() == LeptonType::kele) {
				std::cout << "e" << std::endl;
				std::cout << "isGsfCtfScPixChargeConsistent : " << local.e_fakeable[0].isGsfCtfScPixChargeConsistent() << std::endl;
				std::cout << "isGsfScPixChargeConsistent : " << local.e_fakeable[0].isGsfScPixChargeConsistent() << std::endl;
			}			
		}
		return false;
	}
	
	//////////////////////////
	/// Categorize
	bool passMetLD = false;
	bool passZmassVeto = false;
	
	if (local.leptons_fakeable[0].Type() == LeptonType::kmu and
		local.leptons_fakeable[1].Type() == LeptonType::kmu) {  // mumu

		ilep = 0;
		
		passMetLD = true;
		passZmassVeto = true;
	}
	else if (local.leptons_fakeable[0].Type() == LeptonType::kele and
			 local.leptons_fakeable[1].Type() == LeptonType::kele) {  // ee

		ilep = 1;

		//////////////////////////
		/// MetLD cut (ee only)
		passMetLD = local.metLD > 0.2;
		
		//////////////////////////
		/// Zmass Veto: 91.2 +/- 10 (ee only)
		double eeInvMass =
			(local.leptons_fakeable[0].p4() +
			 local.leptons_fakeable[1].p4()).M();

		passZmassVeto = eeInvMass < (91.2 - 10.0) or eeInvMass > (91.2 + 10.0);
	}
	else {  // emu
		ilep = 2;

		passMetLD = true;
		passZmassVeto = true;

		// determine which one is electron
		int ie = local.leptons_fakeable[1].Type() == LeptonType::kele;
		assert(local.leptons_fakeable[ie].Type() == LeptonType::kele);
		assert(local.leptons_fakeable[!ie].Type() == LeptonType::kmu);
	}	

	if (not passMetLD) {
		if (debug) {
			std::cout << "FAIL metLD cut" << std::endl;
			std::cout << "metLD : " << local.metLD << std::endl;
		}
		return false;
	}
	
	if (not passZmassVeto) {
		if (debug) {
			std::cout << "FAIL Zmass cut " << std::endl;
		}
		return false;
	}
	
	//////////////////////////
	/// number of jets and btags
   	int njets = local.jets_selected.size();
	int nbtags_loose = local.jets_selected_btag_loose.size();
	int nbtags_medium = local.jets_selected_btag_medium.size();

	bool passNumJets = njets >= 3;
	bool passNumBtags = nbtags_loose >= 2 or nbtags_medium >= 1;
	
	if (not passNumJets) {
		if (debug) std::cout << "FAIL number of jets requirement" << std::endl;
		return false;
	}
	
	if (not passNumBtags) {
		if (debug) std::cout << "FAIL number of btags requirement" << std::endl;
		return false;
	}
		
	ibtag = nbtags_medium >=2 ? 1 : 0;
	
	//////////////////////////
	///
	//////////////////////////	
	/// Lepton charge
	// same sign
	bool passLepCharge = (local.leptons_fakeable[0].charge() *
					 local.leptons_fakeable[1].charge()) > 0;
	bool passTauCharge = true;
	// for signal region, opposite sign between tau and either lepton
	// = local.leptons_fakeable[0].charge() * local.tau_selected[0].charge() < 0;
	// To save computing time, additional requirement on tau charge applied after ntuple production for either signal region and control region

	if (selection_region == Control_2los1tau) {
		passLepCharge = not passLepCharge;  // two leptons are opposite signs

		// + the lepton that has the same sign as tau has to be a electron
		// (and the charge flip rate is only applied to this electron)
		int tauCharge = local.tau_selected[0].charge();
		if (local.leptons_fakeable[0].charge() == tauCharge)
			passTauCharge = local.leptons_fakeable[0].Type()==LeptonType::kele;
		else if (local.leptons_fakeable[1].charge() == tauCharge)
			passTauCharge = local.leptons_fakeable[1].Type()==LeptonType::kele;
		else
			passTauCharge = false;	
	}

	if (not (passLepCharge and passTauCharge)) {
		if (debug) {
			std::cout << "FAIL lepton charge requirement" << std::endl;
		}
		return false;
	}

	//////////////////////////
	/// Leptons WP
	bool passLepSel = false;
	// signal region: two leading leptons are tight
	passLepSel =
		local.leptons_fakeable[0].passTightSel() and
		local.leptons_fakeable[1].passTightSel();

	if (selection_region == Control_1lfakeable) {
		// at least one lepton fails tight selection
		passLepSel = not passLepSel;
	} 
	
	if (not passLepSel) {
		if (debug) {
			std::cout << "FAIL lepton WP requirements" << std::endl;
		}
		return false;
	}

	//////////////////////////
	/// MC truth matching
	if (not isdata) {

		assert(local.leptons_fakeable.size()>=2);
		bool matchGenLeps = true;
		int lep_cnt = 0;
		
		for (const auto & lep : local.leptons_fakeable) {
			if (lep.MCMatchType==0) {  // initial value is set to 0
				std::cerr << "WARNING!! MC match type is not set!!" << std::endl;
			}
			assert(lep.MCMatchType!=0);
			
			if (lep.Type() == LeptonType::kele) {
				if (debug) {
					std::cout << "ele mc match: " << lep.MCMatchType << std::endl;
				}
				// require prompt or from prompt tau decay
				if (not (lep.MCMatchType == 1 or lep.MCMatchType == 3))
					matchGenLeps = false;
			}
			if (lep.Type() == LeptonType::kmu) {
				if (debug) {
					std::cout << "mu mc match: " << lep.MCMatchType << std::endl;
				}
				// require prompt or from prompt tau decay
				if (not (lep.MCMatchType == 2 or lep.MCMatchType == 4))
					matchGenLeps = false;
			}

			// match only the two leading leptons
			if (++lep_cnt > 1) break;
		}

		// Tau genMatch not used in event selection
		// split MC sample into _gentau and _faketau based on tau genMatch
		
		bool matchGenTau = false; 
		// should have at least one tau at this stage
		// Only match the leading one if there are more
		assert(local.tau_selected.size()>=1);
		int mtype = local.tau_selected[0].userFloat("MCMatchType");
		if (debug)
			std::cout << "tau mc match: " << mtype << std::endl;
		if (mtype==1 or mtype==2 or mtype==3 or mtype==4 or mtype==5) {
			matchGenTau = true;
		}
		
		// set Gen Matching flag
		local.isGenMatched = matchGenLeps and matchGenTau;
		
		// tauID scale factor
		local.tauID_sf =
			sf_helper->Get_TauIDSF(local.tau_selected[0],matchGenTau);
		// ToDo: move this upper level together with other weights and SF
		
		// for signal region
		if (selection_region == Signal_2lss1tau) {
			bool passMCMatching = matchGenLeps;
			
			assert(local.leptons_tight.size()==2);			
			if (not passMCMatching) {
				if (debug) std::cout << "FAIL MC truth matching" << std::endl;
				return false;
			}
		}
		
	}
	
	//////////////////////////
	if (debug) std::cout << "PASSED event seletion!" << std::endl;
	
	return true;
}

bool CU_ttH_EDA::pass_event_sel_3l(CU_ttH_EDA_event_vars &local,
								   Selection_types selection_type)
{
	//////////////////////////
    /// Lepton number
	// at least three fakeable leptons
	if (!(local.leptons_fakeable.size() >= 3)) return false;
	// signal region: three leading leptons are tight
    bool passLepSel =
		local.leptons_fakeable[0].passTightSel() and
		local.leptons_fakeable[1].passTightSel() and
		local.leptons_fakeable[2].passTightSel();

	if (not passLepSel) return false;

	//////////////////////////
	/// Lepton pt
	float minpt_ldg = 20.;
	float minpt_subldg = 10.;
	
	bool passLeptonPt =
		local.leptons_tight[0].conePt() > minpt_ldg and
		local.leptons_tight[1].conePt() > minpt_subldg  and
		local.leptons_tight[2].conePt() > minpt_subldg;

	if (not passLeptonPt) return false;

	//////////////////////////
	/// veto two loose leptons with invariant mass < 12 GeV
	bool passPairMassVeto = true;
	
	for (auto it = local.leptons_loose.begin(); it != local.leptons_loose.end()-1; ++it) {
		for (auto it2 = it+1; it2 != local.leptons_loose.end(); ++it2) {
			if ((it->p4() + it2->p4()).M() < 12.) {
				passPairMassVeto = false;
				break;
			}
		}
		if (not passPairMassVeto) break;
	}
	if (not passPairMassVeto) return false;

	
	bool hasSFOSpair = false;
	//////////////////////////
	/// veto loose lepton pair with invariant mass closer than 10 GeV to Z mass
	bool passZmassVeto = true;
		
	for (auto it = local.leptons_loose.begin(); it != local.leptons_loose.end()-1; ++it) {
		for (auto it2 = it+1; it2 != local.leptons_loose.end(); ++it2) {
			// same flavor
			if (it->Type() != it2->Type()) continue;
			// opposite charge
			if ((it->charge()) * (it2->charge()) >= 0) continue;

			hasSFOSpair = true;		
			float invMass = (it->p4() + it2->p4()).M();
			
			if (invMass > (91.2-10.0) and invMass < (91.2+10.0)) {
				passZmassVeto = false;
				break;
			}
		}
		if (not passZmassVeto) break;
	}

	if (selection_type == Control_WZ)
		passZmassVeto = not passZmassVeto;

	if (not passZmassVeto) return false;

	
	/// number of jets and btags
	int njets = local.jets_selected.size();
	int nbtags_loose = local.jets_selected_btag_loose.size();
	int nbtags_medium = local.jets_selected_btag_medium.size();
	
	//////////////////////////
	/// metLD cut
	float metLDcut = hasSFOSpair ? 0.3 : 0.2;
	bool passMetLD = (njets >= 4) or (local.metLD >= metLDcut);

	if (not passMetLD) return false;

	//////////////////////////
	/// charge
	bool passCharge = false;
	
	int total_charge =
		local.leptons_tight[0].charge() +
		local.leptons_tight[1].charge() +
		local.leptons_tight[2].charge();

	if (abs(total_charge) == 1) passCharge = true;

	if (not passCharge) return false;

	//////////////////////////
	/// To suppress electron from photon conversion
	bool passPhotonVeto = true;
	for (auto & l : local.leptons_fakeable) {
		if (l.Type() == LeptonType::kele) {
			passPhotonVeto = passPhotonVeto and
				l.conversionVeto() and l.noMissingHits();
		}
	}

	if (not passPhotonVeto) return false;
	
	//////////////////////////
	/// number of jets and btags
	bool passNumJets = njets >= 2;
	bool passNumBtags = nbtags_loose >= 2 or nbtags_medium >= 1;
	if (selection_type == Control_WZ) {
		passNumBtags = not passNumBtags;  // inconsistent between AN and twiki
	}

	if (not passNumJets) return false;
	if (not passNumBtags) return false;

	return true;
}

int CU_ttH_EDA::partition2DBDT(double ttbar, double ttV)
/*
  ICHEP binning
 */
{
	if (ttbar > -1.0 and ttbar <= -0.2 and ttV > -1.0 and ttV <= 1.0)
		return 1;
	else if (ttbar > -0.2 and ttbar <= 0.1 and ttV > -1.0 and ttV <= 1.0)
		return 2;
	else if (ttbar > 0.1 and ttbar <= 0.4) {
		if (ttV > -1.0 and ttV <= 0.3)
			return 3;
		else if (ttV > 0.3 and ttV <= 1.0)
			return 4;
	}
	else if (ttbar > 0.4 and ttbar <= 1.0) {
		if (ttV > -1.0 and ttV <= 0.1)
			return 5;
		else if (ttV > 0.1 and ttV <= 0.4)
			return 6;
		else if (ttV > 0.4 and ttV <= 1.0)
			return 7;
	}

	return 0;
}

std::vector<pat::Jet> CU_ttH_EDA::GetCorrectedJets(const std::vector<pat::Jet>& input_jets, JetCorrectionUncertainty* jecUnc, const std::string JECType)
{
	float shift = 0.;

	if (JECType=="JESUp")
		shift = 1.;
	else if (JECType=="JESDown")
		shift = -1.;
	else
		return input_jets;

	std::vector<pat::Jet> output_jets = input_jets;
	
	for (auto & jet : output_jets) {
		jecUnc->setJetEta(jet.eta());
		jecUnc->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
		float unc = jecUnc->getUncertainty(true);
	    float jes = 1 + shift*unc;
		jet.scaleEnergy(jes);
	}

	return output_jets;
}

std::vector<pat::Tau> CU_ttH_EDA::GetCorrectedTaus(const std::vector<pat::Tau>& input_taus, float tauES_Unc, const std::string TESType)
{
	float shift = 0.;

	if (TESType=="tauESUp")
		shift = 1.;
	else if (TESType=="tauESDown")
		shift = -1.;
	else
		return input_taus;

	std::vector<pat::Tau> output_taus = input_taus;

	for (auto & tau : output_taus) {
		auto corrP4 = tau.p4();
		corrP4 *= 1 + shift*tauES_Unc;
		tau.setP4(corrP4);
	}

	return output_taus;
}


int CU_ttH_EDA::HiggsDaughterPdgId(const std::vector<reco::GenParticle>& genParticles)
{
	for (auto & p : genParticles) {
		if (p.pdgId() != 25) continue;

		int ndaugs = p.numberOfDaughters();
		if (ndaugs != 2) continue;

		const reco::Candidate *d1 = p.daughter(0);
		const reco::Candidate *d2 = p.daughter(1);

		if ( abs(d1->pdgId()) != abs(d2->pdgId()) ) continue;

		//assert(p.statusFlags().isLastCopy());
		
		return d1->pdgId();
	}

	return -9999;
}

bool CU_ttH_EDA::HiggsDecayFilter(const std::vector<reco::GenParticle>& genParticles, const TString& decayMode)
{
	int required_id = 0;
	if (decayMode == "ttH_htt")
		required_id = 15;
	else if (decayMode == "ttH_hww")
		required_id = 24;
	else if (decayMode == "ttH_hzz")
		required_id = 23;

	int daug_id = HiggsDaughterPdgId(genParticles);
	
	return abs(daug_id) == required_id;

}

// MC Matching type encoding: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_Matching
const reco::GenParticle* CU_ttH_EDA::getMatchedGenParticle(const pat::Muon& patMu, const std::vector<reco::GenParticle>& gen_particles)
{
	const reco::GenParticle* out = NULL;
	float dRmin = 666.;
	
	// loop over genParticle collections to find match
	for (auto& gen : gen_particles) {
		if (abs(gen.pdgId()) == 13) {
			
			auto genStatus = gen.statusFlags();

			if (not (genStatus.isPrompt() or
					 genStatus.isDirectPromptTauDecayProduct()))
				continue;

			float dR = reco::deltaR(gen.eta(),gen.phi(),patMu.eta(),patMu.phi());
			
			if (dR > 0.2) continue;
			if (gen.pt() < 8.) continue;

			if (dR > dRmin) continue;  // find the closest in dR

			dRmin = dR;
			out = &gen;
		}
	}

	return out;
}

const reco::GenParticle* CU_ttH_EDA::getMatchedGenParticle(const pat::Electron& patEle, const std::vector<reco::GenParticle>& gen_particles)
{
	const reco::GenParticle* out = NULL;
	float dRmin = 666.;
	
	// loop over genParticle collections to find match
	for (auto& gen : gen_particles) {
		if (abs(gen.pdgId()) == 11) {
			auto genStatus = gen.statusFlags();

			if (not (genStatus.isPrompt() or
					 genStatus.isDirectPromptTauDecayProduct()))
				continue;

			float dR = reco::deltaR(gen.eta(),gen.phi(),patEle.eta(),patEle.phi());
			if (dR > 0.2) continue;
			if (gen.pt() < 8.) continue;

			if (dR > dRmin) continue;  // find the closest in dR
			
			dRmin = dR;
			out = &gen;
		}
	}

	return out;
}

const reco::GenParticle* CU_ttH_EDA::getMatchedGenParticle(const pat::Tau& patTau, const std::vector<reco::GenParticle>& gen_particles)
{
	const reco::GenParticle* out = NULL;
	float dRmin = 666.;
	
	// loop over genParticle collections to find match
	for (auto& gen : gen_particles) {
		if ( abs(gen.pdgId()) == 11 or abs(gen.pdgId()) == 13 ) {
			auto genStatus = gen.statusFlags();

			if (not (genStatus.isPrompt() or
					 genStatus.isDirectPromptTauDecayProduct()))
				continue;

			float dR = reco::deltaR(gen.eta(),gen.phi(),patTau.eta(),patTau.phi());

			if (dR > 0.2) continue;
			if (gen.pt() < 8.) continue;

		    if (dR > dRmin) continue;

			dRmin = dR;
			out = &gen;
		}

		if ( abs(gen.pdgId()) == 15 ) {
			auto genStatus = gen.statusFlags();

			if (not genStatus.isPrompt()) continue;

			reco::Candidate::LorentzVector visP4;
			for (unsigned int i = 0; i < gen.numberOfDaughters(); ++i) {
				auto daug = gen.daughter(i);
				int id = abs(daug->pdgId());
				if (id == 11 or id == 12 or id == 13 or id == 14 or id == 16)
					continue;
				visP4 += daug->p4();
			}

			float dR = reco::deltaR(visP4.eta(),visP4.phi(),patTau.eta(),patTau.phi());

			if (dR > 0.2) continue;
			// if (visP4.pt() < 15.) continue;
			if (gen.pt() < 15.) continue;

			if (dR > dRmin) continue;
			
			dRmin = dR;
			out = &gen;
		}
	}

	return out;
}


template <typename T>
int CU_ttH_EDA::MatchGenParticle_Type(const T& reco_particle, const std::vector<reco::GenParticle>& gen_particles)
{
	auto matchedGen =
		getMatchedGenParticle(reco_particle, gen_particles);
		//reco_particle.genParticle();

	if (matchedGen == NULL) return 6;

	if (debug) {
		std::cout << "gen pt eta phi pdgid: "<< matchedGen->pt() << " " << matchedGen->eta()<< " "<< matchedGen->phi()<<" "<< matchedGen->pdgId() << std::endl;
	}
	
	auto genStatus = matchedGen->statusFlags();

	int mtype = 6;
	
	if (abs(matchedGen->pdgId()) == 11) {
		if (genStatus.isPrompt()) mtype = 1;
		if (genStatus.isDirectPromptTauDecayProduct()) mtype = 3;
	}

	if (abs(matchedGen->pdgId()) == 13) {
		if (genStatus.isPrompt()) mtype = 2;
		if (genStatus.isDirectPromptTauDecayProduct())	mtype = 4;
	}

	if ( abs(matchedGen->pdgId()) == 15 and genStatus.isPrompt() ) {
		mtype = 5;
	}
	
	return mtype;
}

template int CU_ttH_EDA::MatchGenParticle_Type<pat::Electron>(const pat::Electron&, const std::vector<reco::GenParticle>&);
template int CU_ttH_EDA::MatchGenParticle_Type<pat::Muon>(const pat::Muon&, const std::vector<reco::GenParticle>&);
template int CU_ttH_EDA::MatchGenParticle_Type<pat::Tau>(const pat::Tau&, const std::vector<reco::GenParticle>&);

#endif
