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
			std::cout << "FAIL lepton number requirements" << std::endl;
		}
		return false;
	}
	
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
			std::cout << "FAIL lepton number requirements" << std::endl;
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
		local.leptons_fakeable[0].conePt() > minpt_ldg and
		local.leptons_fakeable[1].conePt() > minpt_subldg;

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
	/// Lepton charge
	bool passLepCharge = false;
	// same sign
	if (local.leptons_fakeable[0].charge() *
		local.leptons_fakeable[1].charge() > 0)
		passLepCharge = true;

	if (selection_region == Control_2los1tau)
		passLepCharge = not passLepCharge;

	if (not passLepCharge) {
		if (debug) {
			std::cout << "FAIL lepton charge requirement" << std::endl;
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
	bool passPhotonVeto = false;
	
	if (local.leptons_fakeable[0].Type() == LeptonType::kmu and
		local.leptons_fakeable[1].Type() == LeptonType::kmu) {  // mumu

		ilep = 0;
		
		passMetLD = true;
		passZmassVeto = true;
		passPhotonVeto = true;
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

		//////////////////////////
		/// To suppress electrons from photon conversions
		passPhotonVeto = local.leptons_fakeable[0].conversionVeto() and
			local.leptons_fakeable[0].noMissingHits() and
			local.leptons_fakeable[1].conversionVeto() and
			local.leptons_fakeable[1].noMissingHits();
	}
	else {  // emu
		ilep = 2;

		passMetLD = true;
		passZmassVeto = true;

		// determine which one is electron
		int ie = local.leptons_fakeable[1].Type() == LeptonType::kele;
		assert(local.leptons_fakeable[ie].Type() == LeptonType::kele);
		assert(local.leptons_fakeable[!ie].Type() == LeptonType::kmu);

		passPhotonVeto = local.leptons_fakeable[ie].conversionVeto() and
			local.leptons_fakeable[ie].noMissingHits();
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
	
	if (not passPhotonVeto) {
		if (debug) {
			std::cout << "FAIL photon conversion veto" << std::endl;
		}
		return false;
	}

	//////////////////////////
	/// number of taus
	bool passNumTaus = local.n_taus >= 1;
	if (not passNumTaus) {
		if (debug)
			std::cout << "FAIL number of taus requirement" << std::endl;
		return false;
	}
	
	//////////////////////////
	/// number of jets and btags
   	int njets = local.jets_selected.size();
	int nbtags_loose = local.jets_selected_btag_loose.size();
	int nbtags_medium = local.jets_selected_btag_medium.size();

	bool passNumJets = njets >= 4;
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

	if (debug) std::cout << "PASSED event seletion!" << std::endl;
	
	return true;
}

bool CU_ttH_EDA::pass_event_sel_3l(CU_ttH_EDA_event_vars &local,
								   Selection_types selection_region)
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

	if (selection_region == Control_WZ)
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
	if (selection_region == Control_WZ) {
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

	
/*
             bin 1       bin 2       bin 3       bin 4      bin 5       bin 6
2lss(ttbar) (-1.0,-0.2] (-1.0,-0.2] (-0.2,0.3]  (-0.2,0.3] (0.3,1.0]   (0.3,1.0]
2lss(ttV)   (-1.0,-0.1] (-0.1,1.0]  (-1.0,-0.1] (-0.1,1.0] (-1.0,-0.1] (-0.1,1.0]
 
{
	int x = -99;
	int y = -99;
	
	if (ttbar > -1.0 and ttbar <= -0.2)
		x = 1;
	else if (ttbar > -0.2 and ttbar <= 0.3)
		x = 2;
	else if (ttbar > 0.3 and ttbar <= 1.0)
		x = 3;

	if (ttV > -1.0  and ttV <= -0.1)
		y = 1;
	else if (ttV > -0.1 and ttV <= 1.0)
		y = 2;
	
	return 2*(x-1)+y;
}
*/

double CU_ttH_EDA::getEvtCSVWeight(std::vector<pat::Jet> & jets, int iSys)//, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF)
{
	int iSysHF = 0;
	switch(iSys){
	case 7:  iSysHF=1; break; //JESUp
	case 8:  iSysHF=2; break; //JESDown
	case 9:  iSysHF=3; break; //LFUp
	case 10: iSysHF=4; break; //LFDown
	case 13: iSysHF=5; break; //Stats1Up
	case 14: iSysHF=6; break; //Stats1Down
	case 15: iSysHF=7; break; //Stats2Up
	case 16: iSysHF=8; break; //Stats2Down
	default : iSysHF = 0; break; //NoSys
	}

	int iSysC = 0;
	switch(iSys){
	case 21: iSysC=1; break;
	case 22: iSysC=2; break;
	case 23: iSysC=3; break;
	case 24: iSysC=4; break;
	default : iSysC = 0; break;
	}

	int iSysLF = 0;
	switch(iSys){
	case 7:  iSysLF=1; break; //JESUp
	case 8:  iSysLF=2; break; //JESDown
	case 11: iSysLF=3; break; //HFUp
	case 12: iSysLF=4; break; //HFDown
	case 17: iSysLF=5; break; //Stats1Up
	case 18: iSysLF=6; break; //Stats1Down
	case 19: iSysLF=7; break; //Stats2Up
	case 20: iSysLF=8; break; //Stats2Down
	default : iSysLF = 0; break; //NoSys
	}

	double csvWgthf = 1.;
	double csvWgtC  = 1.;
	double csvWgtlf = 1.;

	for (auto & jet : jets) {
		double csv = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		double jetPt = jet.pt();
		double jetAbsEta = std::abs(jet.eta());
		int flavor = jet.hadronFlavour();

		int iPt = -1; int iEta = -1;
		if (jetPt >=19.99 && jetPt<30) iPt = 0;
		else if (jetPt >=30 && jetPt<40) iPt = 1;
		else if (jetPt >=40 && jetPt<60) iPt = 2;
		else if (jetPt >=60 && jetPt<100) iPt = 3;
		else if (jetPt >=100) iPt = 4;

		if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
		else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
		else if ( jetAbsEta>=1.6 && jetAbsEta<2.41 ) iEta = 2;

		 if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;

		 if (abs(flavor) == 5 ){
			 int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
			 double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
			 if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;
			 assert(csvWgthf > 0.);
		 }
		 else if( abs(flavor) == 4 ){
			 int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
			 double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
			 if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;
			 assert(csvWgtC > 0.);
		 }
		 else {
			 if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
			 int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
			 double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
			 if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;
			 assert(csvWgtlf);
		 }
		 
	}  // end of jet loop

	double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

	//csvWgtHF = csvWgthf;
	//csvWgtLF = csvWgtlf;
	//csvWgtCF = csvWgtC;
	
	return csvWgtTotal;
}

/*
double CU_ttH_EDA::getEvtCSVWeight(std::vector<pat::Jet> & jets, std::string & sys)
{
	double weight_evt = 1.;

	for (auto & j : jets) {		
		double w = getJetCSVWeight(j, sys);
		weight_evt *= w;
	}

	return weight_evt;
}

double CU_ttH_EDA::getJetCSVWeight(pat::Jet & jet, std::string sys)
{
	double pt = jet.pt();
	double eta = jet.eta();
	double csv = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	int flavor = jet.hadronFlavour();

	BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;	
	if ( abs(flavor) == 5 )
		jf = BTagEntry::FLAV_B;
	else if ( abs(flavor) == 4 )
		jf = BTagEntry::FLAV_C;

	double weight_jet = 1.;

	if (sys == "JESUp")
		sys = "up_jes";
	else if (sys == "JESDown")
		sys = "down_jes";
	else if (sys == "LFUp" and jf == BTagEntry::FLAV_B)
		sys = "up_lf";
	else if (sys == "LFDown" and jf == BTagEntry::FLAV_B)
		sys = "down_lf";
	else if (sys == "HFStats1Up" and jf == BTagEntry::FLAV_B)
		sys = "up_hfstats1";
	else if (sys == "HFStats1Down" and jf == BTagEntry::FLAV_B)
		sys = "down_hfstats1";
	else if (sys == "HFStats2Up" and jf == BTagEntry::FLAV_B)
		sys = "up_hfstats2";
	else if (sys == "HFStats2Down" and jf == BTagEntry::FLAV_B)
		sys = "down_hfstats2";
	else if (sys == "HFUp" and jf == BTagEntry::FLAV_UDSG)
		sys = "up_hf";
	else if (sys == "HFDown" and jf == BTagEntry::FLAV_UDSG)
		sys = "down_hf";
	else if (sys == "LFStats1Up" and jf == BTagEntry::FLAV_UDSG)
		sys = "up_lfstats1";
	else if (sys == "LFStats1Down" and jf == BTagEntry::FLAV_UDSG)
		sys = "down_lfstats1";
	else if (sys == "LFStats2Up" and jf == BTagEntry::FLAV_UDSG)
		sys = "up_lfstats2";
	else if (sys == "LFStats2Down" and jf == BTagEntry::FLAV_UDSG)
		sys = "down_lfstats2";
	else if (sys == "cErr1Up" and jf == BTagEntry::FLAV_C)
		sys = "up_cferr1";
	else if (sys == "cErr1Down" and jf == BTagEntry::FLAV_C)
		sys = "down_cferr1";
	else if (sys == "cErr2Up" and jf == BTagEntry::FLAV_C)
		sys = "up_cferr2";
	else if (sys == "cErr2Down" and jf == BTagEntry::FLAV_C)
		sys = "down_cferr2";
	else
		sys = "central";

	weight_jet = BTagCaliReader->eval_auto_bounds(sys, jf, eta, pt, csv);

	assert(weight_jet > 0.);
	// problems in CSVv2 file? negative weight e.g. line 1735 and 1736
	// for now:
	//if (weight_jet <= 0.) {
	//	weight_jet = 1.;
	//}

	return weight_jet;
}
*/

float CU_ttH_EDA::getEleChargeMisIDProb(const miniLepton& lepton, bool isdata)
{
	// muon
	if (lepton.Type() == LeptonType::kmu) return 0.;

	// electron
	assert(lepton.Type() == LeptonType::kele);

	if (isdata) {
		if (abs(lepton.eta()) < 1.479) {
			if (lepton.pt() >= 10 and lepton.pt() < 25)
				return 0.000337;
			else if (lepton.pt() < 50)
				return 0.000259;
			else
				return 0.000403;
		}
		else if (abs(lepton.eta()) < 2.5) {
			if (lepton.pt() >= 10 and lepton.pt() < 25)
				return 0.001476;
			else if (lepton.pt() < 50)
				return 0.002599;
			else
				return 0.003963;
		}

		return 0.;
	}
	else {
		if (abs(lepton.eta()) < 1.479) {
			if (lepton.pt() >= 10 and lepton.pt() < 25)
				return 0.000608;
			else if (lepton.pt() < 50)
				return 0.000296;
			else
				return 0.000177;
		}
		else if (abs(lepton.eta()) < 2.5) {
			if (lepton.pt() >= 10 and lepton.pt() < 25)
				return 0.001047;
			else if (lepton.pt() < 50)
				return 0.002376;
			else
				return 0.003101;
		}

		return 0.;
	}

	return 0.;
	
}

float CU_ttH_EDA::read2DHist(TH2* h2d, float x, float y)
{
	TAxis* xaxis = h2d->GetXaxis();
	int nbinx = xaxis->GetNbins();
	int xbin = xaxis->FindBin(x);
    if (xbin < 1) xbin = 1;
	if (xbin > nbinx) xbin = nbinx;

	TAxis* yaxis = h2d->GetYaxis();
	int nbiny = yaxis->GetNbins();
	int ybin = yaxis->FindBin(abs(y));
    if (ybin < 1) ybin = 1;
	if (ybin > nbiny) ybin = nbiny;

	float result = h2d->GetBinContent(xbin, ybin);

	return result;
}

float CU_ttH_EDA::readTGraph(TGraphAsymmErrors* graph, float x)
{
	float x1 = std::max(float(graph->GetXaxis()->GetXmin()+1e-5),
						std::min(float(graph->GetXaxis()->GetXmax()-1e-5), x)
						);
	return graph->Eval(x1);
}

float CU_ttH_EDA::getFakeRate(const miniLepton& lepton)
{
	float fakerate = 0;
	
	if (lepton.Type() == LeptonType::kele)
		fakerate = read2DHist(h_fakerate_el, lepton.conePt(), lepton.eta());
	else if (lepton.Type() == LeptonType::kmu)
		fakerate = read2DHist(h_fakerate_mu, lepton.conePt(), lepton.eta());

	return fakerate;
}

float CU_ttH_EDA::getLeptonSF(const miniLepton& lepton)
{
	float sf = 1.;
	assert(lepton.passLooseSel());
	
	sf *= getLeptonSF_loose(lepton);
	// SF fakeable to loose is currently assumed to be 1.
	if (lepton.passTightSel())
		sf *= getLeptonSF_tight_vs_loose(lepton);

	return sf;
}

float CU_ttH_EDA::getLeptonSF_loose(const miniLepton& lepton)
{
	// TODO: implement systematic uncertainty
	
	float sf = 1.;
	
	if (lepton.Type()==LeptonType::kmu) {
		if (abs(lepton.eta())<1.2) {
			sf *= readTGraph(h_recoToLoose_leptonSF_mu1_b, lepton.pt());
		}
		else {
			sf *= readTGraph(h_recoToLoose_leptonSF_mu1_e, lepton.pt());
		}

		sf *= read2DHist(h_recoToLoose_leptonSF_mu2,
						 lepton.pt(), abs(lepton.eta()));

		sf *= readTGraph(h_recoToLoose_leptonSF_mu3, lepton.eta());
	}
	else if (lepton.Type()==LeptonType::kele) {
		sf *= read2DHist(h_recoToLoose_leptonSF_el1,
						 lepton.pt(), abs(lepton.eta()));
		sf *= read2DHist(h_recoToLoose_leptonSF_el2,
						 lepton.pt(), abs(lepton.eta()));
		sf *= read2DHist(h_recoToLoose_leptonSF_el3,
						 lepton.pt(), abs(lepton.eta()));
		// !! different pt eta xaxis
		sf *= read2DHist(h_recoToLoose_leptonSF_gsf, lepton.eta(), lepton.pt());
	}

	return sf;
}

float CU_ttH_EDA::getLeptonSF_tight_vs_loose(const miniLepton& lepton)
{
	float sf = 1.;

	assert(lepton.passTightSel());
	
	// TODO: implement systematic uncertainty
	if (analysis_type == Analyze_2lss1tau) {
		if (lepton.Type()==LeptonType::kmu) {
			sf = read2DHist(h_looseToTight_leptonSF_mu_2lss,
							lepton.pt(), abs(lepton.eta()));
		}
		else if (lepton.Type()==LeptonType::kele) {
			sf = read2DHist(h_looseToTight_leptonSF_el_2lss,
							lepton.pt(), abs(lepton.eta()));
		}
	}
	else if (analysis_type == Analyze_3l) {
		if (lepton.Type()==LeptonType::kmu) {
			sf = read2DHist(h_looseToTight_leptonSF_mu_3l,
							lepton.pt(), abs(lepton.eta()));
		}
		else if (lepton.Type()==LeptonType::kele) {
			sf = read2DHist(h_looseToTight_leptonSF_el_3l,
							lepton.pt(), abs(lepton.eta()));
		}
	}

	return sf;
}

float CU_ttH_EDA::getLepHLTSF(int ilep)
{
	if (ilep == 0) // mumu
		return 1.01;
	else if (ilep == 1) // ee
		return 1.02;
	else if (ilep == 2) // emu
		return 1.02;
	
	std::cerr << "not valid lepton category !" << std::endl;
	assert(0);
	return 0.;
}

bool CU_ttH_EDA::HiggsDecayFilter(const std::vector<reco::GenParticle>& genParticles, const TString& decayMode)
{
	for (auto & p : genParticles) {
		if (p.pdgId() != 25) continue;

		int ndaugs = p.numberOfDaughters();
		if (ndaugs != 2) continue;

		const reco::Candidate *d1 = p.daughter(0);
		const reco::Candidate *d2 = p.daughter(1);

		if (abs(d1->pdgId()) != abs(d2->pdgId()) ) continue;

		int daug_id = abs(d1->pdgId());

		int required_id = 0;
		if (decayMode == "ttH_htt")
			required_id = 15;
		else if (decayMode == "ttH_hww")
			required_id = 24;
		else if (decayMode == "ttH_hzz")
			required_id = 23;

		return daug_id == required_id;
	}

	return false;
}

#endif
