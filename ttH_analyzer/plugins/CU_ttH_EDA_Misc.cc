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

/*
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
*/

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

	for (auto & mu : local.mu_preselected_sorted) {
		MHT_x -= mu.px();
		MHT_y -= mu.py();
	}

	for (auto & ele : local.e_preselected_sorted) {
		MHT_x -= ele.px();
		MHT_y -= ele.py();
	}

	for (auto & tau : local.tau_selected_sorted) {
		MHT_x -= tau.px();
		MHT_y -= tau.py();
	}

	for (auto & jet : local.jets_selected_sorted) {
		MHT_x -= jet.px();
		MHT_y -= jet.py();
	}

	return sqrt(MHT_x * MHT_x + MHT_y * MHT_y);
}

/*
bool CU_ttH_EDA::passExtraforTight(pat::Muon mu)
{
	if (mu.innerTrack().isAvailable()) {
		if (mu.innerTrack()->ptError()/mu.innerTrack()->pt() < 0.2)
			return true;
	}

	return false;
}

bool CU_ttH_EDA::passExtraforTight(pat::Electron ele)
{
	return (
			ele.userFloat("numMissingHits") == 0 and
			ele.passConversionVeto() and
			ele.isGsfCtfScPixChargeConsistent()
			);
}
*/

bool CU_ttH_EDA::pass_event_sel_2lss1tauh(CU_ttH_EDA_event_vars &local,
										  int jecType,
										  const std::string & selection_region)
{
	//////////////////////////
    /// Lepton number
	if (local.leptons_selected_sorted.size() != 2) return false;
	
	bool passLepSel = false;	 
	// exactly two tight leptons
	passLepSel =
		local.leptons_selected_sorted[0].passTightSel() and
		local.leptons_selected_sorted[1].passTightSel();
	
	if (selection_region == "control_1lfakeable") {
		// exactly one lepton fails tight selection
		passLepSel =
			(not local.leptons_selected_sorted[0].passTightSel() and
			 local.leptons_selected_sorted[1].passTightSel() )
			or (local.leptons_selected_sorted[0].passTightSel() and
				not local.leptons_selected_sorted[1].passTightSel());
	} 
	
	if (not passLepSel) return false;
	
	//////////////////////////
	/// veto two loose leptons with invariant mass < 12 GeV
	bool passPairMassVeto = true;
	std::vector<math::XYZTLorentzVector> leptons_p4;
	for (auto & mu : local.mu_preselected) {
		leptons_p4.push_back(mu.p4());
	}
	for (auto & ele : local.e_preselected) {
		leptons_p4.push_back(ele.p4());
	}

	for (auto it = leptons_p4.begin(); it != leptons_p4.end()-1; ++it) {
		for (auto it2 = it+1; it2 != leptons_p4.end(); ++it2) {
			if ( (*it + *it2).mass() <  12. )
				passPairMassVeto = false;
		}
	}

	if (not passPairMassVeto) return false;

	//////////////////////////
	/// Lepton charge
	bool passLepCharge = false;
	// same sign
	if (local.leptons_selected_sorted[0].charge() *
		local.leptons_selected_sorted[1].charge() > 0)
		passLepCharge = true;

	if (selection_region == "control_2los")
		passLepCharge = not passLepCharge;

	if (not passLepCharge) return false;

	//////////////////////////
	/// Lepton pt
	float minpt_ldg = 20.;
	float minpt_subldg = 10;
	if (local.leptons_selected_sorted[1].Type() == LeptonType::kele)
		minpt_subldg = 15.;

	bool passLeptonPt =
		local.leptons_selected_sorted[0].conePt() > minpt_ldg and
		local.leptons_selected_sorted[1].conePt() > minpt_subldg;

	if (not passLeptonPt) return false;


	//////////////////////////
	/// ee only cuts	
	bool passMetLD = true;
	bool passZmassVeto = true;
	if (local.leptons_selected_sorted[0].Type() == LeptonType::kele and
		local.leptons_selected_sorted[1].Type() == LeptonType::kele) {
		//////////////////////////
		/// MetLD cut (ee only)
		passMetLD = local.metLD > 0.2;
		
		//////////////////////////
		/// Zmass Veto: 91.2 +/- 10
		double eeInvMass =
			(local.leptons_selected_sorted[0].p4() +
			 local.leptons_selected_sorted[1].p4()).M();
		passZmassVeto = eeInvMass < (91.2 - 10.0) or eeInvMass > (91.2 + 10.0);
	}

	if (not passMetLD) return false;
	if (not passZmassVeto) return false;

	//////////////////////////
	/// number of taus
	bool passNumTaus = local.n_taus >= 1;
	if (not passNumTaus) return false;

	//////////////////////////
	/// number of jets and btags
   	int njets = 0;
	int nbtags_loose = 0;
	int nbtags_medium = 0;

	if (jecType == 1) {       // JESUp
		njets = local.jets_selected_jesup.size();
		nbtags_loose = local.jets_selected_btag_loose_jesup.size();
		nbtags_medium = local.jets_selected_btag_medium_jesup.size();
	}
	else if (jecType == -1) { // JESDown
		njets = local.jets_selected_jesdown.size();
		nbtags_loose = local.jets_selected_btag_loose_jesdown.size();
		nbtags_medium = local.jets_selected_btag_medium_jesdown.size();
	}
	else {                    // NA
		njets = local.jets_selected.size();
		nbtags_loose = local.jets_selected_btag_loose.size();
		nbtags_medium = local.jets_selected_btag_medium.size();
	}
	
	bool passNumJets = njets >= 4;
	bool passNumBtags = nbtags_loose >= 2 or nbtags_medium >= 1;
	
	if (selection_region == "control_2los") {
		passNumJets = njets >= 2;
		passNumBtags = nbtags_medium >= 1;
	}
	
	if (not passNumJets) return false;	
	if (not passNumBtags) return false;


	return true;
}

bool CU_ttH_EDA::pass_event_sel_1l2tauh(CU_ttH_EDA_event_vars &local, int sys, const std::string& selection_region)
{
	return false;
}

int CU_ttH_EDA::partition2DBDT(double ttbar, double ttV)
/*
             bin 1       bin 2       bin 3       bin 4      bin 5       bin 6
2lss(ttbar) (-1.0,-0.2] (-1.0,-0.2] (-0.2,0.3]  (-0.2,0.3] (0.3,1.0]   (0.3,1.0]
2lss(ttV)   (-1.0,-0.1] (-0.1,1.0]  (-1.0,-0.1] (-0.1,1.0] (-1.0,-0.1] (-0.1,1.0]
 */
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

	if (!(pt > 20. and fabs(eta) < 2.4)) return 1.;

	if (pt > 1000) pt = 999.;
	if (csv < 0.) csv = -0.05;
	if (csv > 1.) csv = 1.0;

	BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;	
	if ( abs(flavor) == 5 )
		jf = BTagEntry::FLAV_B;
	else if ( abs(flavor) == 4 )
		jf = BTagEntry::FLAV_C;

	double weight_jet = 1.;
	
	if (sys == "LFUp" or sys == "LFDown" or 
	    sys == "HFStats1Up" or sys == "HFStats1Down" or
	    sys == "HFStats2Up" or sys == "HFStats2Down") {
	  
	  if ( jf != BTagEntry::FLAV_B) sys = "NA";
	}

	if (sys == "HFUp" or sys == "HFDown" or
	    sys == "LFStats1Up" or sys == "LFStats1Down" or
	    sys == "LFStats2Up" or sys == "LFStats2Down") {

	  if (jf != BTagEntry::FLAV_UDSG) sys = "NA";
	}

	if (sys == "cErr1Up" or sys == "cErr1Down" or
	    sys == "cErr2Up" or sys == "cErr2Down") {

	  if (jf != BTagEntry::FLAV_C) sys = "NA";
	}

	weight_jet = BTagCaliReaders[sys]->eval(jf, eta, pt, csv);

	//assert(weight_jet > 0.);
	// problems in CSVv2 file? negative weight e.g. line 1735 and 1736
	// use csv_rwt_fit_hf/lf_76x_2016_02_08.root instead?
	// for now:
	if (weight_jet <= 0.) {
		weight_jet = 1.;
	}

	return weight_jet;
}

float CU_ttH_EDA::getEleChargeMisIDProb(const miniLepton& lepton, bool isdata)
{
	// muon
	if (lepton.Type() == LeptonType::kmu) return 0.;

	// electron
	assert(lepton.Type() == LeptonType::kele);

	if (isdata) {
		if (abs(lepton.eta()) < 1.479) {
			if (lepton.pt() >= 10 and lepton.pt() < 25)
				return 0.000301;
			else if (lepton.pt() < 50)
				return 0.000287;
			else
				return 0.000293;
		}
		else if (abs(lepton.eta()) < 2.5) {
			if (lepton.pt() >= 10 and lepton.pt() < 25)
				return 0.001728;
			else if (lepton.pt() < 50)
				return 0.001974;
			else
				return 0.003457;
		}

		return 0.;
	}
	else {
		if (abs(lepton.eta()) < 1.479) {
			if (lepton.pt() >= 10 and lepton.pt() < 25)
				return 0.000131;
			else if (lepton.pt() < 50)
				return 0.000255;
			else
				return 0.000340;
		}
		else if (abs(lepton.eta()) < 2.5) {
			if (lepton.pt() >= 10 and lepton.pt() < 25)
				return 0.000966;
			else if (lepton.pt() < 50)
				return 0.002160;
			else
				return 0.004170;
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

float CU_ttH_EDA::getFakeRate(const miniLepton& lepton)
{
	float fakerate = 0;
	
	if (lepton.Type() == LeptonType::kele)
		fakerate = read2DHist(h_fakerate_el, lepton.conePt(), lepton.eta());
	else if (lepton.Type() == LeptonType::kmu)
		fakerate = read2DHist(h_fakerate_mu, lepton.conePt(), lepton.eta());

	return fakerate;
}

void CU_ttH_EDA::Delete_BTagCalibration_Readers()
{
	delete BTagCaliReaders["NA"];
	delete BTagCaliReaders["JESUp"];
	delete BTagCaliReaders["JESDown"];
	delete BTagCaliReaders["LFUp"];
	delete BTagCaliReaders["LFDown"];
	delete BTagCaliReaders["HFUp"];
	delete BTagCaliReaders["HFDown"];
	delete BTagCaliReaders["HFStats1Up"];
	delete BTagCaliReaders["HFStats1Down"];
	delete BTagCaliReaders["HFStats2Up"];
	delete BTagCaliReaders["HFStats2Down"];
	delete BTagCaliReaders["LFStats1Up"];
	delete BTagCaliReaders["LFStats1Down"];
	delete BTagCaliReaders["LFStats2Up"];
	delete BTagCaliReaders["LFStats2Down"];
	delete BTagCaliReaders["cErr1Up"];
	delete BTagCaliReaders["cErr1Down"];
	delete BTagCaliReaders["cErr2Up"];
	delete BTagCaliReaders["cErr2Down"];
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
