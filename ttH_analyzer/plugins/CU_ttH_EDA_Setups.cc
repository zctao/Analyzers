#ifndef CU_ttH_EDA_Setups_cc
#define CU_ttH_EDA_Setups_cc

/// Includes
#include "CU_ttH_EDA.h"

void CU_ttH_EDA::Close_output_files()
{
	if (analysis_type == Analyze_lepton_jet) {
		fclose(events_e_cut1);
		fclose(events_e_cut2);
		fclose(events_e_cut3);
		fclose(events_e_cut4);
		fclose(events_e_cut5);
		fclose(events_e_cut6);
		fclose(events_e_cut7);

		fclose(events_mu_cut1);
		fclose(events_mu_cut2);
		fclose(events_mu_cut3);
		fclose(events_mu_cut4);
		fclose(events_mu_cut5);
		fclose(events_mu_cut6);
		fclose(events_mu_cut7);
	}

	if (analysis_type == Analyze_dilepton) {
		fclose(events_dimu_cut1);
		fclose(events_dimu_cut2);
		fclose(events_dimu_cut3);
		fclose(events_dimu_cut4);
		fclose(events_dimu_cut5);
		fclose(events_dimu_cut6);
		fclose(events_dimu_cut7);

		fclose(events_diele_cut1);
		fclose(events_diele_cut2);
		fclose(events_diele_cut3);
		fclose(events_diele_cut4);
		fclose(events_diele_cut5);
		fclose(events_diele_cut6);
		fclose(events_diele_cut7);

		fclose(events_elemu_cut1);
		fclose(events_elemu_cut2);
		fclose(events_elemu_cut3);
		fclose(events_elemu_cut4);
		fclose(events_elemu_cut5);
	}
}

void CU_ttH_EDA::Set_up_histograms()
{
	// 	h_electron_selection = fs_->make<TH1D>("h_electron_selection",
	// ";electron cut", 12, 0 , 12 );
	// 	h_muon_selection = fs_->make<TH1D>("h_muon_selection", ";muon cut", 12,
	// 0 , 12 );
	//
	// 	h_electron_selection->GetXaxis()->SetBinLabel(1, "All");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(2, "p_{T}>20, |#eta|<2.4,
	// !inCrack");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(3, "full5x5 #sigma_{i#eta
	// i#eta}");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(4, "|#Delta #eta_{in}|");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(5, "|#Delta #phi_{in}|");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(6, "hOverE");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(7, "1/E - 1/p");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(8, "d0");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(9, "dZ");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(10,
	// "expectedMissingInnerHits");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(11, "passConversionVeto");
	// 	h_electron_selection->GetXaxis()->SetBinLabel(12, "relIso (#Delta Beta,
	// 0.3)");
	//
	// 	h_muon_selection->GetXaxis()->SetBinLabel(1, "All");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(2, "p_{T}>20, |#eta|<2.4");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(3, "GlobalMuon ||
	// TrackerMuon");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(4, "PFMuon");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(5, "#Chi^{2}");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(6, "validMuonHit");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(7, "validPixelHit");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(8, "trk layers w/meas");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(9, "matched stations");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(10, "d0");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(11, "dZ");
	// 	h_muon_selection->GetXaxis()->SetBinLabel(12, "relIso < 0.1");

	
	if (analysis_type == Analyze_lepton_jet) {
		h_tth_syncex1_ele =
			fs_->make<TH1D>("h_tth_syncex1_ele", ";cut", 8, 0, 8);
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(2, "Single ele trig");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(3, "==1 electron");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(4, "==0 muons");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(5, ">=4 jets");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(6, ">=2 b-tags");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(7, ">=1 top-tags");
		h_tth_syncex1_ele->GetXaxis()->SetBinLabel(8, ">=1 Higgs-tags");

		h_tth_syncex1_mu = fs_->make<TH1D>("h_tth_syncex1_mu", ";cut", 8, 0, 8);
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(2, "Single mu trig");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(3, "==1 muon");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(4, "==0 electrons");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(5, ">=4 jets");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(6, ">=2 b-tags");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(7, ">=1 top-tags");
		h_tth_syncex1_mu->GetXaxis()->SetBinLabel(8, ">=1 Higgs-tags");
	}

	if (analysis_type == Analyze_dilepton) {
		h_tth_syncex1_dimu =
			fs_->make<TH1D>("h_tth_syncex1_dimu", ";cut", 8, 0, 8);
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(2, "Double mu trig");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(3, ">=2 muons");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(4, "Mll > 20");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(5, "Z Veto   ");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(6, ">=2 jets");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(7, "MET > 40");
		h_tth_syncex1_dimu->GetXaxis()->SetBinLabel(8, ">=1 b-tags");

		h_tth_syncex1_diele =
			fs_->make<TH1D>("h_tth_syncex1_diele", ";cut", 8, 0, 8);
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(2, "Double ele trig");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(3, ">=2 electrons");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(4, "Mll > 20");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(5, "Z Veto   ");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(6, ">=2 jets");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(7, "MET > 40");
		h_tth_syncex1_diele->GetXaxis()->SetBinLabel(8, ">=1 b-tags");

		h_tth_syncex1_elemu =
			fs_->make<TH1D>("h_tth_syncex1_elemu", ";cut", 6, 0, 6);
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(1, "All events");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(2, "Ele-mu trig");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(3, ">=2 leptons");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(4, "Mll > 20");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(5, ">=2 jets");
		h_tth_syncex1_elemu->GetXaxis()->SetBinLabel(6, ">=1 b-tags");
	}

	if (analysis_type == Analyze_taus_dilepton) {
		h_tth_syncex_dileptauh =
			fs_->make<TH1D>("h_tth_syncex_dileptauh", ";cut", 8, 0, 8);
		h_tth_syncex_dileptauh->GetXaxis()->SetBinLabel(1, "All events");
		//h_tth_syncex_dileptauh->GetXaxis()->SetBinLabel(2, "Weighted");
		h_tth_syncex_dileptauh->GetXaxis()->SetBinLabel(2, "HLT");
		h_tth_syncex_dileptauh->GetXaxis()->SetBinLabel(3, ">= 2 lep");
		h_tth_syncex_dileptauh->GetXaxis()->SetBinLabel(4, "1 tauh");
		h_tth_syncex_dileptauh->GetXaxis()->SetBinLabel(5, "Same sign lep");
		h_tth_syncex_dileptauh->GetXaxis()->SetBinLabel(6, "Opposite sign tau");
		h_tth_syncex_dileptauh->GetXaxis()->SetBinLabel(7, "min_njets");
		h_tth_syncex_dileptauh->GetXaxis()->SetBinLabel(8, "min_nbtags");

		// Jet multiplicity
		h_njets = fs_->make<TH1D>("h_njets", "", 11, -0.5, 10.5);
		h_nbtags = fs_->make<TH1D>("h_nbtags", "", 6, -0.5, 5.5);

		// Tau
		h_ntauID = fs_->make<TH1D>("h_ntauID", "", 4, -0.5, 3.5);
		h_ntauID->GetXaxis()->SetBinLabel(1, "nonIso");
		h_ntauID->GetXaxis()->SetBinLabel(2, "loose");
		h_ntauID->GetXaxis()->SetBinLabel(3, "medium");
		h_ntauID->GetXaxis()->SetBinLabel(4, "tight");
		// Tau selection efficiency histograms
		h_num_genHadTau = fs_->make<TH1D>("h_num_genHadTau","", 4, 0, 4);
		h_genHadTau_pt = fs_->make<TH1D>("h_genHadTau_pt", "", 80, 0, 400);
		h_genHadTau_eta = fs_->make<TH1D>("h_genHadTau_eta","", 50,-2.5,2.5);
		h_genHadTau_phi = fs_->make<TH1D>("h_genHadTau_phi","", 64, -3.2, 3.2);
		//h_num_FakeTau = fs_make<TH1D>("h_num_FakeTau","", 4, 0, 4);

		h_num_selectedTau_noniso =
			fs_->make<TH1D>("h_num_selectedTau_noniso","", 4, 0, 4);
		h_selectedTau_noniso_genpt =
			fs_->make<TH1D>("h_selectedTau_noniso_genpt", "", 80, 0, 400);
		h_selectedTau_noniso_geneta =
			fs_->make<TH1D>("h_selectedTau_noniso_geneta","", 50,-2.5,2.5);
		h_selectedTau_noniso_genphi =
			fs_->make<TH1D>("h_selectedTau_noniso_genphi","", 64, -3.2, 3.2);

		h_num_selectedTau_loose =
			fs_->make<TH1D>("h_num_selectedTau_loose","", 4, 0, 4);
		h_selectedTau_loose_genpt =
			fs_->make<TH1D>("h_selectedTau_loose_genpt", "", 80, 0, 400);
		h_selectedTau_loose_geneta =
			fs_->make<TH1D>("h_selectedTau_loose_geneta","", 50,-2.5,2.5);
		h_selectedTau_loose_genphi =
			fs_->make<TH1D>("h_selectedTau_loose_genphi","", 64, -3.2, 3.2);

		h_num_selectedTau_medium =
			fs_->make<TH1D>("h_num_selectedTau_medium","", 4, 0, 4);
		h_selectedTau_medium_genpt =
			fs_->make<TH1D>("h_selectedTau_medium_genpt", "", 80, 0, 400);
		h_selectedTau_medium_geneta =
			fs_->make<TH1D>("h_selectedTau_medium_geneta","", 50,-2.5,2.5);
		h_selectedTau_medium_genphi =
			fs_->make<TH1D>("h_selectedTau_medium_genphi","", 64, -3.2, 3.2);

		h_num_selectedTau_tight =
			fs_->make<TH1D>("h_num_selectedTau_tight","", 4, 0, 4);
		h_selectedTau_tight_genpt =
			fs_->make<TH1D>("h_selectedTau_tight_genpt", "", 80, 0, 400);
		h_selectedTau_tight_geneta =
			fs_->make<TH1D>("h_selectedTau_tight_geneta","", 50,-2.5,2.5);
		h_selectedTau_tight_genphi =
			fs_->make<TH1D>("h_selectedTau_tight_genphi","", 64, -3.2, 3.2);
	}
	
	if (analysis_type == Analyze_taus_lepton_jet) {
		auto nbins = cuts.size()+1;
		h_tth_syncex_eleditauh =
			fs_->make<TH1D>("h_tth_syncex_eleditauh", ";cut", nbins, 0, nbins);
		h_tth_syncex_eleditauh->GetXaxis()->SetBinLabel(1, "All events");
		//h_tth_syncex_eleditauh->GetXaxis()->SetBinLabel(2, "Weighted");
		int itr_e = 2;
		for (auto & icut : cuts) {
			h_tth_syncex_eleditauh->GetXaxis()->SetBinLabel(itr_e++,icut.c_str());
		}

		h_tth_syncex_muditauh =
			fs_->make<TH1D>("h_tth_syncex_muditauh", ";cut", nbins, 0, nbins);
		h_tth_syncex_muditauh->GetXaxis()->SetBinLabel(1, "All events");
		//h_tth_syncex_muditauh->GetXaxis()->SetBinLabel(2, "Weighted");
		int itr_mu = 2;
		for (auto & icut : cuts) {
			h_tth_syncex_muditauh->GetXaxis()->SetBinLabel(itr_mu++,icut.c_str());
		}
	}
}

/// Prepare histograms for trigger/filter counts
int CU_ttH_EDA::Set_up_Run_histograms_triggers()
{
	unsigned int numHLT = trigger_names_no_ver.size();
	h_hlt = fs_->make<TH1D>("h_hlt", ";HLT path", numHLT, 0, numHLT);
	if (!h_hlt)
		return 1;

	TAxis *axis = h_hlt->GetXaxis();
	if (!axis)
		return 1;

	for (unsigned int i = 0; i < numHLT; ++i)
		axis->SetBinLabel(i + 1, trigger_names_no_ver[i].c_str());

	unsigned int numFLT = filter_names_no_ver.size();
	h_flt = fs_->make<TH1D>("h_flt", ";Filter path", numFLT, 0, numFLT);
	if (!h_flt)
		return 1;

	axis = h_flt->GetXaxis();
	if (!axis)
		return 1;

	for (unsigned int i = 0; i < numFLT; ++i)
		axis->SetBinLabel(i + 1, filter_names_no_ver[i].c_str());

	return 0;
}

void CU_ttH_EDA::Set_up_trigger_name_vectors()
{
	/// Fill trigger name vectors and counters
	trigger_names = hlt_config.triggerNames();

	trigger_names_no_ver.clear();
	trigger_names_no_ver.push_back("All");
	std::string prefix = "HLT_";
	for (unsigned int i = 0; i < trigger_names.size(); ++i) {
		std::string pathNameNoVer = hlt_config.removeVersion(trigger_names[i]);

		n_trigger_fired[pathNameNoVer] = 0;

		if (trigger_names[i].compare(0, prefix.length(), prefix) == 0)
			trigger_names_no_ver.push_back(pathNameNoVer);
	}

	/// Fill filter name vectors and counters
	filter_names = filter_config.triggerNames();

	filter_names_no_ver.clear();
	filter_names_no_ver.push_back("All");
	for (unsigned int i = 0; i < filter_names.size(); ++i) {
		std::string pathNameNoVer =
			filter_config.removeVersion(filter_names[i]);

		n_filter_fired[pathNameNoVer] = 0;

		filter_names_no_ver.push_back(pathNameNoVer);
	}
}

/// Make and open write-out files
void CU_ttH_EDA::Set_up_output_files()
{
	if (analysis_type == Analyze_lepton_jet) {
		events_e_cut1 = fopen("Outputs/CU_events_e_cut1.dat", "w");
		events_e_cut2 = fopen("Outputs/CU_events_e_cut2.dat", "w");
		events_e_cut3 = fopen("Outputs/CU_events_e_cut3.dat", "w");
		events_e_cut4 = fopen("Outputs/CU_events_e_cut4.dat", "w");
		events_e_cut5 = fopen("Outputs/CU_events_e_cut5.dat", "w");
		events_e_cut6 = fopen("Outputs/CU_events_e_cut6.dat", "w");
		events_e_cut7 = fopen("Outputs/CU_events_e_cut7.dat", "w");

		events_mu_cut1 = fopen("Outputs/CU_events_mu_cut1.dat", "w");
		events_mu_cut2 = fopen("Outputs/CU_events_mu_cut2.dat", "w");
		events_mu_cut3 = fopen("Outputs/CU_events_mu_cut3.dat", "w");
		events_mu_cut4 = fopen("Outputs/CU_events_mu_cut4.dat", "w");
		events_mu_cut5 = fopen("Outputs/CU_events_mu_cut5.dat", "w");
		events_mu_cut6 = fopen("Outputs/CU_events_mu_cut6.dat", "w");
		events_mu_cut7 = fopen("Outputs/CU_events_mu_cut7.dat", "w");
	}

	if (analysis_type == Analyze_dilepton) {
		events_dimu_cut1 = fopen("Outputs/CU_events_dimu_cut1.dat", "w");
		events_dimu_cut2 = fopen("Outputs/CU_events_dimu_cut2.dat", "w");
		events_dimu_cut3 = fopen("Outputs/CU_events_dimu_cut3.dat", "w");
		events_dimu_cut4 = fopen("Outputs/CU_events_dimu_cut4.dat", "w");
		events_dimu_cut5 = fopen("Outputs/CU_events_dimu_cut5.dat", "w");
		events_dimu_cut6 = fopen("Outputs/CU_events_dimu_cut6.dat", "w");
		events_dimu_cut7 = fopen("Outputs/CU_events_dimu_cut7.dat", "w");

		events_diele_cut1 = fopen("Outputs/CU_events_diele_cut1.dat", "w");
		events_diele_cut2 = fopen("Outputs/CU_events_diele_cut2.dat", "w");
		events_diele_cut3 = fopen("Outputs/CU_events_diele_cut3.dat", "w");
		events_diele_cut4 = fopen("Outputs/CU_events_diele_cut4.dat", "w");
		events_diele_cut5 = fopen("Outputs/CU_events_diele_cut5.dat", "w");
		events_diele_cut6 = fopen("Outputs/CU_events_diele_cut6.dat", "w");
		events_diele_cut7 = fopen("Outputs/CU_events_diele_cut7.dat", "w");

		events_elemu_cut1 = fopen("Outputs/CU_events_elemu_cut1.dat", "w");
		events_elemu_cut2 = fopen("Outputs/CU_events_elemu_cut2.dat", "w");
		events_elemu_cut3 = fopen("Outputs/CU_events_elemu_cut3.dat", "w");
		events_elemu_cut4 = fopen("Outputs/CU_events_elemu_cut4.dat", "w");
		events_elemu_cut5 = fopen("Outputs/CU_events_elemu_cut5.dat", "w");
	}
}

void CU_ttH_EDA::Set_up_tokens()
{
	token.event_gen_info =
		consumes<GenEventInfoProduct>(edm::InputTag(std::string("generator")));
	token.triggerResults = consumes<edm::TriggerResults>(
		edm::InputTag(std::string("TriggerResults"), std::string(""), hltTag));
	token.filterResults = consumes<edm::TriggerResults>(edm::InputTag(
		std::string("TriggerResults"), std::string(""), filterTag));

	token.vertices = consumes<reco::VertexCollection>(
		edm::InputTag(std::string("offlineSlimmedPrimaryVertices")));
	token.sec_vertices = consumes<reco::VertexCompositePtrCandidateCollection>(
		edm::InputTag(std::string("slimmedSecondaryVertices")));
	token.PU_info = consumes<std::vector<PileupSummaryInfo>>(
		edm::InputTag(std::string("addPileupInfo")));
	token.srcRho = consumes<double>(
		edm::InputTag("fixedGridRhoFastjetAll"));

	token.electrons = consumes<pat::ElectronCollection>(
		edm::InputTag(std::string("slimmedElectrons")));
	token.muons = consumes<pat::MuonCollection>(
		edm::InputTag(std::string("slimmedMuons")));
	token.taus = consumes<pat::TauCollection>(
		edm::InputTag(std::string("slimmedTaus")));

	token.jets =
		consumes<pat::JetCollection>(edm::InputTag(std::string("slimmedJets")));
	token.METs =
		consumes<pat::METCollection>(edm::InputTag(std::string("slimmedMETs")));
	token.METSig =
		consumes<double>(edm::InputTag(std::string("METSignificance"),std::string("METSignificance")));
	token.METCov =
		consumes<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>>(edm::InputTag(std::string("METSignificance"),std::string("METCovariance")));
	
	token.MC_particles = consumes<reco::GenParticleCollection>(
		edm::InputTag(std::string("prunedGenParticles")));
	token.PF_candidates = consumes<pat::PackedCandidateCollection>(
		edm::InputTag(std::string("packedPFCandidates")));

	token.BS =
		consumes<reco::BeamSpot>(edm::InputTag(std::string("offlineBeamSpot")));

	token.top_jets = consumes<boosted::HTTTopJetCollection>(
		edm::InputTag("HTTTopJetsPFMatcher", "heptopjets", "p"));
	token.subfilter_jets = consumes<boosted::SubFilterJetCollection>(
		edm::InputTag("CA12JetsCA3FilterjetsPFMatcher", "subfilterjets", "p"));

	token.MC_packed = consumes<pat::PackedGenParticleCollection>(
		edm::InputTag(std::string("packedGenParticles")));
}

void CU_ttH_EDA::Setup_Tree()
{

	eventTree = fs_->make<TTree>("EventTree", "Event tree");

	eventTree->Branch("n_electrons", &n_electrons);
	eventTree->Branch("n_muons", &n_muons);
	eventTree->Branch("n_loose_taus", &n_loose_taus);
	eventTree->Branch("n_medium_taus", &n_medium_taus);
	eventTree->Branch("n_tight_taus", &n_tight_taus);
	eventTree->Branch("n_jets", &n_jets);
	eventTree->Branch("n_btags", &n_btags);

	eventTree->Branch("pv_x", &pv_x);
	eventTree->Branch("pv_y", &pv_y);
	eventTree->Branch("pv_z", &pv_z);
	
	//eventTree->Bronch("electrons", "std::vector<TLorentzVector>", &electrons);
	eventTree->Branch("e_pt", &e_pt);
	eventTree->Branch("e_eta", &e_eta);
	eventTree->Branch("e_phi", &e_phi);
	eventTree->Branch("e_mass", &e_mass);
	eventTree->Branch("e_charges", &e_charges);
	eventTree->Branch("e_vtx_dz", &e_vtx_dz);
	eventTree->Branch("e_vtx_dxy", &e_vtx_dxy);
	eventTree->Branch("e_vz", &e_vz);
	eventTree->Branch("e_vx", &e_vx);
	eventTree->Branch("e_vy", &e_vy);
	eventTree->Branch("e_isGsfCtfScPixChargeConsistent", &e_isGsfCtfScPixChargeConsistent);
	
	//eventTree->Bronch("muons", "TLorentzVector", &muons);
	eventTree->Branch("mu_pt", &mu_pt);
	eventTree->Branch("mu_eta", &mu_eta);
	eventTree->Branch("mu_phi", &mu_phi);
	eventTree->Branch("mu_mass", &mu_mass);
	eventTree->Branch("mu_charges", &mu_charges);
	eventTree->Branch("mu_vtx_dz", &mu_vtx_dz);
	eventTree->Branch("mu_vtx_dxy", &mu_vtx_dxy);
	eventTree->Branch("mu_vz", &mu_vz);
	eventTree->Branch("mu_vx", &mu_vx);
	eventTree->Branch("mu_vy", &mu_vy);
	eventTree->Branch("mu_relTrkPtError", &mu_relTrkPtError);
	
	//eventTree->Bronch("loose_taus", "TLorentzVector", &loose_taus);
	eventTree->Branch("loose_tau_pt", &loose_tau_pt);
	eventTree->Branch("loose_tau_eta", &loose_tau_eta);
	eventTree->Branch("loose_tau_phi", &loose_tau_phi);
	eventTree->Branch("loose_tau_mass", &loose_tau_mass);
	eventTree->Branch("loose_tau_charges", &Ltau_charges);
	eventTree->Branch("loose_tau_vz", &Ltau_vz);
	eventTree->Branch("loose_tau_vx", &Ltau_vx);
	eventTree->Branch("loose_tau_vy", &Ltau_vy);
	eventTree->Branch("loose_tau_vtx_dz", &Ltau_vtx_dz);
	eventTree->Branch("loose_tau_vtx_dxy", &Ltau_vtx_dxy);
	eventTree->Branch("loose_tau_decaymode", &Ltau_decaymode);
	eventTree->Branch("loose_tau_ntrk", &Ltau_ntrk);
	//eventTree->Bronch("medium_taus", "TLorentzVector", &medium_taus);
	eventTree->Branch("medium_tau_pt", &medium_tau_pt);
	eventTree->Branch("medium_tau_eta", &medium_tau_eta);
	eventTree->Branch("medium_tau_phi", &medium_tau_phi);
	eventTree->Branch("medium_tau_mass", &medium_tau_mass);
	eventTree->Branch("medium_tau_charges", &Mtau_charges);
	eventTree->Branch("medium_tau_vz", &Mtau_vz);
	eventTree->Branch("medium_tau_vx", &Mtau_vx);
	eventTree->Branch("medium_tau_vy", &Mtau_vy);
	eventTree->Branch("medium_tau_vtx_dz", &Mtau_vtx_dz);
	eventTree->Branch("medium_tau_vtx_dxy", &Mtau_vtx_dxy);
	eventTree->Branch("medium_tau_decaymode", &Mtau_decaymode);
	eventTree->Branch("medium_tau_ntrk", &Mtau_ntrk);
	//eventTree->Bronch("tight_taus", "TLorentzVector", &tight_taus);
	eventTree->Branch("tight_tau_pt", &tight_tau_pt);
	eventTree->Branch("tight_tau_eta", &tight_tau_eta);
	eventTree->Branch("tight_tau_phi", &tight_tau_phi);
	eventTree->Branch("tight_tau_mass", &tight_tau_mass);
	eventTree->Branch("tight_tau_charges", &Ttau_charges);
	eventTree->Branch("tight_tau_vz", &Ttau_vz);
	eventTree->Branch("tight_tau_vx", &Ttau_vx);
	eventTree->Branch("tight_tau_vy", &Ttau_vy);
	eventTree->Branch("tight_tau_vtx_dz", &Ttau_vtx_dz);
	eventTree->Branch("tight_tau_vtx_dxy", &Ttau_vtx_dxy);
	eventTree->Branch("tight_tau_decaymode", &Ttau_decaymode);
	eventTree->Branch("tight_tau_ntrk", &Ttau_ntrk);
	
	//eventTree->Bronch("jets", "TLorentzVector", &jets);
	eventTree->Branch("jet_pt", &jet_pt);
	eventTree->Branch("jet_eta", &jet_eta);
	eventTree->Branch("jet_phi", &jet_phi);
	eventTree->Branch("jet_mass", &jet_mass);
	eventTree->Branch("jet_charges", &jet_charges);
	eventTree->Branch("jet_vz", &jet_vz);
	eventTree->Branch("jet_vx", &jet_vx);
	eventTree->Branch("jet_vy", &jet_vy);

	//eventTree->Bronch("bjets", "TLorentzVector", &bjets);
	eventTree->Branch("bjet_pt", &bjet_pt);
	eventTree->Branch("bjet_eta", &bjet_eta);
	eventTree->Branch("bjet_phi", &bjet_phi);
	eventTree->Branch("bjet_mass", &bjet_mass);
	eventTree->Branch("bjet_charges", &bjet_charges);
	eventTree->Branch("bjet_vz", &bjet_vz);
	eventTree->Branch("bjet_vx", &bjet_vx);
	eventTree->Branch("bjet_vy", &bjet_vy);

	// MET
	eventTree->Branch("MET_x", &MET_x);
	eventTree->Branch("MET_y", &MET_y);
	eventTree->Branch("METSignificance", &METSig);
	eventTree->Branch("METCovariance00", &METCov00);
	eventTree->Branch("METCovariance01", &METCov01);
	eventTree->Branch("METCovariance10", &METCov10);
	eventTree->Branch("METCovariance11", &METCov11);
	
	eventTree->Branch("gen_x_pdgId", &gen_x_pdgId);
	eventTree->Branch("gen_x_status", &gen_x_status);
	//eventTree->Bronch("gen_x", "TLorentzVector", &gen_x);
	eventTree->Branch("gen_x_pt", &gen_x_pt);
	eventTree->Branch("gen_x_eta", &gen_x_eta);
	eventTree->Branch("gen_x_phi", &gen_x_phi);
	eventTree->Branch("gen_x_mass", &gen_x_mass);
	eventTree->Branch("gen_x_vx", &gen_x_vx);
	eventTree->Branch("gen_x_vy", &gen_x_vy);
	eventTree->Branch("gen_x_vz", &gen_x_vz);

	eventTree->Branch("gen_top_pdgId", &gen_top_pdgId);
	eventTree->Branch("gen_top_status", &gen_top_status);
	//eventTree->Bronch("gen_top", "TLorentzVector", &gen_top);
	eventTree->Branch("gen_top_pt", &gen_top_pt);
	eventTree->Branch("gen_top_eta", &gen_top_eta);
	eventTree->Branch("gen_top_phi", &gen_top_phi);
	eventTree->Branch("gen_top_mass", &gen_top_mass);
	eventTree->Branch("gen_top_vx", &gen_top_vx);
	eventTree->Branch("gen_top_vy", &gen_top_vy);
	eventTree->Branch("gen_top_vz", &gen_top_vz);
	
	eventTree->Branch("gen_x_daughter_pdgId", &gen_xDaug_pdgId);
	eventTree->Branch("gen_x_daughter_status", &gen_xDaug_status);
	//eventTree->Bronch("gen_x_daughter", "TLorentzVector", &gen_xDaug);
	eventTree->Branch("gen_x_daughter_pt", &gen_xDaug_pt);
	eventTree->Branch("gen_x_daughter_eta", &gen_xDaug_eta);
	eventTree->Branch("gen_x_daughter_phi", &gen_xDaug_phi);
	eventTree->Branch("gen_x_daughter_mass", &gen_xDaug_mass);
	eventTree->Branch("gen_x_daughter_vx", &gen_xDaug_vx);
	eventTree->Branch("gen_x_daughter_vy", &gen_xDaug_vy);
	eventTree->Branch("gen_x_daughter_vz", &gen_xDaug_vz);
	
	eventTree->Branch("gen_tau_class", &gen_tau_class);

	eventTree->Branch("gen_top_daughter_pdgId", &gen_topDaug_pdgId);
	eventTree->Branch("gen_top_daughter_status", &gen_topDaug_status);
	//eventTree->Bronch("gen_top_daughter", "TLorentzVector", &gen_topDaug);
	eventTree->Branch("gen_top_daughter_pt", &gen_topDaug_pt);
	eventTree->Branch("gen_top_daughter_eta", &gen_topDaug_eta);
	eventTree->Branch("gen_top_daughter_phi", &gen_topDaug_phi);
	eventTree->Branch("gen_top_daughter_mass", &gen_topDaug_mass);
	eventTree->Branch("gen_top_daughter_vx", &gen_topDaug_vx);
	eventTree->Branch("gen_top_daughter_vy", &gen_topDaug_vy);
	eventTree->Branch("gen_top_daughter_vz", &gen_topDaug_vz);
	
	eventTree->Branch("gen_w_daughter_pdgId", &gen_wDaug_pdgId);
	eventTree->Branch("gen_w_daughter_status", &gen_wDaug_status);
	//eventTree->Bronch("gen_w_daughter", "TLorentzVector", &gen_wDaug);
	eventTree->Branch("gen_w_daughter_pt", &gen_wDaug_pt);
	eventTree->Branch("gen_w_daughter_eta", &gen_wDaug_eta);
	eventTree->Branch("gen_w_daughter_phi", &gen_wDaug_phi);
	eventTree->Branch("gen_w_daughter_mass", &gen_wDaug_mass);
	eventTree->Branch("gen_w_daughter_vx", &gen_wDaug_vx);
	eventTree->Branch("gen_w_daughter_vy", &gen_wDaug_vy);
	eventTree->Branch("gen_w_daughter_vz", &gen_wDaug_vz);
}

#endif
