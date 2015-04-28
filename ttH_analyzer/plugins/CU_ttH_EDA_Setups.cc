#ifndef CU_ttH_EDA_Setups_cc
#define CU_ttH_EDA_Setups_cc

/// Includes
#include "CU_ttH_EDA.h"

void CU_ttH_EDA::Close_output_files()
{
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

	h_tth_syncex1_ele = fs_->make<TH1D>("h_tth_syncex1_ele", ";cut", 8, 0, 8);
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

	h_tth_syncex1_dimu = fs_->make<TH1D>("h_tth_syncex1_dimu", ";cut", 8, 0, 8);
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

/// Prepare histograms for trigger/filter counts
int CU_ttH_EDA::Set_up_Run_histograms()
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

void CU_ttH_EDA::Set_up_name_vectors()
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

	token.electrons = consumes<pat::ElectronCollection>(
		edm::InputTag(std::string("slimmedElectrons")));
	token.muons = consumes<pat::MuonCollection>(
		edm::InputTag(std::string("slimmedMuons")));

	token.jets =
		consumes<pat::JetCollection>(edm::InputTag(std::string("slimmedJets")));
	token.METs =
		consumes<pat::METCollection>(edm::InputTag(std::string("slimmedMETs")));

	token.MC_particles = consumes<reco::GenParticleCollection>(
		edm::InputTag(std::string("prunedGenParticles")));
	token.PF_candidates = consumes<pat::PackedCandidateCollection>(
		edm::InputTag(std::string("packedPFCandidates")));

	token.BS =
		consumes<reco::BeamSpot>(edm::InputTag(std::string("offlineBeamSpot")));

	token.top_jets = consumes<boosted::HEPTopJetCollection>(
		edm::InputTag("HEPTopJetsPFMatcher", "heptopjets", "p"));
	token.subfilter_jets = consumes<boosted::SubFilterJetCollection>(
		edm::InputTag("CA12JetsCA3FilterjetsPFMatcher", "subfilterjets", "p"));
}

#endif
