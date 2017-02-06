#ifndef TreeAnalyzer_h
#define TreeAnalyzer_h

#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TLorentzVector.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "eventSelector.h"

#include <vector>
#include <algorithm>

vector<TH1D*> TreeAnalyzer(TTree* tree, //TString selection,
						   bool isdata, vector<vector<int>>& eventList)
{
	std::vector<TH1D*> out_hists;
	
	TTreeReader reader(tree);

	TTreeReaderValue<int>   rv_run(reader, "run");
	TTreeReaderValue<int>   rv_ls(reader, "ls");
	TTreeReaderValue<int>   rv_nEvent(reader, "nEvent");
	TTreeReaderValue<float> rv_event_weight(reader, "event_weight");
	TTreeReaderValue<int>   rv_matchHLTPath(reader, "matchHLTPath");
	TTreeReaderValue<int>   rv_isGenMatched(reader, "isGenMatched");
	TTreeReaderValue<int>   rv_HiggsDecayType(reader, "HiggsDecayType");
	TTreeReaderValue<int>   rv_lepCategory(reader, "lepCategory");
	// 0: mumu; 1: ee; 2: emu
	TTreeReaderValue<int>   rv_btagCategory(reader, "btagCategory");
	// 0: loose; 1: medium (>=2 medium btags)
	TTreeReaderValue<int>   rv_njet(reader, "n_presel_jet");
	TTreeReaderValue<float> rv_MT_met_lep0(reader, "MT_met_lep0");
	TTreeReaderValue<float> rv_mindr_lep0_jet(reader, "mindr_lep0_jet");
	TTreeReaderValue<float> rv_mindr_lep1_jet(reader, "mindr_lep1_jet");
	TTreeReaderValue<float> rv_lep0_conept(reader, "lep0_conept");
	TTreeReaderValue<float> rv_lep1_conept(reader, "lep1_conept");
	TTreeReaderValue<float> rv_avg_dr_jet(reader, "avg_dr_jet");
	TTreeReaderValue<float> rv_mu0_pt(reader, "mu0_pt");
	TTreeReaderValue<float> rv_mu0_conept(reader, "mu0_conept");
	TTreeReaderValue<float> rv_mu0_eta(reader, "mu0_eta");
	TTreeReaderValue<float> rv_mu0_phi(reader, "mu0_phi");
	TTreeReaderValue<float> rv_mu0_E(reader, "mu0_E");
	TTreeReaderValue<float> rv_mu0_dxy(reader, "mu0_dxy");
	TTreeReaderValue<float> rv_mu0_dz(reader, "mu0_dz");
	TTreeReaderValue<int>   rv_mu0_charge(reader, "mu0_charge");
	TTreeReaderValue<float> rv_mu1_pt(reader, "mu1_pt");
	TTreeReaderValue<float> rv_mu1_conept(reader, "mu1_conept");
	TTreeReaderValue<float> rv_mu1_eta(reader, "mu1_eta");
	TTreeReaderValue<float> rv_mu1_phi(reader, "mu1_phi");
	TTreeReaderValue<float> rv_mu1_E(reader, "mu1_E");
	TTreeReaderValue<float> rv_mu1_dxy(reader, "mu1_dxy");
	TTreeReaderValue<float> rv_mu1_dz(reader, "mu1_dz");
	TTreeReaderValue<int>   rv_mu1_charge(reader, "mu1_charge");
	TTreeReaderValue<float> rv_ele0_pt(reader, "ele0_pt");
	TTreeReaderValue<float> rv_ele0_conept(reader, "ele0_conept");
	TTreeReaderValue<float> rv_ele0_eta(reader, "ele0_eta");
	TTreeReaderValue<float> rv_ele0_phi(reader, "ele0_phi");
	TTreeReaderValue<float> rv_ele0_E(reader, "ele0_E");
	TTreeReaderValue<float> rv_ele0_dxy(reader, "ele0_dxy");
	TTreeReaderValue<float> rv_ele0_dz(reader, "ele0_dz");
	TTreeReaderValue<int>   rv_ele0_charge(reader, "ele0_charge");
	TTreeReaderValue<float> rv_ele1_pt(reader, "ele1_pt");
	TTreeReaderValue<float> rv_ele1_conept(reader, "ele1_conept");
	TTreeReaderValue<float> rv_ele1_eta(reader, "ele1_eta");
	TTreeReaderValue<float> rv_ele1_phi(reader, "ele1_phi");
	TTreeReaderValue<float> rv_ele1_E(reader, "ele1_E");
	TTreeReaderValue<float> rv_ele1_dxy(reader, "ele1_dxy");
	TTreeReaderValue<float> rv_ele1_dz(reader, "ele1_dz");
	TTreeReaderValue<int>   rv_ele1_charge(reader, "ele1_charge");
	TTreeReaderValue<float> rv_tau0_pt(reader, "tau0_pt");
	TTreeReaderValue<float> rv_tau0_eta(reader, "tau0_eta");
	TTreeReaderValue<float> rv_tau0_phi(reader, "tau0_phi");
	TTreeReaderValue<float> rv_tau0_E(reader, "tau0_E");
	TTreeReaderValue<float> rv_tau0_dxy(reader, "tau0_dxy");
	TTreeReaderValue<float> rv_tau0_dz(reader, "tau0_dz");
	TTreeReaderValue<int>   rv_tau0_charge(reader, "tau0_charge");
	TTreeReaderValue<int>   rv_tau0_mcMatchType(reader, "tau0_mcMatchType");
	TTreeReaderValue<int>   rv_njets(reader, "n_presel_jet");
	TTreeReaderValue<float> rv_jet0_pt(reader, "jet0_pt");
	TTreeReaderValue<float> rv_jet0_eta(reader, "jet0_eta");
	TTreeReaderValue<float> rv_jet0_phi(reader, "jet0_phi");
	TTreeReaderValue<float> rv_jet0_E(reader, "jet0_E");
	TTreeReaderValue<float> rv_jet0_CSV(reader, "jet0_CSV");
	TTreeReaderValue<float> rv_jet1_pt(reader, "jet1_pt");
	TTreeReaderValue<float> rv_jet1_eta(reader, "jet1_eta");
	TTreeReaderValue<float> rv_jet1_phi(reader, "jet1_phi");
	TTreeReaderValue<float> rv_jet1_E(reader, "jet1_E");
	TTreeReaderValue<float> rv_jet1_CSV(reader, "jet1_CSV");
	TTreeReaderValue<float> rv_jet2_pt(reader, "jet2_pt");
	TTreeReaderValue<float> rv_jet2_eta(reader, "jet2_eta");
	TTreeReaderValue<float> rv_jet2_phi(reader, "jet2_phi");
	TTreeReaderValue<float> rv_jet2_E(reader, "jet2_E");
	TTreeReaderValue<float> rv_jet2_CSV(reader, "jet2_CSV");
	TTreeReaderValue<float> rv_jet3_pt(reader, "jet3_pt");
	TTreeReaderValue<float> rv_jet3_eta(reader, "jet3_eta");
	TTreeReaderValue<float> rv_jet3_phi(reader, "jet3_phi");
	TTreeReaderValue<float> rv_jet3_E(reader, "jet3_E");
	TTreeReaderValue<float> rv_jet3_CSV(reader, "jet3_CSV");
	TTreeReaderValue<float> rv_met(reader, "PFMET");
	TTreeReaderValue<float> rv_mht(reader, "MHT");

	//////////////////////////////
	// define histograms here
	TH1D* h_mass_2l_tau_leadbjet = new TH1D("massLepsTauBjet","",15,0,750);

	TH1D* h_njet = new TH1D("njet","",10,2.,12.);
	TH1D* h_mindr_lep1_jet = new TH1D("mindr_lep1_jet","",12,0.4,4.0);
	TH1D* h_mindr_lep2_jet = new TH1D("mindr_lep2_jet","",12,0.4,4.0);
	TH1D* h_avg_dr_jet = new TH1D("avg_dr_jet","",10,0.,4.0);
	TH1D* h_max_lep_eta = new TH1D("max_lep_eta","",10,0,4.0);
	TH1D* h_met = new TH1D("met","",10,0.,500.);
	TH1D* h_mT_met_lep1 = new TH1D("mT_met_lep1","",10., 0., 500.);
	TH1D* h_mht = new TH1D("mht","",10,0.,500.);
	TH1D* h_dr_leps = new TH1D("dr_leps","",10,0.,4.0);
	TH1D* h_tau_pt = new TH1D("tau_pt","",11,20.,130.);
	TH1D* h_dr_lep1_tau = new TH1D("dr_lep1_tau","",10,0.,4.0);
	TH1D* h_lep1_pt = new TH1D("lep1_pt","",11,25.,300.);
	TH1D* h_lep2_pt = new TH1D("lep2_pt","",10,10.,210.);
	TH1D* h_mTauTauVis1 = new TH1D("mTauTauVis1","",15,0.,300.);
	TH1D* h_mTauTauVis2 = new TH1D("mTauTauVis2","",15,0.,300.);
	
	h_mass_2l_tau_leadbjet -> Sumw2();
	
	h_njet -> Sumw2();
	h_mindr_lep1_jet -> Sumw2();
	h_mindr_lep2_jet -> Sumw2();
	h_avg_dr_jet -> Sumw2();
	h_max_lep_eta -> Sumw2();
	h_met -> Sumw2();
	h_mT_met_lep1 -> Sumw2();
	h_mht -> Sumw2();
	h_dr_leps -> Sumw2();
	h_tau_pt -> Sumw2();
	h_dr_lep1_tau -> Sumw2();
	h_lep1_pt -> Sumw2();
	h_lep2_pt -> Sumw2();
	h_mTauTauVis1 -> Sumw2();
	h_mTauTauVis2 -> Sumw2();

	//////////////////////////////
	
	// event loop
	while (reader.Next()) {

		// HLT
		if (!(*rv_matchHLTPath)) continue;

		// for data sample only
		if (isdata) {
			std::vector<int>  eventid = {*rv_run, *rv_ls, *rv_nEvent};
			bool alreadyIncluded = find(eventList.begin(), eventList.end(),
										eventid) != eventList.end();
			if (alreadyIncluded)
				continue; // skip this event
			else
				eventList.push_back(eventid);  // update event list				
		}
		
		float weight = *rv_event_weight;

		// before any additional selections
		TLorentzVector lep0, lep1, tau, bjet;

		if (*rv_lepCategory == 0) {// mumu
		    lep0.SetPtEtaPhiE(*rv_mu0_pt,*rv_mu0_eta,*rv_mu0_phi,*rv_mu0_E);
			lep1.SetPtEtaPhiE(*rv_mu1_pt,*rv_mu1_eta,*rv_mu1_phi,*rv_mu1_E);
		}
		else if (*rv_lepCategory == 1) {// ee
			lep0.SetPtEtaPhiE(*rv_ele0_pt,*rv_ele0_eta,*rv_ele0_phi,*rv_ele0_E);
			lep1.SetPtEtaPhiE(*rv_ele1_pt,*rv_ele1_eta,*rv_ele1_phi,*rv_ele1_E);
		}
		else if (*rv_lepCategory == 2) {// emu
			if (*rv_ele0_pt > *rv_mu0_pt) {
				lep0.SetPtEtaPhiE(*rv_ele0_pt,*rv_ele0_eta,*rv_ele0_phi,*rv_ele0_E);
				lep1.SetPtEtaPhiE(*rv_mu0_pt,*rv_mu0_eta,*rv_mu0_phi,*rv_mu0_E);
			}
			else {
				lep0.SetPtEtaPhiE(*rv_mu0_pt,*rv_mu0_eta,*rv_mu0_phi,*rv_mu0_E);
				lep1.SetPtEtaPhiE(*rv_ele0_pt,*rv_ele0_eta,*rv_ele0_phi,*rv_ele0_E);
			}
		}
		else
			assert(0);

		tau.SetPtEtaPhiE(*rv_tau0_pt,*rv_tau0_eta,*rv_tau0_phi,*rv_tau0_E);

		if (*rv_jet0_CSV > 0.5426)
			bjet.SetPtEtaPhiE(*rv_jet0_pt,*rv_jet0_eta,*rv_jet0_phi,*rv_jet0_E);
		else if (*rv_jet1_CSV > 0.5426)
			bjet.SetPtEtaPhiE(*rv_jet1_pt,*rv_jet1_eta,*rv_jet1_phi,*rv_jet1_E);
		else if (*rv_jet2_CSV > 0.5426)
			bjet.SetPtEtaPhiE(*rv_jet2_pt,*rv_jet2_eta,*rv_jet2_phi,*rv_jet2_E);
		else if (*rv_jet3_CSV > 0.5426)
			bjet.SetPtEtaPhiE(*rv_jet3_pt,*rv_jet3_eta,*rv_jet3_phi,*rv_jet3_E);

		h_mass_2l_tau_leadbjet -> Fill((lep0+lep1+tau+bjet).M(),weight);

		//////////////////////////////
		// additional selections
		// Tau charge requirement
		int lepCategory = *rv_lepCategory;
		int ele0_charge = *rv_ele0_charge;
		int mu0_charge = *rv_mu0_charge;
		int tau0_charge = *rv_tau0_charge;
		
		if (!passTauCharge((lepCategory?ele0_charge:mu0_charge),tau0_charge,"control"))
			continue;
		//////////////////////////////
		
		h_njet -> Fill(*rv_njet, weight);
		h_mindr_lep1_jet -> Fill(*rv_mindr_lep0_jet, weight);
		h_mindr_lep2_jet -> Fill(*rv_mindr_lep1_jet, weight);
		h_avg_dr_jet -> Fill(*rv_avg_dr_jet, weight);	
		h_met -> Fill(*rv_met, weight);
		h_mT_met_lep1 -> Fill(*rv_MT_met_lep0, weight);
		h_mht -> Fill(*rv_mht, weight);
		h_tau_pt -> Fill(*rv_tau0_pt, weight);
		h_lep1_pt -> Fill(lep0.Pt(), weight);
		h_lep2_pt -> Fill(lep1.Pt(), weight);

		float max_lep_eta = max(abs(lep0.Eta()),abs(lep1.Eta()));
		h_max_lep_eta -> Fill(max_lep_eta, weight);
		
		float dr_leps =
			reco::deltaR(lep0.Eta(),lep0.Phi(),lep1.Eta(),lep1.Phi());
		float dr_lep1_tau =
			reco::deltaR(lep0.Eta(),lep0.Phi(),tau.Eta(),tau.Phi());	
		h_dr_leps -> Fill(dr_leps, weight);
		h_dr_lep1_tau -> Fill(dr_lep1_tau, weight);
		
		h_mTauTauVis1 -> Fill( (lep0+tau).M(), weight);
		h_mTauTauVis2 -> Fill( (lep1+tau).M(), weight);
		
	} // end of event loop

	// push_back histograms into output vector
	out_hists.push_back(h_mass_2l_tau_leadbjet);

	out_hists.push_back(h_njet);
	out_hists.push_back(h_mindr_lep1_jet);
	out_hists.push_back(h_mindr_lep2_jet);
	out_hists.push_back(h_avg_dr_jet);
	out_hists.push_back(h_max_lep_eta);
	out_hists.push_back(h_met);
	out_hists.push_back(h_mT_met_lep1);
	out_hists.push_back(h_mht);
	out_hists.push_back(h_dr_leps);
	out_hists.push_back(h_tau_pt);
	out_hists.push_back(h_dr_lep1_tau);
	out_hists.push_back(h_lep1_pt);
	out_hists.push_back(h_lep2_pt);
	out_hists.push_back(h_mTauTauVis1);
	out_hists.push_back(h_mTauTauVis2);

	return out_hists;
}

#endif
