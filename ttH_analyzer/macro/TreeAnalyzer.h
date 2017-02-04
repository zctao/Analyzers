#ifndef TreeAnalyzer_h
#define TreeAnalyzer_h

#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TLorentzVector.h"

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

	// define histograms here
	TH1D* h_mass_2l_tau_leadbjet = new TH1D("massLepsTauBjet","",40,0,600);
	TH1D* h_taupt = new TH1D("taupt","",18,20.,200.);
		
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
		
		if (!passTauCharge((lepCategory?ele0_charge:mu0_charge),tau0_charge,"control")) continue;
		
		h_taupt -> Fill(*rv_tau0_pt, weight);
		
	} // end of event loop

	// push_back histograms into output vector
	out_hists.push_back(h_mass_2l_tau_leadbjet);
	
	out_hists.push_back(h_taupt);

	return out_hists;
}

#endif
