#ifndef TreeAnalyzer_h
#define TreeAnalyzer_h

#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TLorentzVector.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "Analyzers/ttH_analyzer/interface/Types_enum.h"
#include "Analyzers/ttH_analyzer/interface/SFHelper.h"
#include "eventSelector.h"

#include <vector>
#include <algorithm>

vector<TH1D*> TreeAnalyzer(TTree* tree, bool isdata,
						   vector<vector<unsigned long long>>& eventList,
						   Analysis_types AnaType, Selection_types SelType)
{
	std::vector<TH1D*> out_hists;
	
	TTreeReader reader(tree);

	// SF Helper
	SFHelper *sf_helper = new SFHelper(AnaType, SelType, isdata);

	TTreeReaderValue<int>   rv_run(reader, "run");
	TTreeReaderValue<int>   rv_ls(reader, "ls");
	TTreeReaderValue<unsigned long long>   rv_nEvent(reader, "nEvent");
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
	TTreeReaderValue<int>   rv_mu0_ismvasel(reader, "mu0_ismvasel");
	TTreeReaderValue<int>   rv_mu0_charge(reader, "mu0_charge");
	TTreeReaderValue<float> rv_mu1_pt(reader, "mu1_pt");
	TTreeReaderValue<float> rv_mu1_conept(reader, "mu1_conept");
	TTreeReaderValue<float> rv_mu1_eta(reader, "mu1_eta");
	TTreeReaderValue<float> rv_mu1_phi(reader, "mu1_phi");
	TTreeReaderValue<float> rv_mu1_E(reader, "mu1_E");
	TTreeReaderValue<float> rv_mu1_dxy(reader, "mu1_dxy");
	TTreeReaderValue<float> rv_mu1_dz(reader, "mu1_dz");
	TTreeReaderValue<int>   rv_mu1_ismvasel(reader, "mu1_ismvasel");
	TTreeReaderValue<int>   rv_mu1_charge(reader, "mu1_charge");
	TTreeReaderValue<float> rv_ele0_pt(reader, "ele0_pt");
	TTreeReaderValue<float> rv_ele0_conept(reader, "ele0_conept");
	TTreeReaderValue<float> rv_ele0_eta(reader, "ele0_eta");
	TTreeReaderValue<float> rv_ele0_phi(reader, "ele0_phi");
	TTreeReaderValue<float> rv_ele0_E(reader, "ele0_E");
	TTreeReaderValue<float> rv_ele0_dxy(reader, "ele0_dxy");
	TTreeReaderValue<float> rv_ele0_dz(reader, "ele0_dz");
	TTreeReaderValue<int>   rv_ele0_ismvasel(reader, "ele0_ismvasel");
	TTreeReaderValue<int>   rv_ele0_charge(reader, "ele0_charge");
	TTreeReaderValue<float> rv_ele1_pt(reader, "ele1_pt");
	TTreeReaderValue<float> rv_ele1_conept(reader, "ele1_conept");
	TTreeReaderValue<float> rv_ele1_eta(reader, "ele1_eta");
	TTreeReaderValue<float> rv_ele1_phi(reader, "ele1_phi");
	TTreeReaderValue<float> rv_ele1_E(reader, "ele1_E");
	TTreeReaderValue<float> rv_ele1_dxy(reader, "ele1_dxy");
	TTreeReaderValue<float> rv_ele1_dz(reader, "ele1_dz");
	TTreeReaderValue<int>   rv_ele1_ismvasel(reader, "ele1_ismvasel");
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
	TTreeReaderValue<float> rv_npuTrue(reader, "npuTrue");
	TTreeReaderValue<float> rv_event_weight(reader, "event_weight");
	TTreeReaderValue<float> rv_PU_weight(reader, "PU_weight");
	TTreeReaderValue<float> rv_MC_weight(reader, "MC_weight");
	TTreeReaderValue<float> rv_bTagSF_weight(reader, "bTagSF_weight");
	TTreeReaderValue<float> rv_leptonSF_weight(reader, "leptonSF_weight");
	TTreeReaderValue<float> rv_tauSF_weight(reader, "tauSF_weight");
	TTreeReaderValue<float> rv_triggerSF_weight(reader, "triggerSF_weight");

	//////////////////////////////
	// define histograms here
	//TH1D* h_mass_2l_tau_leadbjet = new TH1D("massLepsTauBjet","",15,0,750);
	TH1D* h_sum_charges = new TH1D("sum_charges","",7,-3.5,3.5);
	
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
	TH1D* h_tau_eta = new TH1D("tau_eta","", 8,-2.4,2.4);
	TH1D* h_dr_lep1_tau = new TH1D("dr_lep1_tau","",10,0.,4.0);
	TH1D* h_lep1_pt = new TH1D("lep1_pt","",11,25.,300.);
	TH1D* h_lep1_eta = new TH1D("lep1_eta","",8,-2.4,2.4);
	TH1D* h_lep2_pt = new TH1D("lep2_pt","",10,10.,210.);
	TH1D* h_lep2_eta = new TH1D("lep2_eta","",8,-2.4,2.4);
	TH1D* h_mTauTauVis1 = new TH1D("mTauTauVis1","",15,0.,300.);
	TH1D* h_mTauTauVis2 = new TH1D("mTauTauVis2","",15,0.,300.);
	
	//h_mass_2l_tau_leadbjet -> Sumw2();
	h_sum_charges -> Sumw2();
	
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
	h_tau_eta -> Sumw2();
	h_dr_lep1_tau -> Sumw2();
	h_lep1_pt -> Sumw2();
	h_lep1_eta -> Sumw2();
	h_lep2_pt -> Sumw2();
	h_lep2_eta -> Sumw2();
	h_mTauTauVis1 -> Sumw2();
	h_mTauTauVis2 -> Sumw2();

	//////////////////////////////
	
	// event loop
	while (reader.Next()) {

		// HLT
		if (!(*rv_matchHLTPath)) continue;

		// for data sample only
		if (isdata) {
			std::vector<unsigned long long>
				eventid = {static_cast<unsigned long long>(*rv_run), 
						   static_cast<unsigned long long>(*rv_ls), 
						   *rv_nEvent};
			bool alreadyIncluded = find(eventList.begin(), eventList.end(),
										eventid) != eventList.end();
			if (alreadyIncluded)
				continue; // skip this event
			else
				eventList.push_back(eventid);  // update event list				
		}

		// lepton flavor category
		int lepCategory = *rv_lepCategory;
		
		// before any additional selections
		TLorentzVector lep0, lep1;
		int lep0_id = 0;
		int lep1_id = 0;
		bool lep0_istight = false;
		bool lep1_istight = false;

		if (lepCategory == 0) {// mumu
		    lep0.SetPtEtaPhiE(*rv_mu0_pt,*rv_mu0_eta,*rv_mu0_phi,*rv_mu0_E);
			lep1.SetPtEtaPhiE(*rv_mu1_pt,*rv_mu1_eta,*rv_mu1_phi,*rv_mu1_E);
			lep0_id = 13;
			lep1_id = 13;
			lep0_istight = *rv_mu0_ismvasel;
			lep1_istight = *rv_mu1_ismvasel;
		}
		else if (lepCategory == 1) {// ee
			lep0.SetPtEtaPhiE(*rv_ele0_pt,*rv_ele0_eta,*rv_ele0_phi,*rv_ele0_E);
			lep1.SetPtEtaPhiE(*rv_ele1_pt,*rv_ele1_eta,*rv_ele1_phi,*rv_ele1_E);
			lep0_id = 11;
			lep1_id = 11;
			lep0_istight = *rv_ele0_ismvasel;
			lep1_istight = *rv_ele1_ismvasel;
		}
		else if (lepCategory == 2) {// emu
			if (*rv_ele0_pt > *rv_mu0_pt) {
				lep0.SetPtEtaPhiE(*rv_ele0_pt,*rv_ele0_eta,*rv_ele0_phi,*rv_ele0_E);
				lep1.SetPtEtaPhiE(*rv_mu0_pt,*rv_mu0_eta,*rv_mu0_phi,*rv_mu0_E);
				lep0_id = 11;
				lep1_id = 13;
				lep0_istight = *rv_ele0_ismvasel;
				lep1_istight = *rv_mu0_ismvasel;
			}
			else {
				lep0.SetPtEtaPhiE(*rv_mu0_pt,*rv_mu0_eta,*rv_mu0_phi,*rv_mu0_E);
				lep1.SetPtEtaPhiE(*rv_ele0_pt,*rv_ele0_eta,*rv_ele0_phi,*rv_ele0_E);
				lep0_id = 13;
				lep1_id = 11;
				lep0_istight = *rv_mu0_ismvasel;
				lep1_istight = *rv_ele0_ismvasel;
			}
		}
		else
			assert(0);

		TLorentzVector tau;
		tau.SetPtEtaPhiE(*rv_tau0_pt,*rv_tau0_eta,*rv_tau0_phi,*rv_tau0_E);

		float weight = *rv_event_weight;
		//////////////////////////////////////////
		// update weights here if needed
		if (not isdata) {
			// update
			float w_pu = sf_helper->Get_PUWeight(*rv_npuTrue);
			float w_tausf = //*rv_tauSF_weight;
				sf_helper->Get_TauIDSF(*rv_tau0_pt,*rv_tau0_eta,*rv_isGenMatched);

			float w_lepsf = //*rv_leptonSF_weight;
				sf_helper->Get_LeptonIDSF(lep0.Pt(),lep0.Eta(),(lep0_id==11),
										(lep1_id==13), lep0_istight);
			w_lepsf *=
				sf_helper->Get_LeptonIDSF(lep1.Pt(),lep1.Eta(),(lep1_id==11),
										  (lep1_id==13), lep1_istight);
			
			float w_hlt = sf_helper->Get_HLTSF(lepCategory);

			// non trivil to update with current ntuple
			// better to get them right in Analyzer
			float w_mc = *rv_MC_weight;
			float w_btag = *rv_bTagSF_weight;

			weight = w_btag * w_mc * w_hlt * w_lepsf * w_tausf * w_pu;
		}
		//////////////////////////////////////////
		
		//h_mass_2l_tau_leadbjet -> Fill((lep0+lep1+tau+bjet).M(),weight);

		// lepton charges
		int ele0_charge = *rv_ele0_charge;
		int mu0_charge = *rv_mu0_charge;
		int tau0_charge = *rv_tau0_charge;

		// assume two leptons are same sign
		int lepCharge = lepCategory?ele0_charge:mu0_charge;
		h_sum_charges -> Fill( 2*lepCharge + tau0_charge, weight);
		
		//////////////////////////////
		// additional selections
		// Tau charge requirement
		
		if (!passTauCharge(lepCharge,tau0_charge,"control"))
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
		h_tau_eta -> Fill(*rv_tau0_eta, weight);
		h_lep1_pt -> Fill(lep0.Pt(), weight);
		h_lep1_eta -> Fill(lep0.Eta(), weight);
		h_lep2_pt -> Fill(lep1.Pt(), weight);
		h_lep2_eta -> Fill(lep1.Eta(), weight);

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
	out_hists.push_back(h_tau_pt);
	out_hists.push_back(h_tau_eta);
	out_hists.push_back(h_njet);
	out_hists.push_back(h_mindr_lep1_jet);
	out_hists.push_back(h_mindr_lep2_jet);
	out_hists.push_back(h_avg_dr_jet);
	out_hists.push_back(h_max_lep_eta);
	out_hists.push_back(h_met);
	out_hists.push_back(h_mT_met_lep1);
	out_hists.push_back(h_mht);
	out_hists.push_back(h_dr_leps);
	out_hists.push_back(h_dr_lep1_tau);
	out_hists.push_back(h_lep1_pt);
	out_hists.push_back(h_lep1_eta);
	out_hists.push_back(h_lep2_pt);
	out_hists.push_back(h_lep2_eta);
	out_hists.push_back(h_mTauTauVis1);
	out_hists.push_back(h_mTauTauVis2);

	//out_hists.push_back(h_mass_2l_tau_leadbjet);
	out_hists.push_back(h_sum_charges);

	return out_hists;
}

#endif
