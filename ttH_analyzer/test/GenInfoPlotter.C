#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void GenInfoPlotter (const TString input_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_sig.root") {
	// Read ntuples
	TFile* f = new TFile(input_file);
	TTree* tree = (TTree*) f->Get("ttHsyncExercise/EventTree");

	const int nevt = tree->GetEntries();

	if (nevt == 0) {
		cout << "File doesn't exist or is empty, returning..." << endl;
		cout << endl;
		return;
	}
	else
		cout << "Number of events :" << nevt << endl;

	// Define leafs and branches
	vector<float>* gen_top_pt;
	vector<float>* gen_top_eta;
	vector<float>* gen_top_phi;
	vector<float>* gen_top_mass;
	vector<int>* gen_top_pdgid;
	vector<float>* gen_higgs_pt;
	vector<float>* gen_higgs_eta;
	vector<float>* gen_higgs_phi;
	vector<float>* gen_higgs_mass;
	vector<int>* gen_xDaug_pdgId;
	vector<float>* gen_xDaug_pt;
	vector<float>* gen_xDaug_eta;
	vector<float>* gen_xDaug_phi;
	vector<float>* gen_xDaug_mass;

	TBranch* b_gen_top_pt;
	TBranch* b_gen_top_eta;
	TBranch* b_gen_top_phi;
	TBranch* b_gen_top_mass;
	TBranch* b_gen_top_pdgid;
	TBranch* b_gen_higgs_pt;
	TBranch* b_gen_higgs_eta;
	TBranch* b_gen_higgs_phi;
	TBranch* b_gen_higgs_mass;
	TBranch* b_gen_xDaug_pdgId;
	TBranch* b_gen_xDaug_pt;
	TBranch* b_gen_xDaug_eta;
	TBranch* b_gen_xDaug_phi;
	TBranch* b_gen_xDaug_mass;

	gen_top_pt = 0;
	gen_top_eta = 0;
	gen_top_phi = 0;
	gen_top_mass = 0;
	gen_top_pdgid = 0;
	gen_higgs_pt = 0;
	gen_higgs_eta = 0;
	gen_higgs_phi = 0;
	gen_higgs_mass = 0;
	gen_xDaug_pdgId = 0;
	gen_xDaug_pt = 0;
	gen_xDaug_eta = 0;
	gen_xDaug_phi = 0;
	gen_xDaug_mass = 0;

	tree->SetBranchAddress("gen_top_pt", &gen_top_pt, &b_gen_top_pt);
	tree->SetBranchAddress("gen_top_eta", &gen_top_eta, &b_gen_top_eta);
	tree->SetBranchAddress("gen_top_phi", &gen_top_phi, &b_gen_top_phi);
	tree->SetBranchAddress("gen_top_mass", &gen_top_mass, &b_gen_top_mass);
	tree->SetBranchAddress("gen_top_pdgId", &gen_top_pdgid, &b_gen_top_pdgid);
	tree->SetBranchAddress("gen_x_pt", &gen_higgs_pt, &b_gen_higgs_pt);
	tree->SetBranchAddress("gen_x_eta", &gen_higgs_eta, &b_gen_higgs_eta);
	tree->SetBranchAddress("gen_x_phi", &gen_higgs_phi, &b_gen_higgs_phi);
	tree->SetBranchAddress("gen_x_mass", &gen_higgs_mass, &b_gen_higgs_mass);
	tree->SetBranchAddress("gen_x_daughter_pdgId", &gen_xDaug_pdgId, &b_gen_xDaug_pdgId);
	tree->SetBranchAddress("gen_x_daughter_pt", &gen_xDaug_pt, &b_gen_xDaug_pt);
	tree->SetBranchAddress("gen_x_daughter_eta", &gen_xDaug_eta, &b_gen_xDaug_eta);
	tree->SetBranchAddress("gen_x_daughter_phi", &gen_xDaug_phi, &b_gen_xDaug_phi);
	tree->SetBranchAddress("gen_x_daughter_mass", &gen_xDaug_mass, &b_gen_xDaug_mass);

	// Histograms
	TH2F* h_mtautau_mtop1 = new TH2F("mTTmtop1","m(top1) vs m(#tau#tau)",10, 124, 126, 50, 150, 200);	
	TH2F* h_mtautau_mtop2 = new TH2F("mTTmtop2","m(top2) vs m(#tau#tau)", 10, 124, 126, 50, 150, 200);
	TH2F* h_ptautau_mtop1 = new TH2F("pTTmtop1","m(top1) vs |p(#tau#tau)|", 120, 0, 600, 50, 150, 200);	
	TH2F* h_ptautau_mtop2 = new TH2F("pTTmtop2","m(top2) vs |p(#tau#tau)|", 120, 0, 600, 50, 150, 200);

	TH2F* h_mtop1TT_mtop2TT = new TH2F("mtop1TT_mtop2TT", "m(top2+#tau+#tau) vs m(top1+#tau+#tau)", 100, 200, 1200, 100, 200, 1200);
	TH2F* h_mplus_mminus = new TH2F("mplus_mminus", "m_- vs m_+", 100, 200, 1200, 100, 200, 1200);
	TH2F* h_mminus_mplus = new TH2F("mminus_mplus", "m_+ vs m_-", 100, 200, 1200, 100, 200, 1200);

	TH2F* h_mttTT_mplus = new TH2F("mttTT_mplus", "m_+ vs m(top1+top2+#tau+#tau)", 100, 0, 2500, 100, 200, 1200);
	TH2F* h_mttTT_mminus = new TH2F("mttTT_mminus", "m_- vs m(top1+top2+#tau+#tau)", 100, 0, 2500, 100, 200, 1200);

	TH2F* h_costH_costbarH_lab = new TH2F("ctH_ctbarH_lab", "Lab frame", 50, -1, 1, 50, -1, 1);
	TH2F* h_costH_costbarH_com = new TH2F("ctH_ctbarH_com", "COM frame", 50, -1, 1, 50, -1, 1);
	TH2F* h_dRtH_dRtbarH_lab = new TH2F("dRtH_dRtbarH_lab", "Lab frame", 50, 0, 10, 50, 0, 10);
	TH2F* h_dRtH_dRtbarH_com = new TH2F("dRtH_dRtbarH_com", "COM frame", 50, 0, 10, 50, 0, 10);
	
	TH1F* h_top_higgs_dRmin = new TH1F("dRmin", "", 50, 0.0, 10.0);
	TH1F* h_top_higgs_dRmax = new TH1F("dRmax", "", 50, 0.0, 10.0);
	
	// ----------------------------------------------------------------------------------------------------------------
	//        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
	// ----------------------------------------------------------------------------------------------------------------
	
	// event loop
	for (int ievt=0; ievt<nevt; ++ievt) {
		tree -> GetEntry(ievt);
		
		TLorentzVector tau1, tau2;
		TLorentzVector top1, top2;

		int imax, imin;
		if (gen_top_mass->size()!=2) {
			cout << "We've got a problem." << endl;
			return;
		}

		if (gen_top_mass->at(0)>gen_top_mass->at(1)) {
			imax = 0;
			imin = 1;
		}
		else {
			imax = 1;
			imin = 0;
		}
		
		tau1.SetPtEtaPhiM(gen_xDaug_pt->at(0), gen_xDaug_eta->at(0),
						  gen_xDaug_phi->at(0), gen_xDaug_mass->at(0));
		tau2.SetPtEtaPhiM(gen_xDaug_pt->at(1), gen_xDaug_eta->at(1),
						  gen_xDaug_phi->at(1), gen_xDaug_mass->at(1));
		top1.SetPtEtaPhiM(gen_top_pt->at(imax), gen_top_eta->at(imax),
						  gen_top_phi->at(imax), gen_top_mass->at(imax));
		top2.SetPtEtaPhiM(gen_top_pt->at(imin), gen_top_eta->at(imin),
						  gen_top_phi->at(imin), gen_top_mass->at(imin));
		
		float mtautau, ptautau, mtop1, mtop2;	
		
		mtautau	= (tau1+tau2).M();
		ptautau = (tau1+tau2).P();		

		mtop1 = gen_top_mass->at(imax);
		mtop2 = gen_top_mass->at(imin);

		h_mtautau_mtop1 -> Fill(mtautau, mtop1);
		h_mtautau_mtop2 -> Fill(mtautau, mtop2);
		h_ptautau_mtop1 -> Fill(ptautau, mtop1);
		h_ptautau_mtop2 -> Fill(ptautau, mtop2);

		//
		float m_t1TT, m_t2TT, m_plus, m_minus;

		m_t1TT = (top1+tau1+tau2).M();
		m_t2TT = (top2+tau1+tau2).M();

		if (m_t1TT > m_t2TT) {
			m_plus = m_t1TT;
			m_minus = m_t2TT;
		}
		else {
			m_plus = m_t2TT;
			m_minus = m_t1TT;
		}

		h_mtop1TT_mtop2TT->Fill(m_t1TT, m_t2TT);
		h_mplus_mminus->Fill(m_plus, m_minus);
		h_mminus_mplus->Fill(m_minus, m_plus);
		h_mttTT_mplus->Fill((top1+top2+tau1+tau2).M(), m_plus);
		h_mttTT_mminus->Fill((top1+top2+tau1+tau2).M(), m_minus);
		
		// dR(top, Higgs)
		float dR0, dR1;

		/*TLorentzVector top1,top2, higgs;
		top1.SetPtEtaPhiM(gen_top_pt->at(0), gen_top_eta->at(0), gen_top_phi->at(0), gen_top_mass->at(0));
		top2.SetPtEtaPhiM(gen_top_pt->at(1), gen_top_eta->at(1), gen_top_phi->at(1), gen_top_mass->at(1));
		*/
		// LAB frame
		dR0 = TMath::Sqrt(
						  (gen_top_eta->at(0)-gen_higgs_eta->at(0))
						  *(gen_top_eta->at(0)-gen_higgs_eta->at(0))
						  +(gen_top_phi->at(0)-gen_higgs_phi->at(0))
						  *(gen_top_phi->at(0)-gen_higgs_phi->at(0))
						 );
		
		dR1 = TMath::Sqrt(
						  (gen_top_eta->at(1)-gen_higgs_eta->at(0))
						  *(gen_top_eta->at(1)-gen_higgs_eta->at(0))
						  +(gen_top_phi->at(1)-gen_higgs_phi->at(0))
						  *(gen_top_phi->at(1)-gen_higgs_phi->at(0))
						 );
		
		
		if (dR0 < dR1) {
			h_top_higgs_dRmin -> Fill(dR0);
			h_top_higgs_dRmax -> Fill(dR1);
		}
		else {
			h_top_higgs_dRmin -> Fill(dR1);
			h_top_higgs_dRmax -> Fill(dR0);
		}

		TLorentzVector top, antitop, higgs;

		top.SetPtEtaPhiM(gen_top_pt->at(0), gen_top_eta->at(0),
						  gen_top_phi->at(0), gen_top_mass->at(0));
		antitop.SetPtEtaPhiM(gen_top_pt->at(1), gen_top_eta->at(1),
						  gen_top_phi->at(1), gen_top_mass->at(1));

		if (gen_top_pdgid->at(0) < 0) {
			// swap top and antitop
			higgs = top;
			top = antitop;
			antitop = higgs;
		}

		higgs = tau1+tau2;

		TVector3 v_t, v_tbar, v_h;
		v_t = top.Vect();
		v_tbar = antitop.Vect();
		v_h = higgs.Vect();
		
		double cth_th =  v_t.Dot(v_h) / TMath::Abs( v_t.Mag() * v_h.Mag() );
		double cth_tbarh = v_tbar.Dot(v_h) / TMath::Abs( v_tbar.Mag() * v_h.Mag() );
		h_costH_costbarH_lab -> Fill(cth_tbarh,cth_th);

		double dR_th = higgs.DeltaR(top);
		double dR_tbarh = higgs.DeltaR(antitop);
		h_dRtH_dRtbarH_lab -> Fill(dR_th, dR_tbarh);

		// boost into ttH COM frame
		TVector3 v_com = (top+antitop+higgs).BoostVector();
		top.Boost(-v_com);
		antitop.Boost(-v_com);
		higgs.Boost(-v_com);

		v_t = top.Vect();
		v_tbar = antitop.Vect();
		v_h = higgs.Vect();

		cth_th =  v_t.Dot(v_h) / TMath::Abs( v_t.Mag() * v_h.Mag() );
		cth_tbarh = v_tbar.Dot(v_h) / TMath::Abs( v_tbar.Mag() * v_h.Mag() );
		h_costH_costbarH_com -> Fill(cth_tbarh,cth_th);

		dR_th = higgs.DeltaR(top);
		dR_tbarh = higgs.DeltaR(antitop);
		h_dRtH_dRtbarH_com -> Fill(dR_th, dR_tbarh);
		
	} // end of event loop

	TFile histfile("/uscms/home/ztao/work/CU_ttH_WD/Outputs/gen_histograms.root", "RECREATE");
	h_mtautau_mtop1 -> Write();
	h_mtautau_mtop2 -> Write();
	h_ptautau_mtop1 -> Write();
	h_ptautau_mtop2 -> Write();

	h_mtop1TT_mtop2TT -> Write();
	h_mplus_mminus -> Write();
	h_mminus_mplus -> Write();
	h_mttTT_mplus -> Write();
	h_mttTT_mminus -> Write();

	h_costH_costbarH_lab -> Write();
	h_costH_costbarH_com -> Write();
	h_dRtH_dRtbarH_lab -> Write();
	h_dRtH_dRtbarH_com -> Write();
	
	h_top_higgs_dRmin -> Write();
	h_top_higgs_dRmax -> Write();
}
