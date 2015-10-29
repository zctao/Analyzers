#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void GenHisto (const TString);
void HistDrawer(const TString);

void GenInfoPlotter ()
{
	GenHisto("/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_sigtmp.root");
	HistDrawer("/uscms/home/ztao/work/CU_ttH_WD/Outputs/gen_histograms.root");
}

void GenHisto (const TString input_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_sigtmp.root") {
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

	TH2F* h_detatH_detatbarH_lab = new TH2F("detatH_detatbarH_lab", "Lab frame", 50, -5, 5, 50, -5, 5);
	TH2F* h_detatH_detatbarH_com = new TH2F("detatH_detatbarH_com", "COM frame", 50, -5, 5, 50, -5, 5);
	TH2F* h_dphitH_dphitbarH_lab = new TH2F("dphitH_dphitbarH_lab", "Lab frame", 50, -3.5, 3.5, 50, -3.5, 3.5);
	TH2F* h_dphitH_dphitbarH_com = new TH2F("dphitH_dphitbarH_com", "COM frame", 50, -3.5, 3.5, 50, -3.5, 3.5);
	
	TH2F* h_dphi_deta_th_lab = new TH2F("dphi_deta_th_lab","Lab frame", 50, -3.5, 3.5, 50, -5, 5);
	TH2F* h_dphi_deta_tbarh_lab = new TH2F("dphi_deta_tbarh_lab","Lab frame", 50, -3.5, 3.5, 50, -5, 5);
	TH2F* h_dphi_deta_th_com = new TH2F("dphi_deta_th_com","COM frame", 50, -3.5, 3.5, 50, -5, 5);
	TH2F* h_dphi_deta_tbarh_com = new TH2F("dphi_deta_tbarh_com","COM frame", 50, -3.5, 3.5, 50, -5, 5);
	
 	// ----------------------------------------------------------------------------------------------------------------
	//        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
	// ----------------------------------------------------------------------------------------------------------------
	int debug_cnt = 0;
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

		// check if tau
		if (fabs(gen_xDaug_pdgId->at(0))!=15 or fabs(gen_xDaug_pdgId->at(1))!=15) {
			cout << "Oooooopppppssssss" << endl;
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
		

		// Cut of higgs eta
		if (abs(gen_higgs_eta->at(0)) > 2.3) continue;

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

		higgs.SetPtEtaPhiM(gen_higgs_pt->at(0),gen_higgs_eta->at(0),gen_higgs_phi->at(0),gen_higgs_mass->at(0));
		
		TVector3 v_t, v_tbar, v_h;
		v_t = top.Vect();
		v_tbar = antitop.Vect();
		v_h = higgs.Vect();
		
		double cth_th_lab =  v_t.Dot(v_h) / TMath::Abs( v_t.Mag() * v_h.Mag() );
		double cth_tbarh_lab = v_tbar.Dot(v_h) / TMath::Abs( v_tbar.Mag() * v_h.Mag() );
		h_costH_costbarH_lab -> Fill(cth_th_lab,cth_tbarh_lab);

		
		// ---------------------------------------------------------------------
		// Debug
		if (cth_th_lab>0.95 and cth_tbarh_lab>0.95) {
			++debug_cnt;
			cout << "event #" << ievt << endl;
			cout << "cth_th_lab :" << cth_th_lab << "\t" << "cth_tbarh_lab" << cth_tbarh_lab << endl;
			cout << " " << "\t" << "Px" << "\t" << "Py" << "\t" << "Pz" << "\t" << "E" << "\t" << "eta" << "\t" << "phi" << "\t" << "cosTheta" << endl;
			cout << "t " << "\t" << top.Px() << "\t" << top.Py() << "\t" << top.Pz() << "\t" << top.E() << "\t" << top.Eta() << "\t" << top.Phi() << "\t" << top.CosTheta() << endl;
			cout << "tbar " << "\t" << antitop.Px() << "\t" << antitop.Py() << "\t" << antitop.Pz() << "\t" << antitop.E() << "\t" << antitop.Eta() << "\t" << antitop.Phi() << "\t" << antitop.CosTheta() << endl;
			cout << "higgs " << "\t" << higgs.Px() << "\t" << higgs.Py() << "\t" << higgs.Pz() << "\t" << higgs.E() << "\t" << higgs.Eta() <<  "\t" << higgs.Phi() << "\t" << higgs.CosTheta() << endl;
			//cout << "calculated :" << endl;
			//cout << "cosTheta_th_lab :" << (top.Px()*higgs.Px()+top.Py()*higgs.Py()+top.Pz()*higgs.Pz())/(TMath::Sqrt(top.Px()*top.Px()+top.Py()*top.Py()+top.Pz()*top.Pz())*TMath::Sqrt(higgs.Px()*higgs.Px()+higgs.Py()*higgs.Py()+higgs.Pz()*higgs.Pz())) << endl;
			//cout << "cosTheta_tbarh_lab :" << (antitop.Px()*higgs.Px()+antitop.Py()*higgs.Py()+antitop.Pz()*higgs.Pz())/(TMath::Sqrt(antitop.Px()*antitop.Px()+antitop.Py()*antitop.Py()+antitop.Pz()*antitop.Pz())*TMath::Sqrt(higgs.Px()*higgs.Px()+higgs.Py()*higgs.Py()+higgs.Pz()*higgs.Pz())) << endl;
		}
		// ---------------------------------------------------------------------
		
		
		double dR_th_lab = higgs.DeltaR(top);
		double dR_tbarh_lab = higgs.DeltaR(antitop);
		h_dRtH_dRtbarH_lab -> Fill(dR_th_lab, dR_tbarh_lab);

		double deta_th_lab = top.Eta()-higgs.Eta();
		double deta_tbarh_lab = antitop.Eta()-higgs.Eta();
		double dphi_th_lab = TVector2::Phi_mpi_pi(top.Phi()-higgs.Phi());
		double dphi_tbarh_lab = TVector2::Phi_mpi_pi(antitop.Phi()-higgs.Phi());

		h_dphi_deta_th_lab -> Fill(dphi_th_lab, deta_th_lab);
		h_dphi_deta_tbarh_lab -> Fill(dphi_tbarh_lab, deta_tbarh_lab);

	  h_detatH_detatbarH_lab -> Fill(deta_th_lab, deta_tbarh_lab);
		h_dphitH_dphitbarH_lab -> Fill(dphi_th_lab, dphi_tbarh_lab);
		
		// boost into ttH COM frame
		TVector3 v_com = (top+antitop+higgs).BoostVector();
		top.Boost(-v_com);
		antitop.Boost(-v_com);
		higgs.Boost(-v_com);
		
		v_t = top.Vect();
		v_tbar = antitop.Vect();
		v_h = higgs.Vect();

		double cth_th_com =  v_t.Dot(v_h) / TMath::Abs( v_t.Mag() * v_h.Mag() );
		double cth_tbarh_com = v_tbar.Dot(v_h) / TMath::Abs( v_tbar.Mag() * v_h.Mag() );
		h_costH_costbarH_com -> Fill(cth_th_com,cth_tbarh_com);

		double dR_th_com = higgs.DeltaR(top);
		double dR_tbarh_com = higgs.DeltaR(antitop);
		h_dRtH_dRtbarH_com -> Fill(dR_th_com, dR_tbarh_com);

		double deta_th_com = top.Eta()-higgs.Eta();
		double deta_tbarh_com = antitop.Eta()-higgs.Eta();
		double dphi_th_com = TVector2::Phi_mpi_pi(top.Phi()-higgs.Phi());
		double dphi_tbarh_com = TVector2::Phi_mpi_pi(antitop.Phi()-higgs.Phi());
		
		h_dphi_deta_th_com -> Fill(dphi_th_com, deta_th_com);
		h_dphi_deta_tbarh_com -> Fill(dphi_tbarh_com, deta_tbarh_com);

		h_detatH_detatbarH_com -> Fill(deta_th_lab, deta_tbarh_com);
		h_dphitH_dphitbarH_com -> Fill(dphi_th_lab, dphi_tbarh_com);
		
	} // end of event loop

	cout << "debug_cnt :" << debug_cnt << endl;
	
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

	h_detatH_detatbarH_lab -> Write();
	h_detatH_detatbarH_com -> Write();
	h_dphitH_dphitbarH_lab -> Write();
	h_dphitH_dphitbarH_com -> Write();
	
	h_dphi_deta_th_lab -> Write();
	h_dphi_deta_tbarh_lab -> Write();
	h_dphi_deta_th_com -> Write();
	h_dphi_deta_tbarh_com -> Write();	
}

void HistDrawer(const TString input =
					  "/uscms/home/ztao/work/CU_ttH_WD/Outputs/gen_histograms.root" )
{
	// Open root file and get histograms
	TFile* f = new TFile(input);
	
	TH2F* h_mtautau_mtop1 = (TH2F*)f->Get("mTTmtop1");
	TH2F* h_mtautau_mtop2 = (TH2F*)f->Get("mTTmtop2");
	TH2F* h_ptautau_mtop1 = (TH2F*)f->Get("pTTmtop1");
	TH2F* h_ptautau_mtop2 = (TH2F*)f->Get("pTTmtop2");
	TH2F* h_mtop1TT_mtop2TT = (TH2F*)f->Get("mtop1TT_mtop2TT");
	TH2F* h_mplus_mminus = (TH2F*)f->Get("mplus_mminus");
	TH2F* h_mminus_mplus = (TH2F*)f->Get("mminus_mplus");
	TH2F* h_mttTT_mplus = (TH2F*)f->Get("mttTT_mplus");
	TH2F* h_mttTT_mminus = (TH2F*)f->Get("mttTT_mminus");

	TH2F* h_ctH_ctbarH_lab = (TH2F*)f->Get("ctH_ctbarH_lab");
	TH2F* h_ctH_ctbarH_com = (TH2F*)f->Get("ctH_ctbarH_com");
	TH2F* h_dRtH_dRtbarH_lab = (TH2F*)f->Get("dRtH_dRtbarH_lab");
	TH2F* h_dRtH_dRtbarH_com = (TH2F*)f->Get("dRtH_dRtbarH_com");
	TH2F* h_detatH_detatbarH_lab = (TH2F*)f->Get("detatH_detatbarH_lab");
	TH2F* h_detatH_detatbarH_com = (TH2F*)f->Get("detatH_detatbarH_com");
	TH2F* h_dphitH_dphitbarH_lab = (TH2F*)f->Get("dphitH_dphitbarH_lab");
	TH2F* h_dphitH_dphitbarH_com = (TH2F*)f->Get("dphitH_dphitbarH_com");
	
	TH1F* h_top_higgs_dRmin = (TH1F*)f->Get("dRmin");
	TH1F* h_top_higgs_dRmax = (TH1F*)f->Get("dRmax");

	TH2F* h_dphi_deta_th_lab = (TH2F*)f->Get("dphi_deta_th_lab");
	TH2F* h_dphi_deta_tbarh_lab = (TH2F*)f->Get("dphi_deta_tbarh_lab");
	TH2F* h_dphi_deta_th_com = (TH2F*)f->Get("dphi_deta_th_com");
	TH2F* h_dphi_deta_tbarh_com = (TH2F*)f->Get("dphi_deta_tbarh_com");
	
	TCanvas c;
	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(10);
	
	h_mtautau_mtop1->GetXaxis()->SetTitle("m_{#tau#tau} [GeV]");
	h_mtautau_mtop1->GetYaxis()->SetTitle("m_{top1} [GeV]");
	h_mtautau_mtop1->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mtautau_mtop1.pdf");
	
	h_mtautau_mtop2->GetXaxis()->SetTitle("m_{#tau#tau} [GeV]");
	h_mtautau_mtop2->GetYaxis()->SetTitle("m_{top2} [GeV]");
	h_mtautau_mtop2->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mtautau_mtop2.pdf");

	h_ptautau_mtop1->GetXaxis()->SetTitle("|p_{#tau#tau}| [GeV]");
	h_ptautau_mtop1->GetYaxis()->SetTitle("m_{top1} [GeV]");
	h_ptautau_mtop1->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/ptautau_mtop1.pdf");
	
	h_ptautau_mtop2->GetXaxis()->SetTitle("|p_{#tau#tau}| [GeV]");
	h_ptautau_mtop2->GetYaxis()->SetTitle("m_{top2} [GeV]");
	h_ptautau_mtop2->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/ptautau_mtop2.pdf");

	h_mtop1TT_mtop2TT->GetXaxis()->SetTitle("m_{top1+#tau+#tau} [GeV]");
	h_mtop1TT_mtop2TT->GetYaxis()->SetTitle("m_{top2+#tau+#tau} [GeV]");
	h_mtop1TT_mtop2TT->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mtop1TT_mtop2TT.pdf");

	h_mplus_mminus->GetXaxis()->SetTitle("max(m_{top1+#tau+#tau},m_{top2+#tau+#tau}) [GeV]");
	h_mplus_mminus->GetYaxis()->SetTitle("min(m_{top1+#tau+#tau},m_{top2+#tau+#tau}) [GeV]");
	h_mplus_mminus->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mplus_mminus.pdf");

	h_mminus_mplus->GetYaxis()->SetTitle("max(m_{top1+#tau+#tau},m_{top2+#tau+#tau}) [GeV]");
	h_mminus_mplus->GetXaxis()->SetTitle("min(m_{top1+#tau+#tau},m_{top2+#tau+#tau}) [GeV]");
	h_mminus_mplus->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mminus_mplus.pdf");

	h_mttTT_mplus->GetXaxis()->SetTitle("m_{top+top+#tau+#tau}");
	h_mttTT_mplus->GetYaxis()->SetTitle("max(m_{top1+#tau+#tau},m_{top2+#tau+#tau})");
	h_mttTT_mplus->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mttTT_mplus.pdf");

	h_mttTT_mminus->GetXaxis()->SetTitle("m_{top+top+#tau+#tau}");
	h_mttTT_mminus->GetYaxis()->SetTitle("min(m_{top1+#tau+#tau},m_{top2+#tau+#tau})");
	h_mttTT_mminus->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mttTT_mminus.pdf");

	h_ctH_ctbarH_lab->GetXaxis()->SetTitle("cos#theta(t,H)");
	h_ctH_ctbarH_lab->GetYaxis()->SetTitle("cos#theta(tbar,H)");
	h_ctH_ctbarH_lab->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/ctH_ctbarH_lab.pdf");

	h_ctH_ctbarH_com->GetXaxis()->SetTitle("cos#theta(t,H)");
	h_ctH_ctbarH_com->GetYaxis()->SetTitle("cos#theta(tbar,H)");
	h_ctH_ctbarH_com->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/ctH_ctbarH_com.pdf");

	h_dRtH_dRtbarH_lab->GetXaxis()->SetTitle("dR(t,H)");
	h_dRtH_dRtbarH_lab->GetYaxis()->SetTitle("dR(tbar,H)");
	h_dRtH_dRtbarH_lab->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dRtH_dRtbarH_lab.pdf");

	h_dRtH_dRtbarH_com->GetXaxis()->SetTitle("dR(t,H)");
	h_dRtH_dRtbarH_com->GetYaxis()->SetTitle("dR(tbar,H)");
	h_dRtH_dRtbarH_com->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dRtH_dRtbarH_com.pdf");

	h_detatH_detatbarH_lab->GetXaxis()->SetTitle("#Delta#eta(t,H)");
	h_detatH_detatbarH_lab->GetYaxis()->SetTitle("#Delta#eta(tbar,H)");
	h_detatH_detatbarH_lab->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/detatH_detatbarH_lab.pdf");
	
	h_dphitH_dphitbarH_lab->GetXaxis()->SetTitle("#Delta#phi(t,H)");
	h_dphitH_dphitbarH_lab->GetYaxis()->SetTitle("#Delta#phi(tbar,H)");
	h_dphitH_dphitbarH_lab->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dphitH_dphitbarH_lab.pdf");

	h_detatH_detatbarH_com->GetXaxis()->SetTitle("#Delta#eta(t,H)");
	h_detatH_detatbarH_com->GetYaxis()->SetTitle("#Delta#eta(tbar,H)");
	h_detatH_detatbarH_com->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/detatH_detatbarH_com.pdf");
	
	h_dphitH_dphitbarH_com->GetXaxis()->SetTitle("#Delta#phi(t,H)");
	h_dphitH_dphitbarH_com->GetYaxis()->SetTitle("#Delta#phi(tbar,H)");
	h_dphitH_dphitbarH_com->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dphitH_dphitbarH_com.pdf");
	
	h_dphi_deta_th_lab->GetXaxis()->SetTitle("#Delta#phi(t,H)");
	h_dphi_deta_th_lab->GetYaxis()->SetTitle("#Delta#phi(t,H)");
	h_dphi_deta_th_lab->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dphi_deta_th_lab.pdf");

	h_dphi_deta_tbarh_lab->GetXaxis()->SetTitle("#Delta#phi(tbar,H)");
	h_dphi_deta_tbarh_lab->GetYaxis()->SetTitle("#Delta#eta(tbar,H)");
	h_dphi_deta_tbarh_lab->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dphi_deta_tbarh_lab.pdf");

	h_dphi_deta_th_com->GetXaxis()->SetTitle("#Delta#phi(t,H)");
	h_dphi_deta_th_com->GetYaxis()->SetTitle("#Delta#eta(t,H)");
	h_dphi_deta_th_com->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dphi_deta_th_com.pdf");

	h_dphi_deta_tbarh_com->GetXaxis()->SetTitle("#Delta#phi(tbar,H)");
	h_dphi_deta_tbarh_com->GetYaxis()->SetTitle("#Delta#eta(tbar,H)");
	h_dphi_deta_tbarh_com->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dphi_deta_tbarh_com.pdf");
	
	gStyle->SetOptStat(1);
	h_top_higgs_dRmax->GetXaxis()->SetTitle("#DeltaR(top,Higgs)");
	h_top_higgs_dRmax->SetTitle("ttH, H->#tau#tau [GEN]");
	h_top_higgs_dRmax->SetLineColor(1);
	h_top_higgs_dRmax->Draw();
	h_top_higgs_dRmax->SetLineColor(2);
	h_top_higgs_dRmin->Draw("same");
	
	TLegend* leg = new TLegend(0.6,0.58,0.88,0.75);
	leg->AddEntry(h_top_higgs_dRmin,"min #DeltaR","l");
	leg->AddEntry(h_top_higgs_dRmax,"max #DeltaR","l");
	leg->Draw("same");
	
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dR_top_Higgs.pdf");

	
	delete f;
}

