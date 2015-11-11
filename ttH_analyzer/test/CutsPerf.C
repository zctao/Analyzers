#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TPad.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;


const int nstep = 20;

void getEffArray(double ptmin, double ptmax, int nstep, double eff[4][20], TTree* tree);
void MakeROCPlot(TTree *tree_sig, int nevt_sig, TTree* tree_bkg, int nevt_bkg);
int CutHistFiller(TTree* tree, TString label);
void CutHistDrawer(TString histfile1, TString histfile2);
void CutFlowDrawer(TH1F* h_cutflow_sig, TH1F* h_cutflow_TTJets);

int nsig = -99;
int nTTJets = -99;

void CutsPerf(
							//const TString sig_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_sig.root",
							//const TString bkg_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_TTJets.root"
							const TString sig_file = "/eos/uscms/store/user/ztao/ttHToTT_M125_13TeV_powheg_pythia8/ttHToTauTau_Ntuple_signal/151022_185116/0000/CU_ttH_EDA_output.root",
							const TString bkg_file = "/eos/uscms/store/user/ztao/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ttHToTauTau_Ntuple_TTJets/151022_143356/0000/CU_ttH_EDA_output.root"
							)
{
	// Read ntuples
	TFile* f_sig = new TFile(sig_file);
	TFile* f_TTJets = new TFile(bkg_file);
	
	TTree* tree_sig = (TTree*) f_sig->Get("ttHsyncExercise/EventTree");
	TTree* tree_TTJets = (TTree*) f_TTJets->Get("ttHsyncExercise/EventTree");

	nsig = tree_sig->GetEntries();
	nTTJets = tree_TTJets->GetEntries();

	if (nsig == 0 or nTTJets == 0) {
		cout << "File doesn't exist or is empty, returning..." << endl;
		cout << endl;
		return;
	}
	else {
		cout << "Number of entries (signal):" << nsig << endl;
		cout << "Number of entries (TTJets) :" << nTTJets << endl;
	}

	// Total number of events
	// (before basic selection, not the number of events in the tree)
	// Get the total number of events from first bin of cut flow histogram

	TH1F* h_cutflow_sig =
		(TH1F*) f_sig->Get("ttHsyncExercise/h_tth_syncex_dileptauh");
	TH1F* h_cutflow_TTJets =
		(TH1F*) f_TTJets->Get("ttHsyncExercise/h_tth_syncex_dileptauh");
	
	int nevt_sig = h_cutflow_sig -> GetBinContent(1);
	int nevt_TTJets = h_cutflow_TTJets -> GetBinContent(1);

	cout << "Total number of signal sample processed :" << nevt_sig << endl;
	cout << "Total number of TTJets sample processed :" << nevt_TTJets << endl;
	
	
	//MakeROCPlot(tree_sig, nevt_sig, tree_TTJets, nevt_TTJets);
	
	nsig = CutHistFiller(tree_sig, "sig");
	nTTJets = CutHistFiller(tree_TTJets, "TTJets");
	
	// Draw Histograms
	CutHistDrawer("/uscms/home/ztao/work/CU_ttH_WD/Outputs/cuts_hist_sig.root", "/uscms/home/ztao/work/CU_ttH_WD/Outputs/cuts_hist_TTJets.root");
	CutFlowDrawer(h_cutflow_sig, h_cutflow_TTJets);
	
}


// ------------------------------------------------
	
void getEffArray(double ptmin /*=20*/, double ptmax, int nstep, double eff[4][20], TTree* tree) 
{
	// initialize arrays
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 20; ++j)
			eff[i][j] = 0;
	}
	
	const int nEntries = tree->GetEntries();
	
	// Define leafs and branches

	std::vector<float>* tau_pt_noniso;
	std::vector<float>* tau_eta_noniso;
	std::vector<float>* tau_pt_loose;
	std::vector<float>* tau_eta_loose;
	std::vector<float>* tau_pt_medium;
	std::vector<float>* tau_eta_medium;
	std::vector<float>* tau_pt_tight;
	std::vector<float>* tau_eta_tight;

	TBranch* b_tau_pt_noniso;
	TBranch* b_tau_eta_noniso;
	TBranch* b_tau_pt_loose;
	TBranch* b_tau_eta_loose;
	TBranch* b_tau_pt_medium;
	TBranch* b_tau_eta_medium;
	TBranch* b_tau_pt_tight;
	TBranch* b_tau_eta_tight;

	tau_pt_noniso = 0;
	tau_eta_noniso = 0;
	tau_pt_loose = 0;
	tau_eta_loose = 0;
	tau_pt_medium = 0;
	tau_eta_medium = 0;
	tau_pt_tight = 0;
	tau_eta_tight = 0;

	tree->SetBranchAddress("tau_pt_noniso", &tau_pt_noniso, &b_tau_pt_noniso);
	tree->SetBranchAddress("tau_eta_noniso", &tau_eta_noniso, &b_tau_eta_noniso);
	tree->SetBranchAddress("tau_pt_loose", &tau_pt_loose, &b_tau_pt_loose);
	tree->SetBranchAddress("tau_eta_loose", &tau_eta_loose, &b_tau_eta_loose);
	tree->SetBranchAddress("tau_pt_medium", &tau_pt_medium, &b_tau_pt_medium);
	tree->SetBranchAddress("tau_eta_medium", &tau_eta_medium, &b_tau_eta_medium);
	tree->SetBranchAddress("tau_pt_tight", &tau_pt_tight, &b_tau_pt_tight);
	tree->SetBranchAddress("tau_eta_tight", &tau_eta_tight, &b_tau_eta_tight);

	// ----------------------------------------------------------------------------------------------------------------
	//        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
	// ----------------------------------------------------------------------------------------------------------------
	
	// event loop
	for (int ievt=0; ievt<nEntries; ++ievt) {
		tree -> GetEntry(ievt);

		for (int i = 0; i<nstep; ++i) {
			double ptcut = ptmin+(ptmax-ptmin)/nstep*i;

			// tight
			for (auto & taupt : *tau_pt_tight) {
				if (taupt < ptcut)
					continue;
				
				eff[3][i]++;
				break;  // find at least one tau that passes pt cut
				
			}

			// medium
			for (auto & taupt : *tau_pt_medium) {
				if (taupt < ptcut)
					continue;
				
				eff[2][i]++;
				break;  // find at least one tau that passes pt cut
				
			}

			// loose
			for (auto & taupt : *tau_pt_loose) {
				if (taupt < ptcut)
					continue;
				
				eff[1][i]++;
				break;  // find at least one tau that passes pt cut
				
			}

			// nonIso
			for (auto & taupt : *tau_pt_noniso) {
				if (taupt < ptcut)
					continue;
				
				eff[0][i]++;
				break;  // find at least one tau that passes pt cut
				
			}
			
		} // end of pt_cut loop

		
	} // end of event loop
	
	delete tree;
	
}

void MakeROCPlot(TTree *tree_sig, int nevt_sig, TTree* tree_bkg, int nevt_bkg)
{ // FIX ME
	
	// Define efficiency arrays
	//const int nstep = 20;
	double eff_sig[4][20];
	double eff_bkg[4][20];
	
	getEffArray(20, 120, nstep, eff_sig, tree_sig);
	getEffArray(20, 120, nstep, eff_bkg, tree_bkg);
	
	
	
	for (int i = 1; i < 4; ++i) {
		for (int j = 0; j < nstep; ++j){
			eff_sig[i][j]/=nevt_sig;
		}
	}

	for (int i = 1; i < 4; ++i) {
		for (int j = 0; j < nstep; ++j){
			eff_bkg[i][j]/=nevt_bkg;
		}
	}

	TCanvas c;
	gPad->SetGrid();
	gPad->SetLogy();

	//TGraphErrors* gr_noniso =
	//new TGraphErrors(nstep, eff_sig[0], eff_bkg[0], 0, 0);
	TGraphErrors* gr_loose =
		new TGraphErrors(nstep, eff_sig[1], eff_bkg[1], 0, 0);
	TGraphErrors* gr_medium =
		new TGraphErrors(nstep, eff_sig[2], eff_bkg[2], 0, 0);
	TGraphErrors* gr_tight =
		new TGraphErrors(nstep, eff_sig[3], eff_bkg[3], 0, 0);

	//gr_noniso->SetLineColor(4);
	//gr_noniso->SetMarkerColor(4);
	gr_loose->SetLineColor(3);
    gr_loose->SetMarkerColor(3);
	gr_medium->SetLineColor(2);
    gr_medium->SetMarkerColor(2);
	gr_tight->SetLineColor(1);
    gr_tight->SetMarkerColor(1);
	
	gr_loose->SetTitle("Tau ID WP");
	gr_loose->GetXaxis()->SetTitle("Signal efficiency");
	gr_loose->GetYaxis()->SetTitle("Background efficiency");
	gr_loose->SetMarkerStyle(21);
	gr_medium->SetMarkerStyle(20);
	gr_tight->SetMarkerStyle(22);
	//gr_noniso->SetMarkerStyle(7);

	gr_loose->Draw("APZ0");
	gr_medium->Draw("same PZ0");
	gr_tight->Draw("same PZ0");
	//gr_noniso->Draw("same PZ0");
	
	TLegend* leg=new TLegend(0.23,0.7,0.50,0.9);
	//leg->AddEntry(gr_noniso,"nonIso","p");
	leg->AddEntry(gr_loose,"loose","p");
	leg->AddEntry(gr_medium,"medium","p");
	leg->AddEntry(gr_tight,"tight","p");
	leg->Draw("same");

	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/tauWPROC.pdf");
}

int CutHistFiller(TTree* tree, TString label) 
{
	int nevtpass = 0;
	
	const int nEntries = tree->GetEntries();
	
	// Define leafs and branches

	int n_taus_loose;
	int n_taus_medium;
	int n_taus_tight;
	int n_jets;
	int n_btags;
	std::vector<float>* e_pt;
	std::vector<float>* e_eta;
	std::vector<float>* e_phi;
	std::vector<float>* mu_pt;
	std::vector<float>* mu_eta;
	std::vector<float>* mu_phi;
	std::vector<float>* tau_pt_loose;
	std::vector<float>* tau_eta_loose;
	std::vector<float>* tau_phi_loose;
	std::vector<float>* bjet_pt;
	std::vector<float>* bjet_eta;
	std::vector<float>* bjet_phi;
		
	TBranch* b_n_taus_loose;
	TBranch* b_n_taus_medium;
	TBranch* b_n_taus_tight;
	TBranch* b_n_jets;
	TBranch* b_n_btags;
	TBranch* b_e_pt;
	TBranch* b_e_eta;
	TBranch* b_e_phi;
	TBranch* b_mu_pt;
	TBranch* b_mu_eta;
	TBranch* b_mu_phi;
	TBranch* b_tau_pt_loose;
	TBranch* b_tau_eta_loose;
	TBranch* b_tau_phi_loose;
	TBranch* b_bjet_pt;
	TBranch* b_bjet_eta;
	TBranch* b_bjet_phi;

	n_taus_loose = -99;
	n_taus_medium = -99;
	n_taus_tight = -99;
	n_jets = -99;
	n_btags = -99;
	e_pt = 0;
	e_eta = 0;
	e_phi = 0;
	mu_pt = 0;
	mu_eta = 0;
	mu_phi = 0;
	tau_pt_loose = 0;
	tau_eta_loose = 0;
	tau_phi_loose = 0;
	bjet_pt = 0;
	bjet_eta = 0;
	bjet_phi = 0;

	tree->SetBranchAddress("n_taus_loose", &n_taus_loose, &b_n_taus_loose);
	tree->SetBranchAddress("n_taus_medium", &n_taus_medium, &b_n_taus_medium);
	tree->SetBranchAddress("n_taus_tight", &n_taus_tight, &b_n_taus_tight);
	tree->SetBranchAddress("n_jets", &n_jets, &b_n_jets);
	tree->SetBranchAddress("n_btags", &n_btags, &b_n_btags);
	tree->SetBranchAddress("e_pt", &e_pt, &b_e_pt);
	tree->SetBranchAddress("e_eta", &e_eta, &b_e_eta);
	tree->SetBranchAddress("e_phi", &e_phi, &b_e_phi);
	tree->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
	tree->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
	tree->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
	tree->SetBranchAddress("tau_pt_loose", &tau_pt_loose, &b_tau_pt_loose);
	tree->SetBranchAddress("tau_eta_loose", &tau_eta_loose, &b_tau_eta_loose);
	tree->SetBranchAddress("tau_phi_loose", &tau_phi_loose, &b_tau_phi_loose);
	tree->SetBranchAddress("bjet_pt", &bjet_pt, &b_bjet_pt);
	tree->SetBranchAddress("bjet_eta", &bjet_eta, &b_bjet_eta);
	tree->SetBranchAddress("bjet_phi", &bjet_phi, &b_bjet_phi);

	// Histograms
	TH1F* h_njets = new TH1F("h_njets", "", 11, -0.5, 10.5);
	TH1F* h_njetscut = new TH1F("h_njetscut", "", 9, -0.5, 8.5);
	TH1F* h_nbtags = new TH1F("h_nbtags", "", 6, -0.5, 5.5);
	TH1F* h_nbtagscut = new TH1F("h_nbtagscut", "", 5, -0.5, 4.5);
	TH1F* h_ntauID = new TH1F("h_ntauID", "", 4, -0.5, 3.5);
	TH1F* h_lep1pt = new TH1F("h_lep1pt", "", 60, 0., 300.);
	TH1F* h_lep1eta = new TH1F("h_lep1eta", "", 50, -2.5, 2.5);
	TH1F* h_lep2pt = new TH1F("h_lep2pt", "", 60, 0., 300.);
	TH1F* h_lep2eta = new TH1F("h_lep2eta", "", 50, -2.5, 2.5);
	TH1F* h_taupt = new TH1F("h_taupt", "", 60, 0., 300.);
	TH1F* h_taueta = new TH1F("h_taueta", "", 50, -2.5, 2.5);
	TH1F* h_bjet1pt = new TH1F("h_bjet1pt", "", 100, 0., 500.);
	TH1F* h_bjet1eta = new TH1F("h_bjet1eta", "", 50, -2.5, 2.5);
	TH1F* h_dRlep1tau = new TH1F("h_dRlep1tau", "", 50, 0., 5.);
	TH1F* h_dRlep2tau = new TH1F("h_dRlep2tau", "", 50, 0., 5.);
	
	// ----------------------------------------------------------------------------------------------------------------
	//        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
	// ----------------------------------------------------------------------------------------------------------------
	
	// event loop
	for (int ievt=0; ievt<nEntries; ++ievt) {

		//std::cout << "processing event #" << ievt << endl;
		
		tree -> GetEntry(ievt);
		/*
		h_njets -> Fill(n_jets);
		h_nbtags -> Fill(n_btags);
		
		for (int j = 0; j < n_jets+1; ++j) {
			h_njetscut -> Fill(j);
		}
		
	 	for (int j = 0; j < n_btags+1; ++j) {
			h_nbtagscut -> Fill(j);
		}
		
		h_ntauID -> Fill(0);
		if (n_taus_loose < 1) continue;
		h_ntauID -> Fill(1);
		if (n_taus_medium < 1) continue;
		h_ntauID -> Fill(2);
		if (n_taus_tight < 1) continue;
		h_ntauID -> Fill(3);
		*/

		// Apply cuts on number of jets and btags
		if (n_jets < 2 or n_btags < 1)
			continue;
		// Apply cuts on number of loose tau
		if (n_taus_loose < 1)
			continue;

		++nevtpass;
			
		// combine lepton pt
		std::vector<float> *lep_pt, *lep_eta, *lep_phi;
		lep_pt = 0;
		lep_eta = 0;
		lep_phi = 0;
	
		if (e_pt->size() == 0) {
			lep_pt = mu_pt;
			lep_eta = mu_eta;
			lep_phi = mu_phi;
		}
		else if (mu_pt->size() == 0) {
			lep_pt = e_pt;
			lep_eta = e_eta;
			lep_phi = e_phi;
		}
		else {
			lep_pt = e_pt;
			lep_pt->insert(lep_pt->end(),mu_pt->begin(),mu_pt->end());
			lep_eta = e_eta;
			lep_eta->insert(lep_eta->end(),mu_eta->begin(),mu_eta->end());
			lep_phi = e_phi;
			lep_phi->insert(lep_phi->end(),mu_phi->begin(),mu_phi->end());
		}
		
		size_t ilep1=0;
		size_t ilep2=0;
		
		for (size_t i=0; i<lep_pt->size(); ++i) {
			if (lep_pt->at(i) > lep_pt->at(ilep1))
				ilep1 = i;
			if (lep_pt->at(i) > lep_pt->at(ilep2) and lep_pt->at(i) < lep_pt->at(ilep1))
				ilep2 = i;
		}
		h_lep1pt -> Fill(lep_pt->at(ilep1));
		h_lep1eta -> Fill(lep_eta->at(ilep1));
		h_lep2pt -> Fill(lep_pt->at(ilep2));
		h_lep2eta -> Fill(lep_eta->at(ilep2));

		size_t itau1=0;
		for (size_t i=0; i<tau_pt_loose->size(); ++i) {
			if (tau_pt_loose->at(i)>tau_pt_loose->at(itau1))
				itau1 = i;
		}
		h_taupt -> Fill(tau_pt_loose->at(itau1));
		h_taueta -> Fill(tau_eta_loose->at(itau1));

		size_t ibjet1=0;
		for (size_t i=0; i<bjet_pt->size(); ++i) {
			if (bjet_pt->at(i)>bjet_pt->at(ibjet1))
				ibjet1 = i;
		}
		h_bjet1pt -> Fill(bjet_pt->at(ibjet1));
		h_bjet1eta -> Fill(bjet_eta->at(ibjet1));

		float deta_l1T = lep_eta->at(ilep1)-tau_eta_loose->at(itau1);
		float dphi_l1T = TVector2::Phi_mpi_pi(lep_phi->at(ilep1)-tau_phi_loose->at(itau1));
		float dRlep1tau = TMath::Sqrt(deta_l1T*deta_l1T+dphi_l1T*dphi_l1T);
		
		float deta_l2T = lep_eta->at(ilep2)-tau_eta_loose->at(itau1);
		float dphi_l2T = TVector2::Phi_mpi_pi(lep_phi->at(ilep2)-tau_phi_loose->at(itau1));
		float dRlep2tau = TMath::Sqrt(deta_l2T*deta_l2T+dphi_l2T*dphi_l2T);

		h_dRlep1tau -> Fill(dRlep1tau);
		h_dRlep2tau -> Fill(dRlep2tau);
		
		
	} // end of event loop

	TFile *outputfile = new TFile(
				     "/uscms/home/ztao/work/CU_ttH_WD/Outputs/cuts_hist_"+label+".root",
					 "RECREATE");
	
	h_njets -> Write();
	h_njetscut -> Write();
	h_nbtags -> Write();
	h_nbtagscut -> Write();
	h_ntauID -> Write();
	h_lep1pt -> Write();
	h_lep1eta -> Write();
	h_lep2pt -> Write();
	h_lep2eta -> Write();
	h_taupt -> Write();
	h_taueta -> Write();
	h_bjet1pt -> Write();
	h_bjet1eta -> Write();
	h_dRlep1tau -> Write();
	h_dRlep2tau -> Write();

	delete outputfile;

	return nevtpass;
}


void CutHistDrawer(TString histfile1, TString histfile2)
{
	TFile* f_sig = new TFile(histfile1);
	TFile* f_TTJets = new TFile(histfile2);

	TH1F* h_njets_sig = (TH1F*)f_sig->Get("h_njets");
	//TH1F* h_njetscut_sig = (TH1F*)f_sig->Get("h_njetscut");
	TH1F* h_nbtags_sig = (TH1F*)f_sig->Get("h_nbtags");
	//TH1F* h_nbtagscut_sig = (TH1F*)f_sig->Get("h_nbtagscut"); 
	//TH1F* h_ntauID_sig = (TH1F*)f_sig->Get("h_ntauID");
	TH1F* h_lep1pt_sig = (TH1F*)f_sig->Get("h_lep1pt");
	TH1F* h_lep1eta_sig = (TH1F*)f_sig->Get("h_lep1eta");
	TH1F* h_lep2pt_sig = (TH1F*)f_sig->Get("h_lep2pt");
	TH1F* h_lep2eta_sig = (TH1F*)f_sig->Get("h_lep2eta");
	TH1F* h_taupt_sig = (TH1F*)f_sig->Get("h_taupt");
	TH1F* h_taueta_sig = (TH1F*)f_sig->Get("h_taueta");
	TH1F* h_bjet1pt_sig = (TH1F*)f_sig->Get("h_bjet1pt");
	TH1F* h_bjet1eta_sig = (TH1F*)f_sig->Get("h_bjet1eta");
	TH1F* h_dRlep1tau_sig = (TH1F*)f_sig->Get("h_dRlep1tau");
	TH1F* h_dRlep2tau_sig = (TH1F*)f_sig->Get("h_dRlep2tau");
	
	TH1F* h_njets_TTJets = (TH1F*)f_TTJets->Get("h_njets");
	//TH1F* h_njetscut_TTJets = (TH1F*)f_TTJets->Get("h_njetscut"); 
	TH1F* h_nbtags_TTJets = (TH1F*)f_TTJets->Get("h_nbtags");
	//TH1F* h_nbtagscut_TTJets = (TH1F*)f_TTJets->Get("h_nbtagscut"); 
	//TH1F* h_ntauID_TTJets = (TH1F*)f_TTJets->Get("h_ntauID");
	TH1F* h_lep1pt_TTJets = (TH1F*)f_TTJets->Get("h_lep1pt");
	TH1F* h_lep1eta_TTJets = (TH1F*)f_TTJets->Get("h_lep1eta");
	TH1F* h_lep2pt_TTJets = (TH1F*)f_TTJets->Get("h_lep2pt");
	TH1F* h_lep2eta_TTJets = (TH1F*)f_TTJets->Get("h_lep2eta");
	TH1F* h_taupt_TTJets = (TH1F*)f_TTJets->Get("h_taupt");
	TH1F* h_taueta_TTJets = (TH1F*)f_TTJets->Get("h_taueta");
	TH1F* h_bjet1pt_TTJets = (TH1F*)f_TTJets->Get("h_bjet1pt");
	TH1F* h_bjet1eta_TTJets = (TH1F*)f_TTJets->Get("h_bjet1eta");
	TH1F* h_dRlep1tau_TTJets = (TH1F*)f_TTJets->Get("h_dRlep1tau");
	TH1F* h_dRlep2tau_TTJets = (TH1F*)f_TTJets->Get("h_dRlep2tau");
	
	h_njets_sig -> Sumw2();
	//h_njetscut_sig -> Sumw2();
	h_nbtags_sig -> Sumw2();
	//h_nbtagscut_sig -> Sumw2();
	//h_ntauID_sig -> Sumw2();
	/*
	h_lep1pt_sig -> Sumw2();
	h_lep1eta_sig -> Sumw2();
	h_lep2pt_sig -> Sumw2();
	h_lep2eta_sig -> Sumw2();
	h_taupt_sig -> Sumw2();
	h_taueta_sig -> Sumw2();
	h_bjet1pt_sig -> Sumw2();
	h_bjet1eta_sig -> Sumw2();
	h_dRlep1tau_sig -> Sumw2();
	h_dRlep2tau_sig -> Sumw2();
	*/
	h_njets_TTJets -> Sumw2();
	//h_njetscut_TTJets -> Sumw2();
	h_nbtags_TTJets -> Sumw2();
	//h_nbtagscut_TTJets -> Sumw2();
	//h_ntauID_TTJets -> Sumw2();
	/*
	h_lep1pt_TTJets -> Sumw2();
	h_lep1eta_TTJets -> Sumw2();
	h_lep2pt_TTJets -> Sumw2();
	h_lep2eta_TTJets -> Sumw2();
	h_taupt_TTJets -> Sumw2();
	h_taueta_TTJets -> Sumw2();
	h_bjet1pt_TTJets -> Sumw2();
	h_bjet1eta_TTJets -> Sumw2();
	h_dRlep1tau_TTJets -> Sumw2();
	h_dRlep2tau_TTJets -> Sumw2();
	*/
		
	h_njets_sig -> Scale(1.0/nsig);
	h_nbtags_sig -> Scale(1.0/nsig);
	//h_njetscut_sig -> Scale(1.0/nsig);
	//h_nbtagscut_sig -> Scale(1.0/nsig);
	//h_ntauID_sig -> Scale(1.0/nsig);
	h_lep1pt_sig -> Scale(1.0/nsig);
	h_lep1eta_sig -> Scale(1.0/nsig);
	h_lep2pt_sig -> Scale(1.0/nsig);
	h_lep2eta_sig -> Scale(1.0/nsig);
	h_taupt_sig -> Scale(1.0/nsig);
	h_taueta_sig -> Scale(1.0/nsig);
	h_bjet1pt_sig -> Scale(1.0/nsig);
	h_bjet1eta_sig -> Scale(1.0/nsig);
	h_dRlep1tau_sig -> Scale(1.0/nsig);
	h_dRlep2tau_sig -> Scale(1.0/nsig);
	h_njets_TTJets -> Scale(1.0/nTTJets);
	h_nbtags_TTJets -> Scale(1.0/nTTJets);
	//h_njetscut_TTJets -> Scale(1.0/nTTJets);
	//h_nbtagscut_TTJets -> Scale(1.0/nTTJets);
	//h_ntauID_TTJets -> Scale(1.0/nTTJets);
	h_lep1pt_TTJets -> Scale(1.0/nTTJets);
	h_lep1eta_TTJets -> Scale(1.0/nTTJets);
	h_lep2pt_TTJets -> Scale(1.0/nTTJets);
	h_lep2eta_TTJets -> Scale(1.0/nTTJets);
	h_taupt_TTJets -> Scale(1.0/nTTJets);
	h_taueta_TTJets -> Scale(1.0/nTTJets);
	h_bjet1pt_TTJets -> Scale(1.0/nTTJets);
	h_bjet1eta_TTJets -> Scale(1.0/nTTJets);
	h_dRlep1tau_TTJets -> Scale(1.0/nTTJets);
	h_dRlep2tau_TTJets -> Scale(1.0/nTTJets);


	TCanvas c;

	h_lep1pt_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lep1pt_TTJets->GetXaxis()->SetTitle("Leading Lepton p_{T} [GeV]");
	h_lep1pt_TTJets->GetYaxis()->SetTitle("Normalized");
	h_lep1pt_TTJets->SetLineColor(kBlue);
	h_lep1pt_TTJets->Draw("");
	h_lep1pt_sig->SetLineColor(kRed);
	h_lep1pt_sig->Draw("same");
	TLegend *leg_lep1pt = new TLegend(0.72,0.25,0.85,0.38);
	leg_lep1pt -> AddEntry(h_lep1pt_sig, "signal", "L");
	leg_lep1pt -> AddEntry(h_lep1pt_TTJets, "tt+Jets", "L");
	leg_lep1pt -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lep1pt.pdf");

	h_lep1eta_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lep1eta_TTJets->GetXaxis()->SetTitle("Leading Lepton #eta");
	h_lep1eta_TTJets->GetYaxis()->SetTitle("Normalized");
	h_lep1eta_TTJets->SetLineColor(kBlue);
	h_lep1eta_TTJets->Draw("");
	h_lep1eta_sig->SetLineColor(kRed);
	h_lep1eta_sig->Draw("same");
	TLegend *leg_lep1eta = new TLegend(0.72,0.25,0.85,0.38);
	leg_lep1eta -> AddEntry(h_lep1eta_sig, "signal", "L");
	leg_lep1eta -> AddEntry(h_lep1eta_TTJets, "tt+Jets", "L");
	leg_lep1eta -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lep1eta.pdf");
	
	h_lep2pt_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lep2pt_TTJets->GetXaxis()->SetTitle("Subleading Lepton p_{T} [GeV]");
	h_lep2pt_TTJets->GetYaxis()->SetTitle("Normalized");
	h_lep2pt_TTJets->SetLineColor(kBlue);
	h_lep2pt_TTJets->Draw("");
	h_lep2pt_sig->SetLineColor(kRed);
	h_lep2pt_sig->Draw("same");
	TLegend *leg_lep2pt = new TLegend(0.72,0.25,0.85,0.38);
	leg_lep2pt -> AddEntry(h_lep2pt_sig, "signal", "L");
	leg_lep2pt -> AddEntry(h_lep2pt_TTJets, "tt+Jets", "L");
	leg_lep2pt -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lep2pt.pdf");

	h_lep2eta_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lep2eta_TTJets->GetXaxis()->SetTitle("Subleading Lepton #eta");
	h_lep2eta_TTJets->GetYaxis()->SetTitle("Normalized");
	h_lep2eta_TTJets->SetLineColor(kBlue);
	h_lep2eta_TTJets->Draw("");
	h_lep2eta_sig->SetLineColor(kRed);
	h_lep2eta_sig->Draw("same");
	TLegend *leg_lep2eta = new TLegend(0.72,0.25,0.85,0.38);
	leg_lep2eta -> AddEntry(h_lep2eta_sig, "signal", "L");
	leg_lep2eta -> AddEntry(h_lep2eta_TTJets, "tt+Jets", "L");
	leg_lep2eta -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lep2eta.pdf");
	
	h_taupt_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_taupt_TTJets->GetXaxis()->SetTitle("#tau p_{T} [GeV]");
	h_taupt_TTJets->GetYaxis()->SetTitle("Normalized");
	h_taupt_TTJets->SetLineColor(kBlue);
	h_taupt_TTJets->Draw("");
	h_taupt_sig->SetLineColor(kRed);
	h_taupt_sig->Draw("same");
	TLegend *leg_taupt = new TLegend(0.72,0.25,0.85,0.38);
	leg_taupt -> AddEntry(h_taupt_sig, "signal", "L");
	leg_taupt -> AddEntry(h_taupt_TTJets, "tt+Jets", "L");
	leg_taupt -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/taupt.pdf");

	h_taueta_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_taueta_TTJets->GetXaxis()->SetTitle("#tau #eta");
	h_taueta_TTJets->GetYaxis()->SetTitle("Normalized");
	h_taueta_TTJets->SetLineColor(kBlue);
	h_taueta_TTJets->Draw("");
	h_taueta_sig->SetLineColor(kRed);
	h_taueta_sig->Draw("same");
	TLegend *leg_taueta = new TLegend(0.72,0.25,0.85,0.38);
	leg_taueta -> AddEntry(h_taueta_sig, "signal", "L");
	leg_taueta -> AddEntry(h_taueta_TTJets, "tt+Jets", "L");
	leg_taueta -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/taueta.pdf");
	
	h_bjet1pt_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_bjet1pt_TTJets->GetXaxis()->SetTitle("Leading b-jet p_{T} [GeV]");
	h_bjet1pt_TTJets->GetYaxis()->SetTitle("Normalized");
	h_bjet1pt_TTJets->SetLineColor(kBlue);
	h_bjet1pt_TTJets->Draw("");
	h_bjet1pt_sig->SetLineColor(kRed);
	h_bjet1pt_sig->Draw("same");
	TLegend *leg_bjet1pt = new TLegend(0.72,0.25,0.85,0.38);
	leg_bjet1pt -> AddEntry(h_bjet1pt_sig, "signal", "L");
	leg_bjet1pt -> AddEntry(h_bjet1pt_TTJets, "tt+Jets", "L");
	leg_bjet1pt -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/bjet1pt.pdf");

	h_bjet1eta_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_bjet1eta_TTJets->GetXaxis()->SetTitle("Leading b-jet #eta");
	h_bjet1eta_TTJets->GetYaxis()->SetTitle("Normalized");
	h_bjet1eta_TTJets->SetLineColor(kBlue);
	h_bjet1eta_TTJets->Draw("");
	h_bjet1eta_sig->SetLineColor(kRed);
	h_bjet1eta_sig->Draw("same");
	TLegend *leg_bjet1eta = new TLegend(0.72,0.25,0.85,0.38);
	leg_bjet1eta -> AddEntry(h_bjet1eta_sig, "signal", "L");
	leg_bjet1eta -> AddEntry(h_bjet1eta_TTJets, "tt+Jets", "L");
	leg_bjet1eta -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/bjet1eta.pdf");

	h_dRlep1tau_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_dRlep1tau_TTJets->GetXaxis()->SetTitle("#DeltaR(l_{ldg},#tau_{h})");
	h_dRlep1tau_TTJets->GetYaxis()->SetTitle("Normalized");
	h_dRlep1tau_TTJets->SetLineColor(kBlue);
	h_dRlep1tau_TTJets->Draw("");
	h_dRlep1tau_sig->SetLineColor(kRed);
	h_dRlep1tau_sig->Draw("same");
	TLegend *leg_dRlep1tau = new TLegend(0.72,0.25,0.85,0.38);
	leg_dRlep1tau -> AddEntry(h_dRlep1tau_sig, "signal", "L");
	leg_dRlep1tau -> AddEntry(h_dRlep1tau_TTJets, "tt+Jets", "L");
	leg_dRlep1tau -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dRlep1tau.pdf");

	h_dRlep2tau_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_dRlep2tau_TTJets->GetXaxis()->SetTitle("#DeltaR(l_{sub-ldg},#tau_{h})");
	h_dRlep2tau_TTJets->GetYaxis()->SetTitle("Normalized");
	h_dRlep2tau_TTJets->SetLineColor(kBlue);
	h_dRlep2tau_TTJets->Draw("");
	h_dRlep2tau_sig->SetLineColor(kRed);
	h_dRlep2tau_sig->Draw("same");
	TLegend *leg_dRlep2tau = new TLegend(0.72,0.25,0.85,0.38);
	leg_dRlep2tau -> AddEntry(h_dRlep2tau_sig, "signal", "L");
	leg_dRlep2tau -> AddEntry(h_dRlep2tau_TTJets, "tt+Jets", "L");
	leg_dRlep2tau -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dRlep2tau.pdf");
	
	/// jet multiplicity
	TCanvas *c1 = new TCanvas("c1", "", 800, 800);
	TPad *pad1 = new TPad("pad1", "", 0, 0.3, 1, 1.0);
	//pad1->SetTopMargin(0.06);
	pad1->SetBottomMargin(2);
	pad1->SetGridx();
	pad1->Draw();
	pad1->cd();
	h_njets_TTJets->SetStats(0);
	h_njets_TTJets->SetLineColor(kBlue);
	h_njets_TTJets->SetMarkerStyle(20);
	h_njets_TTJets->SetMarkerColor(kBlue);
	h_njets_TTJets->GetXaxis()->SetTitle("number of jets");
	//h_njets_TTJets->SetTitleOffset(0.1);
	h_njets_TTJets->GetYaxis()->SetTitle("(Normalized)");
	h_njets_TTJets->GetYaxis()->SetTitleSize(0.04);
	h_njets_TTJets->GetYaxis()->SetTitleOffset(0.8);
	h_njets_TTJets->GetYaxis()->SetLabelSize(0.03);
	h_njets_TTJets->SetTitle("Number of jets");
	h_njets_TTJets->Draw();	
	h_njets_sig->SetLineColor(kRed);
	h_njets_sig->SetMarkerStyle(20);
	h_njets_sig->SetMarkerColor(kRed);
	h_njets_sig->Draw("same");
	
	TLegend *leg_njets = new TLegend(0.72,0.25,0.85,0.38);
	leg_njets->AddEntry(h_njets_sig, "signal","p");
	leg_njets->AddEntry(h_njets_TTJets, "tt+Jets","p");
	leg_njets->Draw("same");

	c1->cd();
	TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.01);
	pad2->SetGridx();
	pad2->Draw();
	pad2->cd();

	TH1F* h_njets_ratio = (TH1F*) h_njets_TTJets->Clone("h_njets_ratio");
	h_njets_ratio->SetLineColor(kBlack);
	h_njets_ratio->Sumw2();
	h_njets_ratio->SetStats(0);
	h_njets_ratio->Divide(h_njets_sig);
	h_njets_ratio->SetMarkerStyle(21);
	h_njets_ratio->SetMarkerColor(kBlack);
	h_njets_ratio->SetTitle("");
	h_njets_ratio->GetXaxis()->SetLabelSize(0.);
	h_njets_ratio->GetYaxis()->SetTitle("TTJets/Signal");
	h_njets_ratio->GetYaxis()->SetTitleSize(0.1);
	h_njets_ratio->GetYaxis()->SetLabelSize(0.08);
	h_njets_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_njets_ratio->GetYaxis()->SetNdivisions(505);
	h_njets_ratio->Draw("ep");

	c1->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/njets.pdf");
	
	/// number of b-jets
	TCanvas *c2 = new TCanvas("c2", "", 800, 800);
	TPad *pad3 = new TPad("pad3", "", 0, 0.3, 1, 1.0);
	//pad3->SetTopMargin(0.06);
	pad3->SetBottomMargin(2);
	pad3->SetGridx();
	pad3->Draw();
	pad3->cd();
	h_nbtags_TTJets->SetStats(0);
	h_nbtags_TTJets->SetLineColor(kBlue);
	h_nbtags_TTJets->SetMarkerStyle(20);
	h_nbtags_TTJets->SetMarkerColor(kBlue);
	h_nbtags_TTJets->GetXaxis()->SetTitle("n_btags");
	//h_nbtags_TTJets->SetTitleOffset(0.1);
	h_nbtags_TTJets->GetYaxis()->SetTitle("(Normalized)");
	h_nbtags_TTJets->GetYaxis()->SetTitleSize(0.04);
	h_nbtags_TTJets->GetYaxis()->SetTitleOffset(0.8);
	h_nbtags_TTJets->GetYaxis()->SetLabelSize(0.03);
	h_nbtags_TTJets->SetTitle("Number of b-tagged jets");
	h_nbtags_TTJets->Draw();	
	h_nbtags_sig->SetLineColor(kRed);
	h_nbtags_sig->SetMarkerStyle(20);
	h_nbtags_sig->SetMarkerColor(kRed);
	h_nbtags_sig->Draw("same");

	TLegend *leg_nbtags = new TLegend(0.72,0.25,0.85,0.38);
	leg_nbtags->AddEntry(h_nbtags_sig, "signal","p");
	leg_nbtags->AddEntry(h_nbtags_TTJets, "tt+Jets","p");
	leg_nbtags->Draw("same");

	c2->cd();
	TPad *pad4 = new TPad("pad4", "", 0, 0.05, 1, 0.3);
	pad4->SetTopMargin(0);
	pad4->SetBottomMargin(0.01);
	pad4->SetGridx();
	pad4->Draw();
	pad4->cd();

	TH1F* h_nbtags_ratio = (TH1F*) h_nbtags_TTJets->Clone("h_nbtags_ratio");
	h_nbtags_ratio->SetLineColor(kBlack);
	h_nbtags_ratio->Sumw2();
	h_nbtags_ratio->SetStats(0);
	h_nbtags_ratio->Divide(h_nbtags_sig);
	h_nbtags_ratio->SetMarkerStyle(21);
	h_nbtags_ratio->SetMarkerColor(kBlack);
	h_nbtags_ratio->SetTitle("");
	h_nbtags_ratio->GetXaxis()->SetLabelSize(0.);
	h_nbtags_ratio->GetYaxis()->SetTitle("TTJets/Signal");
	h_nbtags_ratio->GetYaxis()->SetTitleSize(0.1);
	h_nbtags_ratio->GetYaxis()->SetLabelSize(0.08);
	h_nbtags_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_nbtags_ratio->GetYaxis()->SetNdivisions(505);
	h_nbtags_ratio->Draw("ep");

	c2->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/nbtags.pdf");

	/*
	/// Tau ID
	TCanvas *c3 = new TCanvas("c3", "", 800, 800);
	TPad *pad5 = new TPad("pad5", "", 0, 0.3, 1, 1.0);
	//pad5->SetTopMargin(0.06);
	pad5->SetBottomMargin(2);
	pad5->SetGridx();
	pad5->Draw();
	pad5->cd();
	h_ntauID_sig->SetStats(0);
	h_ntauID_sig->SetLineColor(kRed);
	h_ntauID_sig->SetMarkerStyle(20);
	h_ntauID_sig->SetMarkerColor(kRed);
	h_ntauID_sig->GetXaxis()->SetTitle("TauID WP");
	h_ntauID_sig->GetXaxis()->SetBinLabel(1,"NonIso");
	h_ntauID_sig->GetXaxis()->SetBinLabel(2,"Loose");
	h_ntauID_sig->GetXaxis()->SetBinLabel(3,"Medium");
	h_ntauID_sig->GetXaxis()->SetBinLabel(4,"Tight");
	//h_ntauID_sig->SetTitleOffset(0.1);
	h_ntauID_sig->GetYaxis()->SetTitle("(Normalized)");
	h_ntauID_sig->GetYaxis()->SetTitleSize(0.04);
	h_ntauID_sig->GetYaxis()->SetTitleOffset(0.8);
	h_ntauID_sig->GetYaxis()->SetLabelSize(0.03);
	h_ntauID_sig->SetTitle(">=1 tau");
	h_ntauID_sig->SetMinimum(0);
	h_ntauID_sig->Draw();
	h_ntauID_TTJets->SetLineColor(kBlue);
	h_ntauID_TTJets->SetMarkerStyle(20);
	h_ntauID_TTJets->SetMarkerColor(kBlue);
	h_ntauID_TTJets->Draw("same");

	TLegend *leg_ntauID = new TLegend(0.72,0.27,0.85,0.40);
	leg_ntauID->AddEntry(h_ntauID_sig, "signal","p");
	leg_ntauID->AddEntry(h_ntauID_TTJets, "tt+Jets","p");
	leg_ntauID->Draw("same");

	c3->cd();
	TPad *pad6 = new TPad("pad6", "", 0, 0.05, 1, 0.3);
	pad6->SetTopMargin(0);
	pad6->SetBottomMargin(0.01);
	pad6->SetGridx();
	pad6->Draw();
	pad6->cd();

	TH1F* h_ntauID_ratio = (TH1F*) h_ntauID_TTJets->Clone("h_ntauID_ratio");
	h_ntauID_ratio->SetLineColor(kBlack);
	h_ntauID_ratio->Sumw2();
	h_ntauID_ratio->SetStats(0);
	h_ntauID_ratio->Divide(h_ntauID_sig);
	h_ntauID_ratio->SetMarkerStyle(21);
	h_ntauID_ratio->SetMarkerColor(kBlack);
	h_ntauID_ratio->SetTitle("");
	h_ntauID_ratio->GetXaxis()->SetLabelSize(0.);
	h_ntauID_ratio->GetYaxis()->SetTitle("TTJets/Signal");
	h_ntauID_ratio->GetYaxis()->SetTitleSize(0.1);
	h_ntauID_ratio->GetYaxis()->SetLabelSize(0.08);
	h_ntauID_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_ntauID_ratio->GetYaxis()->SetNdivisions(505);
	h_ntauID_ratio->Draw("ep");

	c3->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/ntauID.pdf");

	/// cut on number of jets
	TCanvas *c4 = new TCanvas("c4", "", 800, 800);
	TPad *pad7 = new TPad("pad7", "", 0, 0.3, 1, 1.0);
	//pad7->SetTopMargin(0.06);
	pad7->SetBottomMargin(2);
	pad7->SetGridx();
	pad7->Draw();
	pad7->cd();
	h_njetscut_sig->SetStats(0);
	h_njetscut_sig->SetLineColor(kRed);
	h_njetscut_sig->SetMarkerStyle(20);
	h_njetscut_sig->SetMarkerColor(kRed);
	h_njetscut_sig->GetXaxis()->SetTitle(">= n_jets");
	//h_njetscut_sig->SetTitleOffset(0.1);
	h_njetscut_sig->GetYaxis()->SetTitle("(Normalized)");
	h_njetscut_sig->GetYaxis()->SetTitleSize(0.04);
	h_njetscut_sig->GetYaxis()->SetTitleOffset(0.8);
	h_njetscut_sig->GetYaxis()->SetLabelSize(0.03);
	h_njetscut_sig->SetTitle("Cut on number of jets");
	h_njetscut_sig->Draw();
	h_njetscut_TTJets->SetLineColor(kBlue);
	h_njetscut_TTJets->SetMarkerStyle(20);
	h_njetscut_TTJets->SetMarkerColor(kBlue);
	h_njetscut_TTJets->Draw("same");

	TLegend *leg_njetscut = new TLegend(0.72,0.25,0.85,0.38);
	leg_njetscut->AddEntry(h_njetscut_sig, "signal","p");
	leg_njetscut->AddEntry(h_njetscut_TTJets, "tt+Jets","p");
	leg_njetscut->Draw("same");

	c4->cd();
	TPad *pad8 = new TPad("pad8", "", 0, 0.05, 1, 0.3);
	pad8->SetTopMargin(0);
	pad8->SetBottomMargin(0.01);
	pad8->SetGridx();
	pad8->Draw();
	pad8->cd();

	TH1F* h_njetscut_ratio = (TH1F*) h_njetscut_TTJets->Clone("h_njetscut_ratio");
	h_njetscut_ratio->SetLineColor(kBlack);
	h_njetscut_ratio->Sumw2();
	h_njetscut_ratio->SetStats(0);
	h_njetscut_ratio->Divide(h_njetscut_sig);
	h_njetscut_ratio->SetMarkerStyle(21);
	h_njetscut_ratio->SetMarkerColor(kBlack);
	h_njetscut_ratio->SetTitle("");
	h_njetscut_ratio->GetXaxis()->SetLabelSize(0.);
	h_njetscut_ratio->GetYaxis()->SetTitle("TTJets/Signal");
	h_njetscut_ratio->GetYaxis()->SetTitleSize(0.1);
	h_njetscut_ratio->GetYaxis()->SetLabelSize(0.08);
	h_njetscut_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_njetscut_ratio->GetYaxis()->SetNdivisions(505);
	h_njetscut_ratio->Draw("ep");

	c4->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/njetscut.pdf");
	
	/// Cut on number of b-jets
	TCanvas *c5 = new TCanvas("c5", "", 800, 800);
	TPad *pad9 = new TPad("pad9", "", 0, 0.3, 1, 1.0);
	//pad9->SetTopMargin(0.06);
	pad9->SetBottomMargin(2);
	pad9->SetGridx();
	pad9->Draw();
	pad9->cd();
	h_nbtagscut_sig->SetStats(0);
	h_nbtagscut_sig->SetLineColor(kRed);
	h_nbtagscut_sig->SetMarkerStyle(20);
	h_nbtagscut_sig->SetMarkerColor(kRed);
	h_nbtagscut_sig->GetXaxis()->SetTitle(">= n_btags");
	//h_nbtagscut_sig->SetTitleOffset(0.1);
	h_nbtagscut_sig->GetYaxis()->SetTitle("(Normalized)");
	h_nbtagscut_sig->GetYaxis()->SetTitleSize(0.04);
	h_nbtagscut_sig->GetYaxis()->SetTitleOffset(0.8);
	h_nbtagscut_sig->GetYaxis()->SetLabelSize(0.03);
	h_nbtagscut_sig->SetTitle("Cut on number of b-tagged jets");
	h_nbtagscut_sig->Draw();
	h_nbtagscut_TTJets->SetLineColor(kBlue);
	h_nbtagscut_TTJets->SetMarkerStyle(20);
	h_nbtagscut_TTJets->SetMarkerColor(kBlue);
	h_nbtagscut_TTJets->Draw("same");

	TLegend *leg_nbtagscut = new TLegend(0.72,0.25,0.85,0.38);
	leg_nbtagscut->AddEntry(h_nbtagscut_sig, "signal","p");
	leg_nbtagscut->AddEntry(h_nbtagscut_TTJets, "tt+Jets","p");
	leg_nbtagscut->Draw("same");

	c5->cd();
	TPad *pad10 = new TPad("pad10", "", 0, 0.05, 1, 0.3);
	pad10->SetTopMargin(0);
	pad10->SetBottomMargin(0.01);
	pad10->SetGridx();
	pad10->Draw();
	pad10->cd();

	TH1F* h_nbtagscut_ratio = (TH1F*) h_nbtagscut_TTJets->Clone("h_nbtagscut_ratio");
	h_nbtagscut_ratio->SetLineColor(kBlack);
	h_nbtagscut_ratio->Sumw2();
	h_nbtagscut_ratio->SetStats(0);
	h_nbtagscut_ratio->Divide(h_nbtagscut_sig);
	h_nbtagscut_ratio->SetMarkerStyle(21);
	h_nbtagscut_ratio->SetMarkerColor(kBlack);
	h_nbtagscut_ratio->SetTitle("");
	h_nbtagscut_ratio->GetXaxis()->SetLabelSize(0.);
	h_nbtagscut_ratio->GetYaxis()->SetTitle("TTJets/Signal");
	h_nbtagscut_ratio->GetYaxis()->SetTitleSize(0.1);
	h_nbtagscut_ratio->GetYaxis()->SetLabelSize(0.08);
	h_nbtagscut_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_nbtagscut_ratio->GetYaxis()->SetNdivisions(505);
	h_nbtagscut_ratio->Draw("ep");

	c5->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/nbtagscut.pdf");
	*/

}


void CutFlowDrawer(TH1F* h_cutflow_sig, TH1F* h_cutflow_TTJets) {

	int nsig = h_cutflow_sig -> GetBinContent(1);
	int nTTJets = h_cutflow_TTJets -> GetBinContent(1);	
	h_cutflow_sig->Sumw2();
	h_cutflow_TTJets->Sumw2();
	h_cutflow_sig->Scale(1.0/nsig);
	h_cutflow_TTJets->Scale(1.0/nTTJets);
	
	TCanvas *c0 = new TCanvas("c0", "", 800, 800);
	TPad *pad_1 = new TPad("pad_1", "", 0, 0.3, 1, 1.0);
	pad_1->SetBottomMargin(2);
	pad_1->SetGridx();
	pad_1->SetLogy();
	pad_1->Draw();
	pad_1->cd();
	h_cutflow_TTJets->SetStats(0);
	h_cutflow_TTJets->SetLineColor(kBlue);
	h_cutflow_TTJets->SetMarkerStyle(20);
	h_cutflow_TTJets->SetMarkerColor(kBlue);
	h_cutflow_TTJets->GetXaxis()->SetTitle("Cuts");
	h_cutflow_TTJets->GetYaxis()->SetTitle("(Normalized)");
	h_cutflow_TTJets->GetYaxis()->SetTitleSize(0.03);
	h_cutflow_TTJets->GetYaxis()->SetTitleOffset(0.8);
	h_cutflow_TTJets->GetYaxis()->SetLabelSize(0.02);
	h_cutflow_TTJets->SetTitle("Cut Flow");
	h_cutflow_TTJets->GetXaxis()->SetBinLabel(6, ">= 2 jets");
	h_cutflow_TTJets->GetXaxis()->SetBinLabel(7, ">= 1 b-tags");
	h_cutflow_TTJets->Draw();
	h_cutflow_sig->SetLineColor(kRed);
	h_cutflow_sig->SetMarkerStyle(20);
	h_cutflow_sig->SetMarkerColor(kRed);
	h_cutflow_sig->GetXaxis()->SetLabelSize(0.);
	h_cutflow_sig->Draw("same");

	TLegend *leg_cf = new TLegend(0.72, 0.65, 0.85, 0.78);
	leg_cf->AddEntry(h_cutflow_sig, "signal", "p");
	leg_cf->AddEntry(h_cutflow_TTJets, "tt+jets", "p");
	leg_cf->Draw("same");

	c0->cd();
	TPad *pad_2 = new TPad("pad_2", "", 0, 0.05, 1, 0.3);
	pad_2->SetTopMargin(0);
	pad_2->SetBottomMargin(0.01);
	pad_2->SetGridx();
	pad_2->SetLogy();
	pad_2->Draw();
	pad_2->cd();

	TH1F* h_cutflow_ratio = (TH1F*) h_cutflow_TTJets->Clone("h_cutflow_ratio");
	h_cutflow_ratio->SetLineColor(kBlack);
	h_cutflow_ratio->SetStats(0);
	h_cutflow_ratio->Divide(h_cutflow_sig);
	h_cutflow_ratio->SetMarkerStyle(21);
	h_cutflow_ratio->SetMarkerColor(kBlack);
	h_cutflow_ratio->SetTitle("");
	h_cutflow_ratio->GetXaxis()->SetLabelSize(0.);
	h_cutflow_ratio->GetYaxis()->SetTitle("bkg/sig");
	h_cutflow_ratio->GetYaxis()->SetTitleSize(0.1);
	h_cutflow_ratio->GetYaxis()->SetLabelSize(0.08);
	h_cutflow_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_cutflow_ratio->GetYaxis()->SetNdivisions(505);
	h_cutflow_ratio->Draw("p");

	c0->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/cutflow.pdf");
	
}
