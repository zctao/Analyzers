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

#include <iostream>
#include <string>
#include <vector>

using namespace std;


const int nstep = 20;

void getEffArray(double ptmin, double ptmax, int nstep, double eff[4][20], TTree* tree);
void MakeROCPlot(TTree *tree_sig, int nevt_sig, TTree* tree_bkg, int nevt_bkg);
void CutHistFiller(TTree* tree, TH1F* h_njets, TH1F* h_nbtags, TH1F* h_ntauID);
void CutHistDrawer(TString histfile);

int nsig = -99;
int nTTJets = -99;

void CutsPerf(const TString sig_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_sig.root", 
	      const TString bkg_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_TTJets.root")
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
		cout << "Number of selected signal events :" << nsig << endl;
		cout << "Number of selected background events :" << nTTJets << endl;
	}

	// Total number of events
	// (before basic selection, not the number of events in the tree)
	// Get the total number of events from first bin of cut flow histogram

	TH1F* h_cutflow_sig =
		(TH1F*) f_sig->Get("ttHsyncExercise/h_tth_syncex_dileptauh");
	TH1F* h_cutflow_bkg =
		(TH1F*) f_TTJets->Get("ttHsyncExercise/h_tth_syncex_dileptauh");
	
	int nevt_sig = h_cutflow_sig -> GetBinContent(1);
	int nevt_TTJets = h_cutflow_bkg -> GetBinContent(1);

	//MakeROCPlot(tree_sig, nevt_sig, tree_TTJets, nevt_TTJets);
	
	// Histograms
	TH1F* h_njets_sig = new TH1F("h_njets_sig", "", 7, -0.5, 6.5);
	TH1F* h_nbtags_sig = new TH1F("h_nbtags_sig", "", 5, -0.5, 4.5);
	TH1F* h_ntauID_sig = new TH1F("h_ntauID_sig", "", 4, -0.5, 3.5);

	TH1F* h_njets_TTJets = new TH1F("h_njets_TTJets", "", 7, -0.5, 6.5);
	TH1F* h_nbtags_TTJets = new TH1F("h_nbtags_TTJets", "", 5, -0.5, 4.5);
	TH1F* h_ntauID_TTJets = new TH1F("h_ntauID_TTJets", "", 4, -0.5, 3.5);

	CutHistFiller(tree_sig, h_njets_sig, h_nbtags_sig, h_ntauID_sig);
	CutHistFiller(tree_TTJets, h_njets_TTJets, h_nbtags_TTJets, h_ntauID_TTJets);

	TFile *outputfile = new TFile(
					 "/uscms/home/ztao/work/CU_ttH_WD/Outputs/cuts_histograms.root",
					 "RECREATE");
	
	h_njets_sig -> Write();
	h_nbtags_sig -> Write();
	h_ntauID_sig -> Write();
	h_njets_TTJets -> Write();
	h_nbtags_TTJets -> Write();
	h_ntauID_TTJets -> Write();

	delete outputfile;
	
	// Draw Histograms
	CutHistDrawer("/uscms/home/ztao/work/CU_ttH_WD/Outputs/cuts_histograms.root");
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
	
}

void MakeROCPlot(TTree *tree_sig, int nevt_sig, TTree* tree_bkg, int nevt_bkg)
{
	
	// Define efficiency arrays
	//const int nstep = 20;
	double eff_sig[4][20];
	double eff_bkg[4][20];
	
	getEffArray(20, 120, nstep, eff_sig, tree_sig);
	getEffArray(20, 120, nstep, eff_bkg, tree_bkg);
	
	
	
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < nstep; ++j){
			eff_sig[i][j]/=nevt_sig;
		}
	}

	for (int i = 0; i < 4; ++i) {
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

void CutHistFiller(TTree* tree, TH1F* h_njets, TH1F* h_nbtags, TH1F* h_ntauID) 
{
	
	const int nEntries = tree->GetEntries();
	
	// Define leafs and branches

	int n_taus_loose;
	int n_taus_medium;
	int n_taus_tight;
	int n_jets;
	int n_btags;

	TBranch* b_n_taus_loose;
	TBranch* b_n_taus_medium;
	TBranch* b_n_taus_tight;
	TBranch* b_n_jets;
	TBranch* b_n_btags;

	n_taus_loose = -99;
	n_taus_medium = -99;
	n_taus_tight = -99;
	n_jets = -99;
	n_btags = -99;

	tree->SetBranchAddress("n_taus_loose", &n_taus_loose, &b_n_taus_loose);
	tree->SetBranchAddress("n_taus_medium", &n_taus_medium, &b_n_taus_medium);
	tree->SetBranchAddress("n_taus_tight", &n_taus_tight, &b_n_taus_tight);
	tree->SetBranchAddress("n_jets", &n_jets, &b_n_jets);
	tree->SetBranchAddress("n_btags", &n_btags, &b_n_btags);
	
	// ----------------------------------------------------------------------------------------------------------------
	//        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
	// ----------------------------------------------------------------------------------------------------------------
	
	// event loop
	for (int ievt=0; ievt<nEntries; ++ievt) {
		tree -> GetEntry(ievt);

		h_njets -> Fill(n_jets);
		h_nbtags -> Fill(n_btags);

		h_ntauID -> Fill(0);
		if (n_taus_loose < 1) continue;
		h_ntauID -> Fill(1);
		if (n_taus_medium < 1) continue;
		h_ntauID -> Fill(2);
		if (n_taus_tight < 1) continue;
		h_ntauID -> Fill(3);
		
	} // end of event loop
	
}


void CutHistDrawer(TString histfile)
{
	TFile* f = new TFile(histfile);

	TH1F* h_njets_sig = (TH1F*)f->Get("h_njets_sig"); 
	TH1F* h_nbtags_sig = (TH1F*)f->Get("h_nbtags_sig"); 
	TH1F* h_ntauID_sig = (TH1F*)f->Get("h_ntauID_sig"); 
	TH1F* h_njets_TTJets = (TH1F*)f->Get("h_njets_TTJets"); 
	TH1F* h_nbtags_TTJets = (TH1F*)f->Get("h_nbtags_TTJets"); 
	TH1F* h_ntauID_TTJets = (TH1F*)f->Get("h_ntauID_TTJets");

	h_njets_sig -> Sumw2();
	h_nbtags_sig -> Sumw2();
	h_ntauID_sig -> Sumw2();
	h_njets_TTJets -> Sumw2();
	h_nbtags_TTJets -> Sumw2();
	h_ntauID_TTJets -> Sumw2();
	
	h_njets_sig -> Scale(1.0/nsig);
	h_nbtags_sig -> Scale(1.0/nsig);
	h_ntauID_sig -> Scale(1.0/nsig);
	h_njets_TTJets -> Scale(1.0/nTTJets);
	h_nbtags_TTJets -> Scale(1.0/nTTJets);
	h_ntauID_TTJets -> Scale(1.0/nTTJets);

	/// number of jets
	TCanvas *c = new TCanvas("c", "", 800, 800);
	TPad *pad1 = new TPad("pad1", "", 0, 0.3, 1, 1.0);
	//pad1->SetTopMargin(0.06);
	pad1->SetBottomMargin(2);
	pad1->SetGridx();
	pad1->Draw();
	pad1->cd();
	h_njets_sig->SetStats(0);
	h_njets_sig->SetLineColor(kRed);
	h_njets_sig->SetMarkerStyle(20);
	h_njets_sig->SetMarkerColor(kRed);
	h_njets_sig->GetXaxis()->SetTitle("n_jets");
	//h_njets_sig->SetTitleOffset(0.1);
	h_njets_sig->GetYaxis()->SetTitle("(Normalized)");
	h_njets_sig->GetYaxis()->SetTitleSize(0.04);
	h_njets_sig->GetYaxis()->SetTitleOffset(0.8);
	h_njets_sig->GetYaxis()->SetLabelSize(0.03);
	h_njets_sig->SetTitle("Number of jets");
	h_njets_sig->Draw();
	h_njets_TTJets->SetLineColor(kBlue);
	h_njets_TTJets->SetMarkerStyle(20);
	h_njets_TTJets->SetMarkerColor(kBlue);
	h_njets_TTJets->Draw("same");

	TLegend *leg_njets = new TLegend(0.72,0.25,0.85,0.38);
	leg_njets->AddEntry(h_njets_sig, "signal","p");
	leg_njets->AddEntry(h_njets_TTJets, "tt+Jets","p");
	leg_njets->Draw("same");

	c->cd();
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

	c->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/njets.pdf");
	
	/// number of b-jets
	TCanvas *c2 = new TCanvas("c2", "", 800, 800);
	TPad *pad3 = new TPad("pad3", "", 0, 0.3, 1, 1.0);
	//pad3->SetTopMargin(0.06);
	pad3->SetBottomMargin(2);
	pad3->SetGridx();
	pad3->Draw();
	pad3->cd();
	h_nbtags_sig->SetStats(0);
	h_nbtags_sig->SetLineColor(kRed);
	h_nbtags_sig->SetMarkerStyle(20);
	h_nbtags_sig->SetMarkerColor(kRed);
	h_nbtags_sig->GetXaxis()->SetTitle("n_btags");
	//h_nbtags_sig->SetTitleOffset(0.1);
	h_nbtags_sig->GetYaxis()->SetTitle("(Normalized)");
	h_nbtags_sig->GetYaxis()->SetTitleSize(0.04);
	h_nbtags_sig->GetYaxis()->SetTitleOffset(0.8);
	h_nbtags_sig->GetYaxis()->SetLabelSize(0.03);
	h_nbtags_sig->SetTitle("Number of b-tagged jets");
	h_nbtags_sig->Draw();
	h_nbtags_TTJets->SetLineColor(kBlue);
	h_nbtags_TTJets->SetMarkerStyle(20);
	h_nbtags_TTJets->SetMarkerColor(kBlue);
	h_nbtags_TTJets->Draw("same");

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
}
