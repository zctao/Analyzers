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

#include <iostream>
#include <string>
#include <vector>

using namespace std;

//namespace Cuts
//{
//	float pT = 25.0;
//	float eta = 2.3;
//	int nJets = 2;
//	int nbTags = 1;
//}

const int nstep = 20;

void getEff(double ptmin /*=20*/, double ptmax, int nstep, double eff[4][20], TTree* tree) 
{
	// initialize arrays
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 20; ++j)
			eff[i][j] = 0;
	}
	
	const int nevt = tree->GetEntries();
	
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
	for (int ievt=0; ievt<nevt; ++ievt) {
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

	// Divided by total number of events
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < nstep; ++j){
			eff[i][j]/=nevt;
		}
	}
	
}


void CutsPerf(const TString sig_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_sig.root", 
	      const TString bkg_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_TTJets.root")
{
	// Read ntuples
	TFile* f_sig = new TFile(sig_file);
	TFile* f_bkg = new TFile(bkg_file);
	
	TTree* tree_sig = (TTree*) f_sig->Get("ttHsyncExercise/EventTree");
	TTree* tree_bkg = (TTree*) f_bkg->Get("ttHsyncExercise/EventTree");

	const int nsig = tree_sig->GetEntries();
	const int nbkg = tree_bkg->GetEntries();

	if (nsig == 0 or nbkg == 0) {
		cout << "File doesn't exist or is empty, returning..." << endl;
		cout << endl;
		return;
	}
	else {
		cout << "Number of signal events :" << nsig << endl;
		cout << "Number of background events :" << nbkg << endl;
	}

	// Define efficiency arrays
	//const int nstep = 20;
	double eff_sig[4][20];
	double eff_bkg[4][20];
	
	getEff(20, 120, nstep, eff_sig, tree_sig);
	getEff(20, 120, nstep, eff_bkg, tree_bkg);

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
	gr_loose->SetLineColor(1);
    gr_loose->SetMarkerColor(1);
	gr_medium->SetLineColor(2);
    gr_medium->SetMarkerColor(2);
	gr_tight->SetLineColor(3);
    gr_tight->SetMarkerColor(3);
	
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
