#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TCanvas.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

namespace Cuts
{
	const float pT = 30.0;
	const float eta = 2.1;
}

void TauEfficiency (const TString input_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output.root") {
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
	//int n_taus_loose;
	//int n_taus_medium;
	//int n_taus_tight;

	std::vector<float>* tau_pt_loose;
	std::vector<float>* tau_eta_loose;
	std::vector<float>* tau_pt_medium;
	std::vector<float>* tau_eta_medium;
	std::vector<float>* tau_pt_tight;
	std::vector<float>* tau_eta_tight;

	std::vector<float>* gen_tau_pt;
	std::vector<float>* gen_tau_eta;
	std::vector<int>* gen_tauclass;

	//TBranch* b_n_taus_loose;
	//TBranch* b_n_taus_medium;
	//TBranch* b_n_taus_tight;
	TBranch* b_tau_pt_loose;
	TBranch* b_tau_eta_loose;
	TBranch* b_tau_pt_medium;
	TBranch* b_tau_eta_medium;
	TBranch* b_tau_pt_tight;
	TBranch* b_tau_eta_tight;
	TBranch* b_gen_tau_pt;
	TBranch* b_gen_tau_eta;
	TBranch* b_gentauclass;

	//n_taus_loose = -99;
	//n_taus_medium = -99;
	//n_taus_tight = -99;
	tau_pt_loose = 0;
	tau_eta_loose = 0;
	tau_pt_medium = 0;
	tau_eta_medium = 0;
	tau_pt_tight = 0;
	tau_eta_tight = 0;

	gen_tau_pt = 0;
	gen_tau_eta = 0;
	gen_tauclass = 0;
	
	//tree->SetBranchAddress("n_taus_loose", &n_taus_loose, &b_n_taus_loose);
	//tree->SetBranchAddress("n_taus_medium", &n_taus_medium, &b_n_taus_medium);
	//tree->SetBranchAddress("n_taus_tight", &n_taus_tight, &b_n_taus_tight);
	tree->SetBranchAddress("tau_pt_loose", &tau_pt_loose, &b_tau_pt_loose);
	tree->SetBranchAddress("tau_eta_loose", &tau_eta_loose, &b_tau_eta_loose);
	tree->SetBranchAddress("tau_pt_medium", &tau_pt_medium, &b_tau_pt_medium);
	tree->SetBranchAddress("tau_eta_medium", &tau_eta_medium, &b_tau_eta_medium);
	tree->SetBranchAddress("tau_pt_tight", &tau_pt_tight, &b_tau_pt_tight);
	tree->SetBranchAddress("tau_eta_tight", &tau_eta_tight, &b_tau_eta_tight);

	tree->SetBranchAddress("gen_x_daughter_pt", &gen_tau_pt, &b_gen_tau_pt);	
	tree->SetBranchAddress("gen_x_daughter_eta", &gen_tau_eta, &b_gen_tau_eta);
	tree->SetBranchAddress("gen_tau_class", &gen_tauclass, &b_gentauclass);

	// Histograms
	TH1F* h_taus_gen = new TH1F("h_taus_gen","", 4, 0, 4);	
	TH1F* h_taus_selected_loose = new TH1F("h_taus_selected_loose", "", 4, 0, 4);
	TH1F* h_taus_selected_medium = new TH1F("h_taus_selected_medium", "", 4, 0, 4);
	TH1F* h_taus_selected_tight = new TH1F("h_taus_selected_tight", "", 4, 0, 4);

	TH1F* h_taus_gen_eta = new TH1F("h_taus_gen_eta","", 4, 0, 4);	
	TH1F* h_taus_selected_loose_eta = new TH1F("h_taus_selected_loose_eta", "", 4, 0, 4);
	TH1F* h_taus_selected_medium_eta = new TH1F("h_taus_selected_medium_eta", "", 4, 0, 4);
	TH1F* h_taus_selected_tight_eta = new TH1F("h_taus_selected_tight_eta", "", 4, 0, 4);
	
	// ----------------------------------------------------------------------------------------------------------------
	//        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
	// ----------------------------------------------------------------------------------------------------------------
	
	// event loop
	for (int ievt=0; ievt<nevt; ++ievt) {
		tree -> GetEntry(ievt);

		int n_gen = 0;
		int n_gen_eta = 0;
		int n_loose = 0;
		int n_loose_eta = 0;
		int n_medium = 0;
		int n_medium_eta = 0;
		int n_tight = 0;
		int n_tight_eta = 0;
		
		for (size_t i = 0; i < gen_tau_pt->size(); ++i) {
			if (gen_tauclass->at(i) == 0)
				continue;

			if (gen_tau_pt->at(i) > Cuts::pT) {
				n_gen++;
				if (abs(gen_tau_eta->at(i)) < Cuts::eta)
					n_gen_eta++;
			}
		}
		h_taus_gen -> Fill(n_gen);
		h_taus_gen_eta -> Fill(n_gen_eta);

		for (size_t i = 0; i< tau_pt_loose->size(); ++i) {
			if (tau_pt_loose->at(i) > Cuts::pT) {
				n_loose++;
				if (abs(tau_eta_loose->at(i)) < Cuts::eta)
					n_loose_eta++;
			}
		}
		h_taus_selected_loose -> Fill(n_loose);
		h_taus_selected_loose_eta -> Fill(n_loose_eta);

		for (size_t i = 0; i< tau_pt_medium->size(); ++i) {
			if (tau_pt_medium->at(i) > Cuts::pT) {
				n_medium++;
				if (abs(tau_eta_medium->at(i)) < Cuts::eta)
					n_medium_eta++;
			}
		}
		h_taus_selected_medium -> Fill(n_medium);
		h_taus_selected_medium_eta -> Fill(n_medium_eta);

		for (size_t i = 0; i< tau_pt_tight->size(); ++i) {
			if (tau_pt_tight->at(i) > Cuts::pT) {
				n_tight++;
				if (abs(tau_eta_tight->at(i)) < Cuts::eta)
					n_tight_eta++;
			}
		}
		h_taus_selected_tight -> Fill(n_tight);
		h_taus_selected_tight_eta -> Fill(n_tight_eta);

	} // end of event loop

	h_taus_selected_loose -> Sumw2();
	h_taus_selected_medium -> Sumw2();
	h_taus_selected_tight -> Sumw2();
	h_taus_gen -> Sumw2();
	h_taus_selected_loose_eta -> Sumw2();
	h_taus_selected_medium_eta -> Sumw2();
	h_taus_selected_tight_eta -> Sumw2();
	h_taus_gen_eta -> Sumw2();
	
	h_taus_selected_loose -> Scale(1.0/nevt);
	h_taus_selected_medium -> Scale(1.0/nevt);
	h_taus_selected_tight -> Scale(1.0/nevt);
	h_taus_gen -> Scale(1.0/nevt);
	h_taus_selected_loose_eta -> Scale(1.0/nevt);
	h_taus_selected_medium_eta -> Scale(1.0/nevt);
	h_taus_selected_tight_eta -> Scale(1.0/nevt);
	h_taus_gen_eta -> Scale(1.0/nevt);
	
	TFile tauroot("/uscms/home/ztao/work/CU_ttH_WD/Outputs/taueff.root", "RECREATE");
	h_taus_selected_loose -> Write();
	h_taus_selected_medium -> Write();
	h_taus_selected_tight -> Write();
	h_taus_gen -> Write();
	h_taus_selected_loose_eta -> Write();
	h_taus_selected_medium_eta -> Write();
	h_taus_selected_tight_eta -> Write();
	h_taus_gen_eta -> Write();

}
