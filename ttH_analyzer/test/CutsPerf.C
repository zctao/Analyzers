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
	const float pT = 25.0;
	const float eta = 2.3;
}

void CutPerf (const TString sig_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_sig.root", 
	      const TString bkg_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_TTJets.root") 
{
	// Read ntuples
	TFile* f_sig = new TFile(sig_file);
	TFile* f_bkg = new TFile(bkg_file);
	
	TTree* tree = (TTree*) f_sig->Get("ttHsyncExercise/EventTree");
	TTree* tree_bkg = (TTree*) f_bkg->Get("ttHsyncExercise/EventTree");

	const int nevt = tree->GetEntries();
	const int nbkg = tree_bkg->GetEntries();

	if (nevt == 0 or nbkg == 0) {
		cout << "File doesn't exist or is empty, returning..." << endl;
		cout << endl;
		return;
	}
	else {
		cout << "Number of signal events :" << nevt << endl;
		cout << "Number of background events :" << nbkg << endl;
	}
	
	// Define leafs and branches

	std::vector<float>* tau_pt_noniso;
	std::vector<float>* tau_eta_noniso;
	std::vector<float>* tau_pt_loose;
	std::vector<float>* tau_eta_loose;
	std::vector<float>* tau_pt_medium;
	std::vector<float>* tau_eta_medium;
	std::vector<float>* tau_pt_tight;
	std::vector<float>* tau_eta_tight;

	std::vector<float>* gen_tau_pt;
	std::vector<float>* gen_tau_eta;
	std::vector<int>* gen_tauclass;

	// background tree
	std::vector<float>* tau_pt_noniso_bkg;
	std::vector<float>* tau_pt_loose_bkg;
	std::vector<float>* tau_pt_medium_bkg;
	std::vector<float>* tau_pt_tight_bkg;

	TBranch* b_tau_pt_noniso;
	TBranch* b_tau_eta_noniso;
	TBranch* b_tau_pt_loose;
	TBranch* b_tau_eta_loose;
	TBranch* b_tau_pt_medium;
	TBranch* b_tau_eta_medium;
	TBranch* b_tau_pt_tight;
	TBranch* b_tau_eta_tight;
	TBranch* b_gen_tau_pt;
	TBranch* b_gen_tau_eta;
	TBranch* b_gentauclass;

	TBranch* b_tau_pt_noniso_bkg;
	TBranch* b_tau_pt_loose_bkg;
	TBranch* b_tau_pt_medium_bkg;
	TBranch* b_tau_pt_tight_bkg;

	tau_pt_noniso = 0;
	tau_eta_noniso = 0;
	tau_pt_loose = 0;
	tau_eta_loose = 0;
	tau_pt_medium = 0;
	tau_eta_medium = 0;
	tau_pt_tight = 0;
	tau_eta_tight = 0;

	gen_tau_pt = 0;
	gen_tau_eta = 0;
	gen_tauclass = 0;

	tau_pt_noniso_bkg = 0;
	tau_pt_loose_bkg = 0;
	tau_pt_medium_bkg = 0;
	tau_pt_tight_bkg = 0;
	
	tree->SetBranchAddress("tau_pt_noniso", &tau_pt_noniso, &b_tau_pt_noniso);
	tree->SetBranchAddress("tau_eta_noniso", &tau_eta_noniso, &b_tau_eta_noniso);
	tree->SetBranchAddress("tau_pt_loose", &tau_pt_loose, &b_tau_pt_loose);
	tree->SetBranchAddress("tau_eta_loose", &tau_eta_loose, &b_tau_eta_loose);
	tree->SetBranchAddress("tau_pt_medium", &tau_pt_medium, &b_tau_pt_medium);
	tree->SetBranchAddress("tau_eta_medium", &tau_eta_medium, &b_tau_eta_medium);
	tree->SetBranchAddress("tau_pt_tight", &tau_pt_tight, &b_tau_pt_tight);
	tree->SetBranchAddress("tau_eta_tight", &tau_eta_tight, &b_tau_eta_tight);

	tree->SetBranchAddress("gen_x_daughter_pt", &gen_tau_pt, &b_gen_tau_pt);	
	tree->SetBranchAddress("gen_x_daughter_eta", &gen_tau_eta, &b_gen_tau_eta);
	tree->SetBranchAddress("gen_tau_class", &gen_tauclass, &b_gentauclass);

	tree_bkg->SetBranchAddress("tau_pt_noniso", &tau_pt_noniso_bkg, &b_tau_pt_noniso_bkg);
	tree_bkg->SetBranchAddress("tau_pt_loose", &tau_pt_loose_bkg, &b_tau_pt_loose_bkg);
	tree_bkg->SetBranchAddress("tau_pt_medium", &tau_pt_medium_bkg, &b_tau_pt_medium_bkg);
	tree_bkg->SetBranchAddress("tau_pt_tight", &tau_pt_tight_bkg, &b_tau_pt_tight_bkg);
	
	// Histograms
		
	TH1F* h_taus_gen = new TH1F("h_taus_gen","", 4, 0, 4);	
	TH1F* h_taus_selected_noniso = new TH1F("h_taus_selected_noniso", "", 4, 0, 4);
	TH1F* h_taus_selected_loose = new TH1F("h_taus_selected_loose", "", 4, 0, 4);
	TH1F* h_taus_selected_medium = new TH1F("h_taus_selected_medium", "", 4, 0, 4);
	TH1F* h_taus_selected_tight = new TH1F("h_taus_selected_tight", "", 4, 0, 4);

	TH1F* h_taus_gen_etacut = new TH1F("h_taus_gen_etacut","", 4, 0, 4);	
	TH1F* h_taus_selected_noniso_etacut = new TH1F("h_taus_selected_noniso_etacut", "", 4, 0, 4);
	TH1F* h_taus_selected_loose_etacut = new TH1F("h_taus_selected_loose_etacut", "", 4, 0, 4);
	TH1F* h_taus_selected_medium_etacut = new TH1F("h_taus_selected_medium_etacut", "", 4, 0, 4);
	TH1F* h_taus_selected_tight_etacut = new TH1F("h_taus_selected_tight_etacut", "", 4, 0, 4);

	TH1F* h_tauhpt_gen = new TH1F("h_tauhpt_gen","", 100, 0, 500);
	
	TH1F* h_taueff_noniso = new TH1F("h_taueff_noniso","",100, 0, 500);
	TH1F* h_taueff_loose = new TH1F("h_taueff_loose","",100, 0, 500);
	TH1F* h_taueff_medium = new TH1F("h_taueff_medium","",100, 0, 500);
	TH1F* h_taueff_tight = new TH1F("h_taueff_tight","",100, 0, 500);

	TH1F* h_bkgtaupt_noniso = new TH1F("h_bkgtaupt_noniso","",100, 0, 500);
	TH1F* h_bkgtaupt_loose = new TH1F("h_bkgtaupt_loose","",100, 0, 500);
	TH1F* h_bkgtaupt_medium = new TH1F("h_bkgtaupt_medium","",100, 0, 500);
	TH1F* h_bkgtaupt_tight = new TH1F("h_bkgtaupt_tight","",100, 0, 500);
	
	// ----------------------------------------------------------------------------------------------------------------
	//        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
	// ----------------------------------------------------------------------------------------------------------------
	
	// event loop
	for (int ievt=0; ievt<nevt; ++ievt) {
		tree -> GetEntry(ievt);
		tree_bkg -> GetEntry(ievt);
		
		int n_gen = 0;
		int n_gen_etacut = 0;
		int n_noniso = 0;
		int n_noniso_etacut = 0;
		int n_loose = 0;
		int n_loose_etacut = 0;
		int n_medium = 0;
		int n_medium_etacut = 0;
		int n_tight = 0;
		int n_tight_etacut = 0;
		
		for (size_t i = 0; i < gen_tau_pt->size(); ++i) {
			if (gen_tauclass->at(i) == 0)
				continue;

			h_tauhpt_gen->Fill(gen_tau_pt->at(i));
			
			if (gen_tau_pt->at(i) > Cuts::pT) {
				n_gen++;
				if (abs(gen_tau_eta->at(i)) < Cuts::eta)
					n_gen_etacut++;
			}
		}
		h_taus_gen -> Fill(n_gen);
		h_taus_gen_etacut -> Fill(n_gen_etacut);

		for (size_t i = 0; i< tau_pt_noniso->size(); ++i) {
			if (tau_pt_noniso->at(i) > Cuts::pT) {
				n_noniso++;
				
				h_taueff_noniso->Fill(tau_pt_noniso->at(i));
				
				if (abs(tau_eta_noniso->at(i)) < Cuts::eta)
					n_noniso_etacut++;
			}
		}
		h_taus_selected_noniso -> Fill(n_noniso);
		h_taus_selected_noniso_etacut -> Fill(n_noniso_etacut);
		
		for (size_t i = 0; i< tau_pt_loose->size(); ++i) {
			if (tau_pt_loose->at(i) > Cuts::pT) {
				n_loose++;

				h_taueff_loose->Fill(tau_pt_loose->at(i));
				
				if (abs(tau_eta_loose->at(i)) < Cuts::eta)
					n_loose_etacut++;
			}
		}
		h_taus_selected_loose -> Fill(n_loose);
		h_taus_selected_loose_etacut -> Fill(n_loose_etacut);

		for (size_t i = 0; i< tau_pt_medium->size(); ++i) {
			if (tau_pt_medium->at(i) > Cuts::pT) {
				n_medium++;

				h_taueff_medium->Fill(tau_pt_medium->at(i));
				
				if (abs(tau_eta_medium->at(i)) < Cuts::eta)
					n_medium_etacut++;
			}
		}
		h_taus_selected_medium -> Fill(n_medium);
		h_taus_selected_medium_etacut -> Fill(n_medium_etacut);

		for (size_t i = 0; i< tau_pt_tight->size(); ++i) {
			if (tau_pt_tight->at(i) > Cuts::pT) {
				n_tight++;

				h_taueff_tight->Fill(tau_pt_tight->at(i));
				
				if (abs(tau_eta_tight->at(i)) < Cuts::eta)
					n_tight_etacut++;
			}
		}
		h_taus_selected_tight -> Fill(n_tight);
		h_taus_selected_tight_etacut -> Fill(n_tight_etacut);

		// background
		for (size_t i = 0; i < tau_pt_noniso_bkg->size(); ++i) {
			if (tau_pt_noniso_bkg->at(i) > Cuts::pT) 
				h_bkgtaupt_noniso -> Fill(tau_pt_noniso_bkg->at(i));
		}
		
		for (size_t i = 0; i < tau_pt_loose_bkg->size(); ++i) {
			if (tau_pt_loose_bkg->at(i) > Cuts::pT) 
				h_bkgtaupt_loose -> Fill(tau_pt_loose_bkg->at(i));
		}

		for (size_t i = 0; i < tau_pt_medium_bkg->size(); ++i) {
			if (tau_pt_medium_bkg->at(i) > Cuts::pT) 
				h_bkgtaupt_medium -> Fill(tau_pt_medium_bkg->at(i));
		}

	    for (size_t i = 0; i < tau_pt_tight_bkg->size(); ++i) {
			if (tau_pt_tight_bkg->at(i) > Cuts::pT) 
				h_bkgtaupt_tight -> Fill(tau_pt_tight_bkg->at(i));
		}

	} // end of event loop
	
	h_taus_selected_noniso -> Sumw2();
	h_taus_selected_loose -> Sumw2();
	h_taus_selected_medium -> Sumw2();
	h_taus_selected_tight -> Sumw2();
	h_taus_gen -> Sumw2();
	h_taus_selected_noniso_etacut -> Sumw2();
	h_taus_selected_loose_etacut -> Sumw2();
	h_taus_selected_medium_etacut -> Sumw2();
	h_taus_selected_tight_etacut -> Sumw2();
	h_taus_gen_etacut -> Sumw2();
	
	h_taus_selected_noniso -> Scale(1.0/nevt);
	h_taus_selected_loose -> Scale(1.0/nevt);
	h_taus_selected_medium -> Scale(1.0/nevt);
	h_taus_selected_tight -> Scale(1.0/nevt);
	h_taus_gen -> Scale(1.0/nevt);
	h_taus_selected_noniso_etacut -> Scale(1.0/nevt);
	h_taus_selected_loose_etacut -> Scale(1.0/nevt);
	h_taus_selected_medium_etacut -> Scale(1.0/nevt);
	h_taus_selected_tight_etacut -> Scale(1.0/nevt);
	h_taus_gen_etacut -> Scale(1.0/nevt);

	// Efficiency
	h_tauhpt_gen -> Sumw2();
	h_taueff_noniso -> Sumw2();
	h_taueff_loose -> Sumw2();
	h_taueff_medium -> Sumw2();
	h_taueff_tight -> Sumw2();

	h_taueff_noniso -> Divide(h_tauhpt_gen);
	h_taueff_loose -> Divide(h_tauhpt_gen);
	h_taueff_medium -> Divide(h_tauhpt_gen);
	h_taueff_tight -> Divide(h_tauhpt_gen);
	
	TFile tauroot("/uscms/home/ztao/work/CU_ttH_WD/Outputs/taueff.root", "RECREATE");
	h_taus_selected_noniso -> Write();
	h_taus_selected_loose -> Write();
	h_taus_selected_medium -> Write();
	h_taus_selected_tight -> Write();
	h_taus_gen -> Write();
	h_taus_selected_noniso_etacut -> Write();
	h_taus_selected_loose_etacut -> Write();
	h_taus_selected_medium_etacut -> Write();
	h_taus_selected_tight_etacut -> Write();
	h_taus_gen_etacut -> Write();

	h_taueff_noniso->Write();
	h_taueff_loose->Write();
	h_taueff_medium->Write();
	h_taueff_tight->Write();

	h_bkgtaupt_noniso->Write();
	h_bkgtaupt_loose->Write();
	h_bkgtaupt_medium->Write();
	h_bkgtaupt_tight->Write();
	
}
