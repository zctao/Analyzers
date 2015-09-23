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

void GenInfoPlotter (const TString input_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output.root") {
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
	vector<float>* gen_top_mass;
	vector<int>* gen_xDaug_pdgId;
	vector<float>* gen_xDaug_pt;
	vector<float>* gen_xDaug_eta;
	vector<float>* gen_xDaug_phi;
	vector<float>* gen_xDaug_mass;

	TBranch* b_gen_top_mass;
	TBranch* b_gen_xDaug_pdgId;
	TBranch* b_gen_xDaug_pt;
	TBranch* b_gen_xDaug_eta;
	TBranch* b_gen_xDaug_phi;
	TBranch* b_gen_xDaug_mass;

	gen_top_mass = 0;
	gen_xDaug_pdgId = 0;
	gen_xDaug_pt = 0;
	gen_xDaug_eta = 0;
	gen_xDaug_phi = 0;
	gen_xDaug_mass = 0;
	
	tree->SetBranchAddress("gen_top_mass", &gen_top_mass, &b_gen_top_mass);
	tree->SetBranchAddress("gen_x_daughter_pdgId", &gen_xDaug_pdgId, &b_gen_xDaug_pdgId);
	tree->SetBranchAddress("gen_x_daughter_pt", &gen_xDaug_pt, &b_gen_xDaug_pt);
	tree->SetBranchAddress("gen_x_daughter_eta", &gen_xDaug_eta, &b_gen_xDaug_eta);
	tree->SetBranchAddress("gen_x_daughter_phi", &gen_xDaug_phi, &b_gen_xDaug_phi);
	tree->SetBranchAddress("gen_x_daughter_mass", &gen_xDaug_mass, &b_gen_xDaug_mass);

	// Histograms
	TH2F* h_mtautau_mtop1 = new TH2F("mTTmtop1","m(#tau#tau) vs m(top1)", 10, 124, 126, 50, 150, 200);	
	TH2F* h_mtautau_mtop2 = new TH2F("mTTmtop2","m(#tau#tau) vs m(top2)", 10, 124, 126, 50, 150, 200);
	TH2F* h_ptautau_mtop1 = new TH2F("pTTmtop1","|p(#tau#tau)| vs m(top1)", 120, 0, 600, 50, 150, 200);	
	TH2F* h_ptautau_mtop2 = new TH2F("pTTmtop2","|p(#tau#tau)| vs m(top2)", 120, 0, 600, 50, 150, 200);
	
	// ----------------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ----------------------------------------------------------------------------------------------------------------
	
	// event loop
	for (int ievt=0; ievt<nevt; ++ievt) {
		tree -> GetEntry(ievt);
		
		TLorentzVector tau1, tau2;
		tau1.SetPtEtaPhiM(gen_xDaug_pt->at(0), gen_xDaug_eta->at(0), gen_xDaug_phi->at(0), gen_xDaug_mass->at(0));
		tau2.SetPtEtaPhiM(gen_xDaug_pt->at(1), gen_xDaug_eta->at(1), gen_xDaug_phi->at(1), gen_xDaug_mass->at(1));

		float mtautau, ptautau, mtop1, mtop2;

		mtautau	= (tau1+tau2).M();
		ptautau = (tau1+tau2).P();

		if (gen_top_mass->size()!=2) {
			cout << "Err...We've got a problem." << endl;
			return;
		}

		if (gen_top_mass->at(0)>gen_top_mass->at(1)) {
			mtop1 = gen_top_mass->at(0);
			mtop2 = gen_top_mass->at(1);
		}
		else {
			mtop1 = gen_top_mass->at(1);
			mtop2 = gen_top_mass->at(0);
		}

		h_mtautau_mtop1 -> Fill(mtautau, mtop1);
		h_mtautau_mtop2 -> Fill(mtautau, mtop2);
		h_ptautau_mtop1 -> Fill(ptautau, mtop1);
		h_ptautau_mtop2 -> Fill(ptautau, mtop2);
		
	} // end of event loop

	TFile histfile("/uscms/home/ztao/work/CU_ttH_WD/Outputs/histograms.root", "RECREATE");
	h_mtautau_mtop1 -> Write();
	h_mtautau_mtop2 -> Write();
	h_ptautau_mtop1 -> Write();
	h_ptautau_mtop2 -> Write();
}
