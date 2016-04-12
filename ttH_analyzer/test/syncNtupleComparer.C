#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TH1.h"
#include "TAxis.h"

#include <iostream>
#include <vector>
#include <string>

void syncNtupleComparer(TString type)
{
	TString inputFile1, inputFile2, inputFile3, inputFile4;
	if (type == "signal") {
	    inputFile1 = "/afs/cern.ch/work/z/ztao/public/ttHTT_syncNtuple/ttHtausNtuple.root";
		inputFile2 = "/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple_ttH.root";
		inputFile3 = "/afs/cern.ch/user/k/kaehatah/public/ntuples/ttHJetToTT_M125_13TeV_ntuples_sync.root";	
	}
	else if (type == "ttbar") {
		inputFile1 = "/afs/cern.ch/work/z/ztao/public/ttHTT_syncNtuple/ttJetsNtuple.root";
		inputFile2 = "/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple_ttjets.root";
		inputFile3 = "/afs/cern.ch/user/k/kaehatah/public/ntuples/ttJet_13TeV_ntuples_sync.root";
	}
	inputFile4 = "/afs/cern.ch/user/m/matze/public/ttH/sync_ntuple.root";

	TFile* f1 = new TFile(inputFile1);
	TFile* f2 = new TFile(inputFile2);
	TFile* f3 = new TFile(inputFile3);
	TFile* f4 = new TFile(inputFile4);
	
	if (!(f1->IsOpen() and f2->IsOpen() and f3->IsOpen() and f4->IsOpen() )) {
		std::cout << "Cannot open the file ... " << std::endl;
		return;
	}

	TTree* tree1 = (TTree*) f1->Get("ttHtaus/eventTree");
	TTree* tree2 = (TTree*) f2->Get("syncTree");
	TTree* tree3 = (TTree*) f3->Get("tree");
	TTree* tree4;
	if (type == "signal")
		tree4 = (TTree*) f4->Get("ttH");
	if (type == "ttbar")
		tree4 = (TTree*) f4->Get("ttjets");

	tree1->SetFillColor(5);
	tree1->SetLineColor(0);
	tree2->SetLineColor(4);
	tree2->SetLineWidth(2);
	tree3->SetLineColor(2);
	tree4->SetLineColor(8);

	TCanvas c;
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TObjArray *branches1 = tree1->GetListOfBranches ();

	//std::vector<string> gname;
	std::vector<int> nevt_mu;
	std::vector<int> nevt_ele;
	std::vector<int> nevt_tau;
	std::vector<int> nevt_jet;

	for (const auto &branch : *branches1) {
		
		TString bname = branch->GetName();

		TLegend *l = new TLegend(0.86,0.35,0.98,0.47);

		vector<TTree*> forest;

		if (tree1->GetBranch(bname) != nullptr) {
			
			l->AddEntry(tree1, "Cornell","f");
			
			if ( tree1->GetEntries(bname+">-233") != 0)
				forest.push_back(tree1);
		}

		if (tree2->GetBranch(bname) != nullptr) {
			
			l->AddEntry(tree2, "LLR", "l");
			
			if ( tree2->GetEntries(bname+">-233") != 0)
				forest.push_back(tree2);
		}

		if (tree3->GetBranch(bname) != nullptr) {
			
			l->AddEntry(tree3, "Tallinn", "l");
			
			if ( tree3->GetEntries(bname+">-233") != 0)
				forest.push_back(tree3);
		}

		if (tree4->GetBranch(bname) != nullptr) {
			
			l->AddEntry(tree4, "ND", "l");
			
			if ( tree4->GetEntries(bname+">-233") != 0)
				forest.push_back(tree4);
		}

		for (const auto &tree : forest) {
			
			if (tree == *forest.begin())
				tree->Draw(bname, bname+">-233");
			else
				tree->Draw(bname, bname+">-233", "same");
			
			gPad->Update();

			if (bname.EqualTo("n_presel_mu"))
				nevt_mu.push_back(tree->GetEntries(bname+">0"));
			if (bname.EqualTo("n_presel_ele"))
				nevt_ele.push_back(tree->GetEntries(bname+">0"));
			if (bname.EqualTo("n_presel_tau"))
				nevt_tau.push_back(tree->GetEntries(bname+">0"));
			if (bname.EqualTo("n_presel_jet"))
				nevt_jet.push_back(tree->GetEntries(bname+">0"));
		}

		l->Draw("same");

		//c.SaveAs("./syncPlots/"+bname+".pdf");
		c.SaveAs("~ztao/www/"+bname+type+".png");

		delete l;
	}

	std::cout << "muon: ";
	for (auto n: nevt_mu)
		std::cout << n << "\t";
	std::cout << std::endl;
	
	std::cout << "ele: ";
	for (auto n: nevt_ele)
		std::cout << n << "\t";
	std::cout << std::endl;
	
	std::cout << "tau: ";
	for (auto n: nevt_tau)
		std::cout << n << "\t";
	std::cout << std::endl;
	
	std::cout << "jet: ";
	for (auto n: nevt_jet)
		std::cout << n << "\t";
	std::cout << std::endl;

}
