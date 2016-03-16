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

void syncNtupleComparer(
						bool xsel = true,
						const TString inputFile1 = "/afs/cern.ch/work/z/ztao/private/ttH/CMSSW_7_6_3_patch2/src/Analyzers/ttH_analyzer/test/ttHtausNtuple.root", 
						const TString inputFile2 = "/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple.root",
						const TString inputFile3 = "/afs/cern.ch/user/k/kaehatah/public/ntuples/ttHJetToTT_M125_13TeV_ntuples_sync.root",
						const TString inputFile4 = "/afs/cern.ch/user/m/matze/public/ttH/sync_ntuple.root"
						//const TString inputFile1 = "ttHtausNtuple.root",
						//const TString inputFile2 = "~/Desktop/Scratch/syncNtuple.root",
						//const TString inputFile3 = "~/Desktop/Scratch/ttHJetToTT_M125_13TeV_ntuples_sync.root",
						//const TString inputFile4 = "~/Desktop/Scratch/sync_ntuple.root"
)
{	
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
	TTree* tree4 = (TTree*) f4->Get("sync");

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
	
	for (const auto &branch : *branches1) {

		TString presel_cut = "";
		TString suffix = "";
		
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

		if (tree4->GetBranch(bname) != nullptr and xsel) {
			
			l->AddEntry(tree4, "ND", "l");
			presel_cut = " && n_presel_jet >= 2 && (n_presel_mu + n_presel_ele) == 2 && n_presel_tau >= 1";
			suffix = "_xsel";
			
			if ( tree4->GetEntries(bname+">-233") != 0)
				forest.push_back(tree4);
		}
		
		for (auto const &tree : forest) {
			
			if (bname.Contains("_dxy") or bname.Contains("_dz")) {
				if (tree == *forest.begin()) {
					tree->Draw("abs("+bname+")", bname+">-233"+presel_cut);
					TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
					htemp->GetXaxis()->SetTitle("abs("+bname+")");
				}
				else
					tree->Draw("abs("+bname+")", bname+">-233"+presel_cut, "same");
			}
			else {
				if (tree == *forest.begin())
					tree->Draw(bname, bname+">-233"+presel_cut);
				else 
					tree->Draw(bname, bname+">-233"+presel_cut, "same");
			}
			
			gPad->Update();
		}

		l->Draw("same");

		//c.SaveAs("./syncPlots/"+bname+".pdf");
		c.SaveAs("~ztao/www/"+bname+suffix+".png");

		delete l;
	}

}
