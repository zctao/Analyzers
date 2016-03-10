#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TObjArray.h"

#include <iostream>

void syncNtupleComparer(//const TString inputFile1 = "/afs/cern.ch/work/z/ztao/private/ttH/CMSSW_7_6_3_patch2/src/Analyzers/ttH_analyzer/test/ttHtausNtuple.root", 
						//const TString inputFile2 = "/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple.root"
						const TString inputFile1 = "ttHtausNtuple.root",
						const TString inputFile2 = "~/Desktop/Scratch/syncNtuple.root"
)
{	
	TFile* f1 = new TFile(inputFile1);
	TFile* f2 = new TFile(inputFile2);
	
	if (!(f1->IsOpen() and f2->IsOpen())) {
		std::cout << "Cannot open the file ... " << std::endl;
		return;
	}
	
	TTree* tree1 = (TTree*) f1->Get("ttHtaus/eventTree");
	TTree* tree2 = (TTree*) f2->Get("syncTree");

	tree1->SetFillColor(5);
	tree1->SetLineColor(0);

	TCanvas c;
	gStyle->SetOptStat(0);
	
	TObjArray *branches1 = tree1->GetListOfBranches ();
	for (const auto &branch : *branches1) {
		TString bname = branch->GetName();
		tree1->Draw(bname, bname+">-233");
		gPad->Update();
		tree2->Draw(bname, bname+">-233", "same");
		c->SaveAs("./syncPlots/"+bname+".pdf");
	}
}
