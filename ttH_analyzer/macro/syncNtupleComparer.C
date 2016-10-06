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

void syncNtupleComparer(TString inputFile1="", TString inputFile2="",
						TString inputFile3="", TString inputFile4="")
{
	/*
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
	*/
	
	TFile* f1 = new TFile(inputFile1);
	TFile* f2 = new TFile(inputFile2);
	TFile* f3 = new TFile(inputFile3);
	TFile* f4 = new TFile(inputFile4);

	vector<TTree*> trees;
	vector<TString> names;
	vector<TString> options;
	trees.reserve(4);
	names.reserve(4);
	options.reserve(4);
	
	if (f1->IsOpen()) {
		TTree* tree1 = (TTree*) f1->Get("syncTree");
		tree1->SetFillColor(5);
		tree1->SetLineColor(0);
		
		trees.push_back(tree1);
		names.push_back("Cornell");
		options.push_back("f");
	}
	else
		cout << "Cannot open " << inputFile1 << endl;

	if (f2->IsOpen()) {
		TTree* tree2 = (TTree*) f2->Get("syncTree");
		tree2->SetLineColor(4);
		tree2->SetLineWidth(2);

		trees.push_back(tree2);
		names.push_back("LLR");
		options.push_back("l");
	}
	else
		cout << "Cannot open " << inputFile2 << endl;

	if (f3->IsOpen()) {
		TTree* tree3 = (TTree*) f3->Get("");
		tree3->SetLineColor(2);
		//tree3->SetLineWidth(2);

		trees.push_back(tree3);
		names.push_back("Tallinn");
		options.push_back("l");
	}
	else
		cout << "Cannot open " << inputFile3 << endl;

	if (f4->IsOpen()) {
		TTree* tree4 = (TTree*) f4->Get("");
		tree4->SetLineColor(8);
		//tree4->SetLineWidth(2);

		trees.push_back(tree4);
		names.push_back("ND");
		options.push_back("l");
	}
	else
		cout << "Cannot open " << inputFile4 << endl;
	
	TCanvas c;
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	auto branches1 = trees[0]->GetListOfBranches ();

	//std::vector<string> gname;
	std::vector<int> nevt_mu;
	std::vector<int> nevt_ele;
	std::vector<int> nevt_tau;
	std::vector<int> nevt_jet;

	for (const auto branch : *branches1) {
		
		TString bname = branch->GetName();

		TLegend *l = new TLegend(0.86,0.35,0.98,0.47);

		for (unsigned int it = 0; it < trees.size(); ++it) {
			if (trees[it]->GetBranch(bname) != nullptr) {
				l->AddEntry(trees[it], names[it], options[it]);
			}

			if (trees[it]->GetEntries(bname+">-233") != 0) {
				if (it == 0)
					trees[it]->Draw(bname, bname+">-233");
				else
					trees[it]->Draw(bname, bname+">-233", "same");

				gPad->Update();

				if (bname.EqualTo("n_presel_mu"))
					nevt_mu.push_back(trees[it]->GetEntries(bname+">0"));
				if (bname.EqualTo("n_presel_ele"))
					nevt_ele.push_back(trees[it]->GetEntries(bname+">0"));
				if (bname.EqualTo("n_presel_tau"))
					nevt_tau.push_back(trees[it]->GetEntries(bname+">0"));
				if (bname.EqualTo("n_presel_jet"))
					nevt_jet.push_back(trees[it]->GetEntries(bname+">0"));
			}
		}

		l->Draw("same");

		c.SaveAs("~/Documents/ttH/syncPlots/"+bname+".png");
		//c.SaveAs("~ztao/www/"+bname+"_"+type+".png");

		delete l;
	}

	std::cout << "\t";
	for (auto name : names)
		std::cout << name << "\t";
	std::cout << std::endl;
	
	std::cout << "muon:";
	for (auto n: nevt_mu)
		std::cout << "\t" << n;
	std::cout << std::endl;
	
	std::cout << "ele:";
	for (auto n: nevt_ele)
		std::cout << "\t" << n;
	std::cout << std::endl;
	
	std::cout << "tau:";
	for (auto n: nevt_tau)
		std::cout << "\t" << n;
	std::cout << std::endl;
	
	std::cout << "jet:";
	for (auto n: nevt_jet)
		std::cout << "\t" << n;
	std::cout << std::endl;

}
