#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

//TH1D* combineShapes(vector<TH1D*>);  // produce 1D shapes for combine
//void draw2DBDT(); // 
//void drawStackHists();  // draw stack histograms after event selection
//void drawCutflow();  // draw cutflow histogram taken from tree

void MakeHistograms(vector<TString> sample_names = {"ttH_htt"})
{

	TString fname_prefix = "testNtuple_";
	TString sys_coname = "_CMS_ttHl_";
	
	// List of btag systematics
	TString BTagSysts [16] =
		{"LFUp","LFDown","HFUp","HFDown",
		 "HFStats1Up","HFStats1Down","HFStats2Up","HFStats2Down",
		 "LFStats1Up","LFStats1Down","LFStats2Up","LFStats2Down",
		 "cErr1Up","cErr1Down","cErr2Up","cErr2Down"};

	// Histograms
	vector<TH1D*> hists;
	
	for (auto & sname : sample_names) {
		//Open files
		TFile* f = new TFile(fname_prefix+sname+".root");
		TFile* f_jesup = new TFile(fname_prefix+sname+"_JESUp.root");
		TFile* f_jesdo = new TFile(fname_prefix+sname+"_JESUp.root");
				
		TH1D* h_central = (TH1D*)f->Get("ttHtaus/h_MVA_shape");
		h_central->SetName("x_"+sname);
		hists.push_back(h_central);
				
		TH1D* h_jesup = (TH1D*)f_jesup->Get("ttHtaus/h_MVA_shape");
		h_jesup->SetName("x_" + sname + sys_coname + "JESUp");
		hists.push_back(h_jesup);

		TH1D* h_jesdo = (TH1D*)f_jesdo->Get("ttHtaus/h_MVA_shape");
		h_jesdo->SetName("x_" + sname + sys_coname + "JESDown");
		hists.push_back(h_jesdo);
	
		for (auto & bsys : BTagSysts) {
			TH1D* h_ = (TH1D*)f->Get("ttHtaus/h_MVA_shape_"+bsys);
			h_->SetName("x_" + sname + sys_coname + "_btag_" + bsys);
			hists.push_back(h_);
		}		
		
	}

/*
	// Get Ntuples
	vector<TTree*> woods_sig;
	for (auto f_sig : f_sigs) {
		TTree* tree_sig = (TTree*) f_sig->Get("ttHtaus/eventTree");
		woods_sig.push_back(tree_sig);
	}

	vector<TTree*> woods_bkg;
	for (auto f_bkg : f_bkgs) {
		TTree* tree_bkg = (TTree*) f_bkg->Get("ttHtaus/eventTree");
		woods_bkg.push_back(tree_bkg);
	}	
*/
	
	// Write to output root file
	TFile *output_combine = new TFile("ttH_2lss_1tau_shape.root", "RECREATE");
	//gDirectory->pwd();

	for (auto & h : hists) {
		h->Write();
	}
	
}
