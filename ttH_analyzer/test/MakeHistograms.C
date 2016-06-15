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

void MakeHistograms(
					vector<const TString> sigs = {},
					vector<const TString> bkgs = {},
					vector<const TString> syst = {}
)
{
	// List of histogram names
	vector<TString> name_sigs = {"ttH_htt"};
	vector<TString> name_bkgs = {};
	vector<TString> name_syst = {};
	
	// Open files
	vector<TFile*> f_sigs;
	for (auto sig : sigs) {
		TFile* f_sig = new TFile(sig);
		f_sigs.push_back(f_sig);
	}

	vector<TFile*> f_bkgs;
	for (auto bkg : bkgs) {
	    TFile* f_bkg = new TFile(bkg);
		f_bkgs.push_back(f_bkg);
	}

	vector<TFile*> f_syst;
	for (auto & s : syst) {
	    TFile* f_s = new TFile(s);
		f_syst.push_back(f_s);
	}

	// Get Histograms from root files
	vector<TH1D*> hist_sigs;
	vector<TH1D*> hist_bkgs;
	vector<TH1D*> hist_syst;

	auto it_name_sig = name_sigs.begin();
	for (auto f_sig : f_sigs) {
		TH1D* h = (TH1D*)f_sig->Get("ttHtaus/h_MVA_shape");
		TString hname = *it_name_sig;
		h->SetName("x_"+hname);
		hist_sigs.push_back(h);
		++it_name_sig;
	}

	auto it_name_bkg = name_bkgs.begin();
	for (auto f_bkg : f_bkgs) {
		TH1D* h = (TH1D*)f_bkg->Get("ttHtaus/h_MVA_shape");
		TString hname = *it_name_bkg;
		h->SetName("x_"+hname);
		hist_bkgs.push_back(h);
		++it_name_bkg;
	}

	auto it_name_syst = name_syst.begin();
	for (auto f_s : f_syst) {
		TH1D* h = (TH1D*)f_s->Get("ttHtaus/h_MVA_shape");
		TString hname = *it_name_syst;
		h->SetName("x_"+hname);
		hist_syst.push_back(h);
		++it_name_syst;
	}

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
	
	
	// Write to output root file
	TFile *output_combine = new TFile("ttH_2lss_1tau_shape.root", "RECREATE");
	//gDirectory->pwd();

	for (auto & h : hist_sigs) {
		h->Write();
	}
	
	for (auto & h : hist_bkgs) {
		h->Write();
	}

	for (auto & h : hist_syst) {
		h->Write();
	}
	
}
