#ifndef HistoStackDrawer_H
#define HistoStackDrawer_H

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TList.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TPad.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"

#include <iostream>
#include <vector>
#include <map>

void SetHistStyle(TH1D* h, TString channel)
{
	if (channel=="data_obs" or channel.Contains("Observed")) {
		h->SetMarkerStyle(20);
		h->SetMarkerColor(1);
	}
	else if (channel=="ttH") {
		h->SetFillColor(kRed);
		h->SetLineColor(kBlack);
		//h->SetFillStyle();
	}
	else if (channel=="TTW" or channel=="ttW") {
		h->SetFillColor(kGreen+2);
		h->SetLineColor(kBlack);
		//h->SetFillStyle();
	}
	else if (channel=="TTZ" or channel=="ttZ") {
		h->SetFillColor(kGreen-8);
		h->SetLineColor(kBlack);
		//h->SetFillStyle();
	}
	else if (channel=="EWK" or channel=="Electroweak") {
		h->SetFillColor(kViolet-2);
		h->SetLineColor(kBlack);
		//h->SetFillStyle();
	}
	else if (channel=="Rares") {
		h->SetFillColor(kAzure+1);
		h->SetLineColor(kBlack);
		//h->SetFillStyle();
	}
	else if (channel=="Fakes" or channel=="fakes_data") {
		h->SetFillColor(kBlack);
		h->SetLineColor(kBlack);
		h->SetFillStyle(3335);
	}
	else {
		std::cout << "WARNING: undefined channel name. No style applied."
				  << std::endl;
	}
}

void HistoStackDrawer(std::map<TString, TH1D*>);
// overload
void HistoStackDrawer(std::map<TString, std::vector<TH1D*>>);
void HistoStackDrawer(std::vector<TString>, std::vector<TH1D*>);
void HistoStackDrawer(TString);

std::map<TString, TString> properTitle = {
	{"tau_pt","#tau_{h} p_{T} [GeV]"},
	{"tau_eta","#tau_{h} #eta"},
	{"njet", "N(jet, p_{T}>25 GeV)"},
	{"mindr_lep1_jet","mindr_lep1_jet"},
	{"mindr_lep2_jet","mindr_lep2_jet"},
	{"avg_dr_jet","avg_dr_jet"},
	{"max_lep_eta","max(#eta_{l_{1}},#eta_{l_{2}})"},
	{"met","E^{miss}_{T} [GeV]"},
	{"mT_met_lep1","mT_met_lep1"},
	{"mht","MHT [GeV]"},
	{"dr_leps","dR(l_{1},l_{2})"},
	{"dr_lep1_tau","dR(l_{ldg},#tau)"},
	{"lep1_pt","Leading lepton p_{T} [GeV]"},
	{"lep1_eta","Leading lepton #eta"},
	{"lep2_pt","Subleading lepton p_{T} [GeV]"},
	{"lep2_eta", "Subleading lepton #eta"},
	{"mTauTauVis1","m^{vis}_{#tau,l_{1}}"},
	{"mTauTauVis2","m^{vis}_{#tau,l_{2}}"},
	{"sum_charges","Sum of charge"}
};

#endif
