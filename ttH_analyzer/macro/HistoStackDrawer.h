#ifndef HistoStackDrawer_H
#define HistoStackDrawer_H

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TList.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TCanvas.h"

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
		//h->SetFillStyle();
	}
	else if (channel=="TTW" or channel=="ttW") {
		h->SetFillColor(kGreen+2);
		//h->SetFillStyle();
	}
	else if (channel=="TTZ" or channel=="ttZ") {
		h->SetFillColor(kGreen-8);
		//h->SetFillStyle();
	}
	else if (channel=="EWK" or channel=="Electroweak") {
		h->SetFillColor(kViolet-2);
		//h->SetFillStyle();
	}
	else if (channel=="Rares") {
		h->SetFillColor(kAzure+1);
		//h->SetFillStyle();
	}
	else if (channel=="Fakes" or channel=="fakes_data") {
		h->SetFillColor(kBlack);
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

#endif
