#include "HistoStackDrawer.h"

void HistoStackDrawer(std::map<TString, TH1D*> histoMap)
{
	if (histoMap.size()==0) {
		std::cout << "Input is empty! Nothing's plotted!" << std::endl;
		return;
	}

	TCanvas c;

	THStack* hs = new THStack("", "");
	TLegend* l = new TLegend(0.78,0.60,0.95,0.95);

	TH1D* h_obs = histoMap["data_obs"];
	SetHistStyle(h_obs, "data_obs");
	l->AddEntry(h_obs, "Observed", "p");
	
	bool first = true;
	for (auto & hm : histoMap) {
		
		TString channel = hm.first;	
		if (channel.Contains("data_obs")) continue;
		
		if (first) {
			// define THStack
			TString hname = hm.second->GetName();		
			hs->SetName(hname);
			
			first = false;
		}

		if (channel=="fakes_data")
			channel = "Fakes";
		if (channel=="TTZ")
			channel = "ttZ";
		if (channel=="TTZ")
			channel = "ttW";
		
		SetHistStyle(hm.second, channel);

		hs->Add(hm.second);
		l->AddEntry(hm.second, channel, "f");
	}

	// plotting
	TString pname = hs->GetName();

	hs->Draw("HIST");
	hs->GetXaxis()->SetTitle(pname);
	gPad->Modified();
	
	h_obs->Draw("same");
	l->Draw("same");

	c.SaveAs(pname+".pdf");
}

void HistoStackDrawer(std::vector<TString> channels, std::vector<TH1D*> hists)
{
	assert(channels.size()==hists.size());
	
	std::map<TString, TH1D*> histoMap;

	int i = 0;
	for (auto & ch : channels) {
		histoMap[ch] = hists.at(i++);
	}

	HistoStackDrawer(histoMap);
}


void HistoStackDrawer(std::map<TString, std::vector<TH1D*>> histoCollection)
{
	if (histoCollection.size()==0) {
		std::cout << "Input is empty!" << std::endl;
		return;
	}

	unsigned int hists_size = (histoCollection.begin()->second).size();

	for (unsigned int ih = 0; ih < hists_size; ++ih) {
		std::map<TString, TH1D*> histoMap;
		for (auto & hc : histoCollection) {
			histoMap[hc.first] = hc.second.at(ih);
		}

		HistoStackDrawer(histoMap);
	}
}

void HistoStackDrawer(TString filename)
{
	TFile* f = TFile::Open(filename);

	TIter next(f->GetListOfKeys());
	TKey *key;

	std::vector<TString> channels;
	
	while ((key = (TKey*)next())) {
		channels.push_back(key->GetTitle());
	}


	std::map<TString, std::vector<TH1D*>> histoCollection;

	for (auto & ch : channels) {

		std::cout << ch << std::endl;
		std::vector<TH1D*> hists;
		
		f->cd(ch);
		TIter subnext(gDirectory->GetListOfKeys());
		TKey *subkey;

		while ((subkey = (TKey*)subnext())) {
			TClass *c1 = gROOT->GetClass(subkey->GetClassName());
			if (!c1->InheritsFrom("TH1")) continue;

			TH1D* h = (TH1D*)subkey->ReadObj();
			hists.push_back(h);
			std::cout << h->GetName() << std::endl;
		}

		histoCollection[ch] = hists;
	}

	HistoStackDrawer(histoCollection);	
	
}

