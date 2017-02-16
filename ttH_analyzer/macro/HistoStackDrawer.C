#include "HistoStackDrawer.h"

void HistoStackDrawer(std::map<TString, TH1D*> histoMap)
{
	if (histoMap.size()==0) {
		std::cout << "Input is empty! Nothing's plotted!" << std::endl;
		return;
	}

	THStack* hs = new THStack("", "");
	TLegend* l = new TLegend(0.82,0.30,0.95,0.57);

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
		if (channel=="TTW")
			channel = "ttW";
		
		SetHistStyle(hm.second, channel);

		hs->Add(hm.second);
		l->AddEntry(hm.second, channel, "f");
	}
	
	// plotting
	//gROOT->LoadMacro("tdrstyle.C");
	setTDRStyle();

	// CMS_lumi
	//gROOT->LoadMacro("CMS_lumi.C");
	writeExtraText = true;
	extraText = "Preliminary";
	lumi_13TeV = "36.8 fb^{-1}";
	lumi_sqrtS = "13 TeV";

	int iPeriod = 4; // 4=13TeV, 0=free form (use lumi_sqrtS)
	int iPos = 0;
	if (iPos==0) relPosX = 0.12;

	// Canvas
	TCanvas* c = new TCanvas("c","",0,0,600,600);
	
	// for showing ratio and uncertainties
	float ymax = max(hs->GetMaximum(),h_obs->GetMaximum());

	// upper pad
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetTopMargin(0.08);
	pad1->SetBottomMargin(0.); // Upper and lower plot are joined
	//pad1->SetGridx();         // Vertical grid
	pad1->Draw();             // Draw the upper pad: pad1	
	pad1->cd();

	hs->SetMaximum(ymax*1.2);
	hs->Draw("HIST");
	l->Draw("same");
	h_obs->Draw("same");	
	
	c->cd();

	// lower pad
	TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->Draw();
	pad2->cd();

	// define the ratio plot
	TH1D* h_tot = (TH1D*)(hs->GetStack()->Last())->Clone();
	TH1D* h_ratio = (TH1D*)h_obs->Clone();
	TH1D* h_err = (TH1D*)h_tot->Clone();
	h_ratio->Divide(h_tot);
	h_err->Divide(h_tot);
	
	h_ratio->SetStats(0);
	h_ratio->Draw("ep");
	
	h_err->SetLineColor(1);
	h_err->SetFillColor(1);
	h_err->SetFillStyle(3003);
	h_err->Draw("same E2");

	float xmax = h_err->GetXaxis()->GetXmax();
	float xmin = h_err->GetXaxis()->GetXmin();
	TF1 *f1 = new TF1("f1","1",xmin,xmax);
	f1->SetLineColor(1);
	f1->SetLineWidth(1);
	f1->SetLineStyle(1);
	f1->Draw("same");

	// Upper pad Y axis
	pad1->cd();
	hs->GetYaxis()->SetTitle("Events");
	hs->GetYaxis()->SetTitleOffset(1.25);
	hs->GetYaxis()->SetTitleSize();
	hs->GetYaxis()->SetNdivisions(510);
	hs->GetYaxis()->SetLabelSize(0.035);
	hs->GetYaxis()->ChangeLabel(1,0,0);
	
	// Lower pad Y axis
	pad2->cd();
	h_ratio->GetYaxis()->SetTitle("Data/Pred.");
	h_ratio->GetYaxis()->SetTitleOffset(0.51);
	h_ratio->GetYaxis()->SetTitleSize(0.1);
	h_ratio->GetYaxis()->CenterTitle();
	h_ratio->GetYaxis()->SetNdivisions(505);
	h_ratio->GetYaxis()->SetRangeUser(0.0,2.0);
	h_ratio->GetYaxis()->SetLabelSize(0.10);
	h_ratio->GetYaxis()->ChangeLabel(1,0,0);
	h_ratio->GetYaxis()->ChangeLabel(-1,0,0);
	// Lower pad X axis
	TString pname = hs->GetName();
	h_ratio->GetXaxis()->SetTitle(properTitle[pname]);
	h_ratio->GetXaxis()->SetTitleSize(0.1);
	//h_ratio->GetXaxis()->SetTitleOffset(4.0);
	h_ratio->GetXaxis()->SetNdivisions(510);
	h_ratio->GetXaxis()->SetLabelSize(0.1);
	
	// writing the lumi information and the CMS "logo"
	CMS_lumi( c, iPeriod, iPos );

	c->Update();
	c->RedrawAxis();
	c->GetFrame()->Draw();
	
	c->SaveAs(pname+".pdf");
	delete c;
}

void HistoStackDrawer(TString filename="histograms.root")
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

		//std::cout << ch << std::endl;
		std::vector<TH1D*> hists;
		
		f->cd(ch);
		TIter subnext(gDirectory->GetListOfKeys());
		TKey *subkey;

		while ((subkey = (TKey*)subnext())) {
			TClass *c1 = gROOT->GetClass(subkey->GetClassName());
			if (!c1->InheritsFrom("TH1")) continue;

			TH1D* h = (TH1D*)subkey->ReadObj();
			hists.push_back(h);
			//std::cout << h->GetName() << std::endl;
		}

		histoCollection[ch] = hists;
	}

	HistoStackDrawer(histoCollection);	
	
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

