#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"

#include <vector>

template <typename T>
void drawHistograms(TString plotname,
					std::vector<T>& hists,
					std::vector<TString>& channels,
					std::vector<double> scales = {})
{

	int nhists = hists.size();
	
	assert(hists.size() == channels.size());
	
	if (scales.size()==0) { // default 1
		for (auto& ch : channels)
			scales.push_back(1.);
	}
	assert(scales.size() == hists.size());
	
	// Scale
	for (int ih = 0; ih < nhists; ++ih) {
		double norm = hists[ih]->Integral();
		if (norm > 0)
			hists[ih]->Scale(scales[ih]/norm);
	}
	
	TCanvas c;
	gStyle->SetOptStat(0);
	
	int iColor = 1;
	for (auto& h : hists) {
		h->SetLineColor(iColor++);
	}
	
	// Set axis maximum
	double ymax = 0;

	for (int ih = 0; ih < nhists; ++ih) {
		if (hists[ih]->GetMaximum()>ymax)
			ymax = hists[ih]->GetMaximum();
	}
	
	hists[0]->SetMaximum(ymax*1.2);
	hists[0]->GetXaxis()->SetTitle(plotname);
	hists[0]->Draw("E1");
	
	TLegend *l = new TLegend(0.6,0.7,0.8,0.8);
	for (int ih = 0; ih < nhists; ++ih) {
		l->AddEntry(hists[ih], channels[ih], "l");
		hists[ih]->Draw("same E1");
	}

	l->Draw("same");
	
	c.SaveAs(plotname+".pdf");

	delete l;
   
}

template void drawHistograms<TH1D*>(TString,
									std::vector<TH1D*>&,
									std::vector<TString>&,
									std::vector<double>);
