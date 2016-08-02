#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

using namespace std;

void datacardsComparer(
					   TString rootfile0="ttH_2lss_1tau.input.root",
					   TString rootfile1="/Users/ztao/Documents/ttH/HIG15008_datacards/CERN/ttH_2lss_1tau.input.root",
					   TString rootfile2="/Users/ztao/Documents/ttH/HIG15008_datacards/Tallinn/ttH_2lss_1tau.input.root"
					   )
{
	// List of histograms to compare:
	vector<const TString> histList = {"ttH_htt", "ttH_hww", "ttH_hzz",
									  "TTW", "TTZ", "Rares",
									  "data_obs", "fakes_data", "flips_data"};
	
	TFile* f[3];
	f[0] = new TFile(rootfile0);
	f[1] = new TFile(rootfile1);
	f[2] = new TFile(rootfile2);

	for (auto & hname : histList) {
		TH1D* h[3];
		for (int i = 0; i < 3; i++) {
			h[i] = (TH1D*) f[i]->Get("x_"+hname);
		}

		// Plot
		TCanvas c;
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		gPad->SetLogy(1);

		h[0]->SetLineColor(4);
		h[0]->SetMarkerColor(4);
		h[0]->SetMarkerStyle(4);
		h[1]->SetLineColor(1);
		h[1]->SetMarkerColor(1);
		h[1]->SetMarkerStyle(4);
		h[2]->SetLineColor(2);
		h[2]->SetMarkerColor(2);
		h[2]->SetMarkerStyle(4);

		TLegend *l = new TLegend(0.6,0.7,0.8,0.8);
		l->AddEntry(h[0],"CU","p");
		l->AddEntry(h[1],"HIG-15-008","p");
		l->AddEntry(h[2],"Tallinn", "p");

		h[0]->SetMaximum(100);
		h[0]->SetMinimum(.0001);

		h[0]->Draw();
		h[1]->Draw("same");
		h[2]->Draw("same");
		l->Draw("same");

		c.SaveAs("compare_"+hname+".pdf");
		
	} // end of histList loop
	
}
