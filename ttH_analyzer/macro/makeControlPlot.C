#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <map>

#include "Misc_Constants.h"
#include "Cross_Sections.h"
#include "TreeAnalyzer.h"

using namespace std;

const float LUMI = 12.9 * 1000; // 1/pb

void makeControlPlot(vector<TString> channels =
		{"ttH", "TTW", "TTZ", "EWK", "Rares", "fakes_data", "data_obs"})
{

	map<TString, vector<TH1D*>> histsCollection;
	
	for (auto & channel : channels) {

		vector<TH1D*> vhists;

		bool isdata = channel.Contains("data");
		vector<vector<int>> eventList;
		
		bool first = true;	
		auto samples = SamplesInChannel.at(channel);
		for (auto & sample : samples) {
			
			// open file
			TString fname;
			if (isdata)
				fname = dir_map.at(sample)+"output_"+channel+".root";
			else
				fname = dir_map.at(sample)+"output_"+sample+".root";
			
			TFile* f = TFile::Open(fname);
			// get tree
			TTree* tree = (TTree*) f->Get("ttHtaus/eventTree");
			
			vector<TH1D*> hists = TreeAnalyzer(tree, isdata, eventList);

			// Scale histograms here
			if (not isdata) {
				TH1D* h_SumGenWeight = (TH1D*)f->Get("ttHtaus/h_SumGenWeight");
				//TH1D* h_SumGenWeightxPU = (TH1D*)f->Get("ttHtaus/h_SumGenWeightxPU");
				float nSum = h_SumGenWeight->GetBinContent(1);
				//float nSum = h_SumGenWeightxPU->GetBinContent(1);
				float XS = xsection::xsection[string(sample)];

				for (auto & h : hists) {
					h->Sumw2();
					h->Scale(LUMI * XS / nSum);
				}
			}

			// Add histograms to the channel collection
			int ih = 0;
			for (auto & h : hists) {
				if (first)
					vhists.push_back(h);
				else {
					assert(vhists.at(ih)->GetName()==h->GetName());
					vhists.at(ih++)->Add(h);
				}
			}

			first = false;
		} // end of sample loop

		histsCollection[channel] = vhists;
	} // end of channel loop

	// write and plot histograms
	// output file
	TFile *out_hists = new TFile("histograms.root", "RECREATE");
	//out_hists.WriteObject(&histsCollection, "histsCollection");

	for (auto & channel : channels) {
		 TDirectory *dir = out_hists->mkdir(channel);
		 dir->cd();

		 for (auto & hist : histsCollection[channel])
			 hist->Write();
	}

	return;
}
