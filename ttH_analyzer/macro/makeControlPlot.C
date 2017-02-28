#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <map>

#include "Misc_Constants.h"
#include "Cross_Sections.h"
//#include "TreeAnalyzer.h"
#include "Analyzers/ttH_analyzer/interface/Types_enum.h"
#include "Analyzers/ttH_analyzer/interface/TreeAnalyzer.h"

//const float LUMI = 12.9 * 1000; // 1/pb
const float LUMI = 36.773 * 1000; // 1/pb

void makeControlPlot(vector<TString> channels =
		{"ttH", "TTW", "TTZ", "EWK", "Rares", "fakes_data", "data_obs"},
					 bool verbose=false)
{
	using namespace std;

	gROOT->ProcessLine(".L ../src/SFHelper.cc+");
	gROOT->ProcessLine(".L ../src/TreeAnalyzer.cc+");
	
	map<TString, vector<TH1D*>> histsCollection;
	
	for (auto & channel : channels) {

		vector<TH1D*> vhists;

		bool isdata = channel.Contains("data");
		vector<vector<unsigned long long>> eventList;
		
		// analysis and selection type
		Analysis_types AnaType = Analysis_types::Analyze_2lss1tau;
		Selection_types SelType = Selection_types::Signal_2lss1tau;
		if (channel.Contains("fakes"))
			SelType = Selection_types::Control_1lfakeable;
		if (channel.Contains("flips"))
			SelType = Selection_types::Control_2los1tau;
		
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
			if (not f->IsOpen()) {
				cout << "CANNOT open file " << fname << endl;
				continue;
			}
			
			// get tree
			TTree* tree = (TTree*) f->Get("ttHtaus/eventTree");

			TreeAnalyzer tana(tree, AnaType, SelType, isdata, verbose);		
			vector<TH1D*> hists =
				tana.makeHistograms(true, eventList); // 'true' for control region
				//TreeAnalyzer(tree, isdata, eventList, AnaType, SelType);

			// Scale histograms here
			if (not isdata) {
				TH1D* h_SumGenWeight = (TH1D*)f->Get("ttHtaus/h_SumGenWeight");
				//TH1D* h_SumGenWeightxPU = (TH1D*)f->Get("ttHtaus/h_SumGenWeightxPU");
				float nSum = h_SumGenWeight->GetBinContent(1);
				//float nSum = h_SumGenWeightxPU->GetBinContent(1);
				float XS = xsection::xsection[string(sample)];

				for (auto & h : hists) {
					h->Scale(LUMI * XS / nSum);
				}
			}

			// print out sample yields and number of entries
			assert(hists.size()>0 and string(hists.at(0)->GetName())==string("tau_pt"));
			int nentries = hists.at(0)->GetEntries();
			int nbins = hists.at(0)->GetNbinsX();
			double error = 0.;
			double yields = hists.at(0)->IntegralAndError(1,nbins,error);
			assert(yields == hists.at(0)->Integral());
			
			std::cout << sample << "  yields: " << yields << " +/- " << error
					  << "  nEntries: " << nentries << std::endl;
			
			// Add histograms to the channel collection
			int ih = 0;
			for (auto & h : hists) {
				if (first)
					vhists.push_back(h);
				else {
					assert(string(vhists.at(ih)->GetName())==string(h->GetName()));
					vhists.at(ih++)->Add(h);
				}
			}

			first = false;
		} // end of sample loop

		histsCollection[channel] = vhists;

		// print out channel yields and number of entries
		assert(vhists.size()>0 and string(vhists.at(0)->GetName())==string("tau_pt"));
		int chNentries = vhists.at(0)->GetEntries();
		int chNbins = vhists.at(0)->GetNbinsX();
		double chError = 0.;
		double chYields = vhists.at(0)->IntegralAndError(1,chNbins,chError);
		std::cout << "- - - - - - - - - - - - - - - - - - - - " << std::endl;
		std::cout << channel << "  yields:" << chYields << " +/- " << chError
				  << "  nEntries: " << chNentries << std::endl;
		std::cout << "----------------------------------------" << std::endl;
		
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
