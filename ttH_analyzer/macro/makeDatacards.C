#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

const TString BTagSysts [16] =
		{"LFUp","LFDown","HFUp","HFDown",
		 "HFStats1Up","HFStats1Down","HFStats2Up","HFStats2Down",
		 "LFStats1Up","LFStats1Down","LFStats2Up","LFStats2Down",
		 "cErr1Up","cErr1Down","cErr2Up","cErr2Down"};

const TString fname_prefix = "output_";
const TString sys_coname = "_CMS_ttHl_";

// data directories
map<TString,TString> dir_map =
	{{"ttH_htt", "~/Documents/ttH/Outputs/"},
	 {"ttH_hww", "~/Documents/ttH/Outputs/"},
	 {"ttH_hzz", "~/Documents/ttH/Outputs/"},
	 {"TTW", "~/Documents/ttH/Outputs/"},
	 {"TTZ", "~/Documents/ttH/Outputs/"},
	 {"rares_TTTT", "~/Documents/ttH/Outputs/"},
	 {"rares_WZZ", "~/Documents/ttH/Outputs/"},
	 {"rares_tZq", "~/Documents/ttH/Outputs/"},
	 {"rares_WW", "~/Documents/ttH/Outputs/"},
	 {"data_obs", "~/Documents/ttH/Outputs/"},
	 {"fakes_data", "~/Documents/ttH/Outputs/"},
	 {"flips_data", "~/Documents/ttH/Outputs/"}
	};

const TString dataSamples [5] =
	{"DoubleMuon", "SingleMuon", "DoubleEG", "SingleElectron", "MuonEG"};

//const float lumi = (0.02 + 2.24) * 1000.; // 1/pb
const float lumi = 2.17 * 1000.; // 1/pb

map<TString, float> xs_map =
	{{"ttH_htt", 0.212}, {"ttH_hww", 0.212}, {"ttH_hzz", 0.212}, // ttHnonbb xs
	 {"TTW", 0.2043}, {"TTZ", 0.2529}, {"WZ", 4.102},
	 {"rares_TTTT", 0.009103}, {"rares_tZq", 0.0758}, {"rares_WW", 1.64},
	 {"rares_WZZ", 0.05565}};

vector<TH1D*> getShapesMC(const TString, const TString);
vector<TH1D*> getShapesData(const TString, const TString);
vector<TH1D*> getShapesData(const TString, vector<const TString>&);
void ScaleHist(TH1D*, int, float, float);
void ScaleHist(TH1D*, int, int);
void fillHistoFromTree(TH1D*, vector<vector<int>>&, TTree*);
int partition2DBDT(float, float);


void makeDatacards(
				   vector<TString> channels = {"ttH_htt", "ttH_hww", "ttH_hzz",
						   "TTW", "TTZ",
						   //"WZ",
						   "data_obs", "flips_data", "fakes_data",
						   "Rares"})
{
	vector<TH1D*> datacards;
	
	for (auto & channel : channels) {
		vector<TH1D*> shapes;
		
		if (channel.Contains("data")) {
			
			vector<const TString> data_dir;
			for (auto & ds : dataSamples) {
				data_dir.push_back(dir_map[channel] + ds +"/");
			}
			
		    //shapes = getShapesData(channel, dir_map[channel]);
			shapes = getShapesData(channel, data_dir);
		}
		else if (channel == "Rares") {
			//{"rares_TTTT", "rares_WZZ",/*"rares_WW",*/ "rares_tZq"};
		    vector<TH1D*> rare1 = getShapesMC("rares_TTTT",dir_map["rares_TTTT"]);
			vector<TH1D*> rare2 = getShapesMC("rares_WZZ", dir_map["rares_WZZ"]);
			vector<TH1D*> rare3 = getShapesMC("rares_tZq",dir_map["rares_tZq"]);
			vector<TH1D*> rare4 = getShapesMC("rares_WW", dir_map["rares_WW"]);

			int nhist = rare1.size();
			for (int i = 0; i < nhist; ++i) {
				TH1D* hrare_ = rare1.at(i);
				hrare_->Add(rare2.at(i));
				hrare_->Add(rare3.at(i));

				TString name = hrare_->GetName();
				name.ReplaceAll("rares_TTTT","Rares");
				hrare_->SetName(name);

				shapes.push_back(hrare_);
			}					
		}
		else {
			shapes = getShapesMC(channel, dir_map[channel]);
		}

	    datacards.insert(datacards.end(), shapes.begin(), shapes.end());
	
	}  // end of channel loop

	// Output file
	TFile *output_combine = new TFile("ttH_2lss_1tau.input.root", "RECREATE");
	
	for (auto & h : datacards) {
		h->Write();
	}
}

vector<TH1D*> getShapesMC(const TString channel, const TString dir)
{
	vector<TH1D*> hists;
	// Open files
	TFile* f = new TFile(dir + fname_prefix + channel+".root");
	if (not f->IsOpen()) {
		cout << "Warning: cannot find " << fname_prefix << channel
			 << ".root in " << dir << endl;
		cout << "return an empty vector." << endl;
		return hists;
	}

	// For scaling
	TH1I* h_nProcessed = (TH1I*)f->Get("ttHtaus/h_nProcessed");
	int nsamples = h_nProcessed->GetEntries();
	double xsection = xs_map[channel];
	
	
	TH1D* h_central = (TH1D*)f->Get("ttHtaus/h_MVA_shape");
	if (channel != "WZ") {
		ScaleHist(h_central, nsamples, lumi, xsection);
	}
	else {
		ScaleHist(h_central, nsamples, 10);
	}
	h_central->SetName("x_"+channel);
	hists.push_back(h_central);

	TH1D* h_jesup = (TH1D*)f->Get("ttHtaus/h_MVA_shape_JESUp");
	if (channel != "WZ") {
		ScaleHist(h_jesup, nsamples, lumi, xsection);
	}
	else {
		ScaleHist(h_jesup, nsamples, 10);
	}
	h_jesup->SetName("x_" + channel + sys_coname + "JESUp");
	hists.push_back(h_jesup);

	TH1D* h_jesdo = (TH1D*)f->Get("ttHtaus/h_MVA_shape_JESDown");
	if (channel != "WZ") {
		ScaleHist(h_jesdo, nsamples, lumi, xsection);
	}
	else {
		ScaleHist(h_jesdo, nsamples, 10);
	}
	h_jesdo->SetName("x_" + channel + sys_coname + "JESDown");
	hists.push_back(h_jesdo);

	for (auto & bsys : BTagSysts) {
		TH1D* h_ = (TH1D*)f->Get("ttHtaus/h_MVA_shape_"+bsys);
		ScaleHist(h_, nsamples, lumi, xsection);
		h_->SetName("x_" + channel + sys_coname + "_btag_" + bsys);
		hists.push_back(h_);
	}

	return hists;
}

vector<TH1D*> getShapesData(const TString channel, const TString dir)
{
	vector<TH1D*> hists;
	TFile* f_data = new TFile(dir + fname_prefix + channel + ".root");

	if (not f_data->IsOpen()) {
		cout << "Warning: cannot find " << fname_prefix << channel
			 << ".root in " << dir << endl;
		cout << "return an empty vector." << endl;
		return hists;
	}
	
	TH1D* h_ = (TH1D*)f_data->Get("ttHtaus/h_MVA_shape");
	h_->SetName("x_"+channel);
	hists.push_back(h_);

	return hists;
}

vector<TH1D*> getShapesData(const TString channel, vector<const TString>& dirs)
{
	vector<TH1D*> hists;

	vector<vector<int>> eventList;  // A set of vectors: (run, lumisection, event#)

	TH1D* h_ = new TH1D("h_MVA_shape","", 6, 0.5, 6.5);
	h_ -> Sumw2();
	//TH1D* h_sys = ...
	
	// Open files
	for (auto & dir : dirs) {
		TFile* f_data = new TFile(dir + fname_prefix + channel + ".root");

		TTree* tree = (TTree*) f_data->Get("ttHtaus/eventTree");

	    fillHistoFromTree(h_, eventList, tree);
	}

	h_->SetName("x_"+channel);
	hists.push_back(h_);
	
	return hists;
}

void ScaleHist(TH1D* h, int nsamples, float lumi, float xsection)
{
	//h -> Sumw2();
	h -> Scale(lumi * xsection / nsamples);
}

void ScaleHist(TH1D* h, int nsamples, int nevents) // WZ
{
	//h -> Sumw2();
	h -> Scale( nevents * 0.0371 / nsamples);  // FIX NEEED HERE
}

void fillHistoFromTree(TH1D* h, vector<vector<int>>& eventList, TTree* tree)
// Fill shape histogram from TTree entries.
//Duplicate events characterized by (run, lumisection, event#) are only filled once
{
	
	int run;
	int ls;
	int event;
	double weight;
	double mva_ttbar;
	double mva_ttV;

	tree->SetBranchAddress("nEvent", &event);
	tree->SetBranchAddress("ls", &ls);
	tree->SetBranchAddress("run", &run);
	tree->SetBranchAddress("MVA_2lss_ttbar", &mva_ttbar);
	tree->SetBranchAddress("MVA_2lss_ttV", &mva_ttV);
	tree->SetBranchAddress("evtWeight", &weight);

	int nEntries = tree->GetEntries();

	//int newListSize = eventList.size()+nEntries;
	
	for (int i = 0; i < nEntries; ++i) {
		tree->GetEntry(i);

		vector<int> eventid = {run, ls, event};
		if (run < 0 or ls < 0 or event < 0) continue;

		bool alreadyIncluded =
			find(eventList.begin(), eventList.end(), eventid) != eventList.end();
		// better algorithm?
		
		if (not alreadyIncluded) {
			
			eventList.push_back(eventid);
			//cout << "run, lumisection, event : " << eventid.at(2) << " "
			//	 << eventid.at(1) << " " << eventid.at(0) << endl;

			//fill histogram
			int bin = partition2DBDT(mva_ttbar, mva_ttV);
			h->Fill(bin, weight);
		}
	}

}

int partition2DBDT(float ttbar, float ttV)
/*
             bin 1       bin 2       bin 3       bin 4      bin 5       bin 6
2lss(ttbar) (-1.0,-0.2] (-1.0,-0.2] (-0.2,0.3]  (-0.2,0.3] (0.3,1.0]   (0.3,1.0]
2lss(ttV)   (-1.0,-0.1] (-0.1,1.0]  (-1.0,-0.1] (-0.1,1.0] (-1.0,-0.1] (-0.1,1.0]
 */
{
	int x = -99;
	int y = -99;
	
	if (ttbar > -1.0 and ttbar <= -0.2)
		x = 1;
	else if (ttbar > -0.2 and ttbar <= 0.3)
		x = 2;
	else if (ttbar > 0.3 and ttbar <= 1.0)
		x = 3;

	if (ttV > -1.0  and ttV <= -0.1)
		y = 1;
	else if (ttV > -0.1 and ttV <= 1.0)
		y = 2;
	
	return 2*(x-1)+y;
}
