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

#include "Cross_Sections.h"

using namespace std;

const TString BTagSysts [16] =
		{"_LFUp","_LFDown","_HFUp","_HFDown",
		 "_HFStats1Up","_HFStats1Down","_HFStats2Up","_HFStats2Down",
		 "_LFStats1Up","_LFStats1Down","_LFStats2Up","_LFStats2Down",
		 "_cErr1Up","_cErr1Down","_cErr2Up","_cErr2Down"};
const TString ThuSysts [4] = {"_xUp","_xDown","_yUp","_yDown"};

const TString fname_prefix = "output_";
const TString sys_coname = "_";//"_CMS_ttHl_";

// data directories
map<TString,TString> dir_map =
	{{"ttH_htt", "/eos/uscms/store/user/ztao/ttH_80X/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix/"},
	 {"ttH_hww", "/eos/uscms/store/user/ztao/ttH_80X/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix/"},
	 {"ttH_hzz", "/eos/uscms/store/user/ztao/ttH_80X/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix/"},
	 {"TTW", "/eos/uscms/store/user/ztao/ttH_80X/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/"},
	 {"TTZ", "/eos/uscms/store/user/ztao/ttH_80X/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/"},
	 //FIXME
	 //{"rares_TTTT", "~/Documents/ttH/Outputs/"},
	 //{"rares_WZZ", "~/Documents/ttH/Outputs/"},
	 //{"rares_tZq", "~/Documents/ttH/Outputs/"},
	 //{"rares_WW", "~/Documents/ttH/Outputs/"},
	 {"data_obs", "/eos/uscms/store/user/ztao/ttH_80X/"},
	 {"fakes_data", "/eos/uscms/store/user/ztao/ttH_80X/"},
	 {"flips_data", "/eos/uscms/store/user/ztao/ttH_80X/"}
	};

const TString DataSamples [5] =
	{"DoubleMuon", "SingleMuon", "DoubleEG", "SingleElectron", "MuonEG"};

const TString lepCategory [3] = {"_mumu", "_ee", "_emu"};
const TString bCategory[2] = {"_bloose", "_bmedium"};

const float LUMI = 12.9 * 1000.; // 1/pb

vector<TH1D*> getShapesMC(const TString, bool);
vector<TH1D*> getShapesData(const TString);
void fillHistoFromTree(TH1D*, vector<vector<int>>&, TTree*, vector<float>&);
int partition2DBDT(float, float);
TH1D* combineHistCategories(TFile*, TString, float, float, float, vector<float>&);
TH1D* combineHistCategories(TFile*, TString, float, float, float);

void makeDatacards(
				   vector<TString> channels = {"ttH_htt", "ttH_hww", "ttH_hzz",
						   "TTW", "TTZ",
						   //"WZ",
						   "data_obs", "flips_data", "fakes_data"
						   //"Rares"
						   },
				   bool doSys = false)
{
	vector<TH1D*> datacards;
	
	for (auto & channel : channels) {
		
		cout << "processing channel: " << channel << endl;
		
		vector<TH1D*> shapes;
		
		if (channel.Contains("data")) {
			shapes = getShapesData(channel);
		}
		else if (channel == "Rares") {
			//{"rares_TTTT", "rares_WZZ",/*"rares_WW",*/ "rares_tZq"};
		    vector<TH1D*> rare1 = getShapesMC("rares_TTTT", doSys);
			vector<TH1D*> rare2 = getShapesMC("rares_WZZ", doSys);
			vector<TH1D*> rare3 = getShapesMC("rares_tZq", doSys);
			vector<TH1D*> rare4 = getShapesMC("rares_WW", doSys);

			int nhist = rare1.size();
			for (int i = 0; i < nhist; ++i) {
				TH1D* hrare_ = rare1.at(i);
				hrare_->Add(rare2.at(i));
				hrare_->Add(rare3.at(i));
				hrare_->Add(rare4.at(i));

				TString name = hrare_->GetName();
				name.ReplaceAll("rares_TTTT","Rares");
				hrare_->SetName(name);

				shapes.push_back(hrare_);
			}					
		}
		else {
			shapes = getShapesMC(channel, doSys);
		}

	    datacards.insert(datacards.end(), shapes.begin(), shapes.end());
	
	}  // end of channel loop

	// Output file
	TFile *output_combine = new TFile("ttH_2lss_1tau.input.root", "RECREATE");
	
	for (auto & h : datacards) {
		h->Write();
	}
}

vector<TH1D*> getShapesMC(const TString channel, bool addSysts)
{
	vector<TH1D*> hists;
	
	// Open files
	const TString dir = dir_map[channel];
	const TString fname = fname_prefix + channel + ".root";
	TFile* f = new TFile(dir + fname);
	if (not f->IsOpen()) {
		cout << "Warning: cannot find " << fname  << " in " << dir << endl;
		cout << "return an empty vector" << endl;
		return hists;
	}

	// For scaling
	//TH1I* h_nProcessed = (TH1I*)f->Get("ttHtaus/h_nProcessed");
	//int nsample = h_nProcessed->GetEntries();
	TH1F* h_SumMCWeight = (TH1F*)f->Get("ttHtaus/h_SumGenWeight");
	float nSum = h_SumMCWeight->GetBinContent(1);
	
	float xs = 0.;
	if (channel.Contains("ttH"))
		xs = xsection::xsection["ttH"];
	else
		xs = xsection::xsection[string(channel)];

	vector<float> yields;
	yields.clear();
	yields.reserve(3 * 2);

	TH1D* h_central = combineHistCategories(f, "", nSum, LUMI, xs, yields);
	h_central->SetName("x_"+channel);
	hists.push_back(h_central);

	// print out yields in each category
	int icat = 0;
	float ysum = 0.;
	cout << channel << " : ";
	for (const TString lcat : lepCategory) {
		for (const TString bcat : bCategory) {
			cout << yields.at(icat) << "("<<lcat<<","<<bcat<<") ";
			ysum += yields.at(icat++);
		}
	}
	cout << ysum << "(sum)" << endl;

	if (not addSysts) return hists;
	
	// open jesup and jesdown files
	const TString fname_jesup = fname_prefix + channel + "_jesup.root";
	const TString fname_jesdo = fname_prefix + channel + "_jesdown.root";
	
	TFile* f_jesup = new TFile(dir + fname_jesup);
	if (not f_jesup -> IsOpen()) {
		cout << "Warning: cannot find " << fname_jesup << " in " << dir << endl;
	}
	
	TFile* f_jesdo = new TFile(dir + fname_jesdo);
	if (not f_jesdo -> IsOpen()) {
		cout << "Warning: cannot find " << fname_jesdo << " in " << dir << endl;
	}
	
	TH1D* h_jesup = combineHistCategories(f_jesup, "", nSum, LUMI, xs);
	TH1D* h_jesdo = combineHistCategories(f_jesdo, "", nSum, LUMI, xs);
	
	h_jesup->SetName("x_" + channel + sys_coname + "JESUp");
	hists.push_back(h_jesup);
	
	h_jesdo->SetName("x_" + channel + sys_coname + "JESDown");
	hists.push_back(h_jesdo);
	
	// btag 
	for (auto & bsys : BTagSysts) {
		TH1D* h_ = combineHistCategories(f, bsys, nSum, LUMI, xs);
		h_->SetName("x_" + channel + sys_coname + "btag" + bsys);
		hists.push_back(h_);
	}
	
	// theoretical
	for (auto & thu : ThuSysts) {
		TH1D* h_ = combineHistCategories(f, thu, nSum, LUMI, xs);
		h_->SetName("x_" + channel + sys_coname + "thu_shape_" + thu);
		hists.push_back(h_);
	}

	return hists;
}

vector<TH1D*> getShapesData(const TString channel)
{
	vector<TH1D*> hists;

	vector<vector<int>> eventList; // A set of vectors: (run, lumisection, event#)
	eventList.clear();

	TH1D* h_ = new TH1D("h_MVA_shape","", 7, 0.5, 7.5);
	h_ -> Sumw2();
	//TH1D* h_sys = ...

	vector<float> yields;
	yields.resize(3*2);
	// initialize yields
	for (size_t iy = 0; iy < yields.size(); ++iy) {
		yields.at(iy) = 0.;
	}
	
	// Open files
	for (auto & ds : DataSamples) {
		
		const TString dir = dir_map[channel] + ds + "/";
		const TString fname = fname_prefix + channel + ".root";		
		TFile* f_data = new TFile(dir + fname);
		if (not f_data -> IsOpen()) {
			cout << "Warning: cannot find " << fname << " in " << dir << endl;
		}
		
		TTree* tree = (TTree*) f_data->Get("ttHtaus/eventTree");

	    fillHistoFromTree(h_, eventList, tree, yields);
	}
	
	h_->SetName("x_"+channel);
	hists.push_back(h_);

	// print out yields
	int icat = 0;
	float ysum = 0.;
	cout << channel << " : ";
	for (const TString lcat : lepCategory) {
		for (const TString bcat : bCategory) {
			cout << yields.at(icat) << "("<<lcat<<","<<bcat<<") ";
			ysum += yields.at(icat++);
		}
	}
	cout << ysum << "(sum)" << endl;
	
	return hists;
}

void fillHistoFromTree(TH1D* h, vector<vector<int>>& eventList, TTree* tree, vector<float>& yields)
// Fill shape histogram from TTree entries.
//Duplicate events characterized by (run, lumisection, event#) are only filled once
{
	
	int run;
	int ls;
	int event;
	float weight;
	float mva_ttbar;
	float mva_ttV;
	int lepCat = 0;
	int btagCat = 0;

	tree->SetBranchAddress("nEvent", &event);
	tree->SetBranchAddress("ls", &ls);
	tree->SetBranchAddress("run", &run);
	tree->SetBranchAddress("MVA_2lss_ttbar", &mva_ttbar);
	tree->SetBranchAddress("MVA_2lss_ttV", &mva_ttV);
	tree->SetBranchAddress("event_weight", &weight);
	//tree->SetBranchAddress("btagCategory", &btagCat);
	//tree->SetBranchAddress("lepCategory", &lepCat);

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

			
			// yields
			yields[2*lepCat+btagCat] += weight;
		}
	}

}

int partition2DBDT(float ttbar, float ttV)
/*
  ICHEP binning
 */
{
	if (ttbar > -1.0 and ttbar <= -0.2 and ttV > -1.0 and ttV <= 1.0)
		return 1;
	else if (ttbar > -0.2 and ttbar <= 0.1 and ttV > -1.0 and ttV <= 1.0)
		return 2;
	else if (ttbar > 0.1 and ttbar <= 0.4) {
		if (ttV > -1.0 and ttV <= 0.3)
			return 3;
		else if (ttV > 0.3 and ttV <= 1.0)
			return 4;
	}
	else if (ttbar > 0.4 and ttbar <= 1.0) {
		if (ttV > -1.0 and ttV <= 0.1)
			return 5;
		else if (ttV > 0.1 and ttV <= 0.4)
			return 6;
		else if (ttV > 0.4 and ttV <= 1.0)
			return 7;
	}

	return 0;
}

TH1D* combineHistCategories(TFile* f, TString name, float nSum,
						   float lumi, float xsection, vector<float>& yields)
{
	TH1D* hist_out = 0;
	bool first = true;

	for (const TString lcat : lepCategory) {
		for (const TString bcat : bCategory) {
			TString hname = "ttHtaus/h_MVA_shape" + lcat + bcat + name;
			TH1D* h_ = (TH1D*)f->Get(hname);
			
			yields.push_back(h_->Integral());

			if (first) {
				hist_out = (TH1D*)h_->Clone();
				first = false;
			}
			else {
				hist_out->Add(h_);
			}

			delete h_;
		}
	}

	// scale
	//hist_out -> Sumw2();
	hist_out -> Scale(lumi * xsection / nSum);

	for (float & y : yields) {
		y *= lumi * xsection / nSum;
	}
	
	return hist_out;
}

TH1D* combineHistCategories(TFile* f, TString name,
						   int nsample, float lumi, float xsection)
{
	vector<float> yields;
	return combineHistCategories(f, name, nsample, lumi, xsection, yields);
}
