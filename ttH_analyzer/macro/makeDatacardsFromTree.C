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
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

#include "Cross_Sections.h"
#include "Misc_Constants.h"

using namespace std;

const float LUMI = 12.9 * 1000.;  // 1/pb

vector<TH1D*> getShapesMC(TString, vector<TString>);
vector<TH1D*> getShapesData(TString, vector<TString>);
map<TString, TH1D*> setupHistoMap(map<TString, TH1D*>&, TString, bool, TString);
void fillHistoFromTreeMC(map<TString, TH1D*>&, TTree*);
void fillHistoFromTreeData(TH1D*, TTree*, vector<vector<int>>&);
void updateHistoName(pair<const TString, TH1D*>&, TString);
vector<TH1D*> getClosureTestShapes(TH1D*);

// hists map
// "", "_gentau", "_faketau", "_"+syts, "_gentau"+"_"+syst, "_faketau"+"_"+syst,
// "_JESUp", "_JESDown"; jesup jesdown on seperate loop

void makeDatacardsFromTree(TString outfile_suffix = "12_9")
{
	vector<TString> channels =
		{"ttH", "TTW", "TTZ", "EWK", "Rares",
		 "fakes_data", "flips_data", "data_obs"
		};
	
	vector<TH1D*> datacards;
	
	for (const auto & channel : channels) {
		cout << "=========================================" << endl;
		cout << "processing channel: " << channel << endl;
		cout << "-----------------------------------------" << endl;

		vector<TH1D*> shapes;
		
		if (channel.Contains("data"))
			shapes = getShapesData(channel, SamplesInChannel.at(channel));
		else {
			cout << "sample" << "\t" << "yields" << "\t" << "yields(gentau)"
				 << "\t" << "yields(faketau)" << endl;
			cout << "-----------------------------------------" << endl;
			shapes = getShapesMC(channel, SamplesInChannel.at(channel));
		}

		datacards.insert(datacards.end(), shapes.begin(), shapes.end());
	} // end of channel loop
	
	// output file
	TString outname = "ttH_2lss_1tau_"+outfile_suffix+".root";
	TFile output_combine(outname, "RECREATE");
	
	for (auto & h : datacards)
		h->Write();

	cout << "Output file: " << outname << endl;

	return;
}

map<TString, TH1D*> setupHistoMap(TString sample, bool addSyst, TString jec_suffix)
{
	map<TString, TH1D*> hists;

	assert(jec_suffix == "" or jec_suffix == "_JESUp" or
		   jec_suffix == "_JESDown");
	
	vector<TString> pre_keys;

	pre_keys.push_back("");
	pre_keys.push_back("_gentau");
	pre_keys.push_back("_faketau");

	if (sample=="ttH") {
		for (TString mode : HiggsDecayMode) {
			pre_keys.push_back("_"+mode);
			pre_keys.push_back("_"+mode+"_gentau");
			pre_keys.push_back("_"+mode+"_faketau");
		}
	}	
	
	for (auto & pre : pre_keys) {

		vector<TString> keys_sys;
		
		keys_sys.push_back("" + jec_suffix);

		if (addSyst and jec_suffix == "") {  // only for the central one
			if (sample=="ttH" or sample=="TTW" or sample=="TTZ") {
				for (auto & thu : ThSysts)
					keys_sys.push_back("_"+thu);
			}
			
			for (auto & btagsys : BTagSysts)
				keys_sys.push_back("_"+btagsys);
		}

		for (auto & ks : keys_sys) {
			TH1D* h_s = new TH1D("x_"+sample+pre+ks, "", 7, 0.5, 7.5);
			h_s->Sumw2();
			hists[pre+ks] = h_s;
		}
	}
	
	return hists;
}

vector<TH1D*> getShapesMC( TString channel, vector<TString> samples)
{
	vector<TH1D*> shapes;
    float yields[3] = {0., 0., 0.};

	bool first = true;
	
	for (const TString & sample : samples) {
		
		cout << sample << "\t";
		
		// open files
		TFile* f_central =
			new TFile(dir_map.at(sample)+"output_"+sample+".root");
		TFile* f_jesup =
			new TFile(dir_map.at(sample)+"output_"+sample+"_jesup.root");
		TFile* f_jesdown =
			new TFile(dir_map.at(sample)+"output_"+sample+"_jesdown.root");

		// Check if files are open
		if (not f_central->IsOpen())
			cout << "CANNOT open file " << "output_" << sample << ".root";
		if (not f_jesup->IsOpen())
			cout << "CANNOT open file " << "output_" << sample << "_jesup.root";
		if (not f_jesdown->IsOpen())
			cout << "CANNOT open file " << "output_" << sample << "_jesdown.root";

		// get trees
		TTree* tree_central = (TTree*) f_central->Get("ttHtaus/eventTree");
		TTree* tree_jesup = (TTree*) f_jesup->Get("ttHtaus/eventTree");
		TTree* tree_jesdown = (TTree*) f_jesdown->Get("ttHtaus/eventTree");

		// set up histograms
		map<TString,TH1D*> hists = setupHistoMap(sample, true, "");
		map<TString,TH1D*> hists_jesup = setupHistoMap(sample, false, "_JESUp");
		map<TString,TH1D*> hists_jesdown = setupHistoMap(sample, false, "_JESDown");

		fillHistoFromTreeMC(hists, tree_central);
		fillHistoFromTreeMC(hists_jesup, tree_jesup);
		fillHistoFromTreeMC(hists_jesdown, tree_jesdown);
		
		// append maps
		hists.insert(hists_jesup.begin(), hists_jesup.end());
		hists.insert(hists_jesdown.begin(), hists_jesdown.end());
		
		// scale histograms
		//TH1D* h_SumGenWeight = (TH1D*)f_central->Get("ttHtaus/h_SumGenWeight");
		TH1D* h_SumGenWeightxPU = (TH1D*)f_central->Get("ttHtaus/h_SumGenWeightxPU");
		//float nSum = h_SumGenWeight->GetBinContent(1);
		float nSum = h_SumGenWeightxPU->GetBinContent(1);
		float XS = xsection::xsection[string(sample)];

		int ih = 0;
		for (auto & h : hists) {
			//h.second->Sumw2();
			h.second->Scale(LUMI * XS / nSum);
			
			// update histogram names and add to the output histograms
			if (first) {
				updateHistoName(h, channel);
				shapes.push_back(h.second);
			}
			else {
				shapes.at(ih++)->Add(h.second);
			}
		}		

		first = false;
		
		// print out yields
		if (sample == "ttH") {
			cout << endl;
			for (TString mode : HiggsDecayMode) {
				cout << "ttH_" << mode << "\t";
				cout << setw(4) << hists["_"+mode]->Integral() << "\t";
				cout << setw(4) << hists["_"+mode+"_gentau"]->Integral() << "\t";
				cout << setw(4) << hists["_"+mode+"_faketau"]->Integral() << endl;
			}
			cout << "\t";
		}
		
		cout << setw(4) << hists[""]->Integral() << "\t";
		cout << setw(4) << hists["_gentau"]->Integral() << "\t";
		cout << setw(4) << hists["_faketau"]->Integral() << endl;

		yields[0] += hists[""]->Integral();
		yields[1] += hists["_gentau"]->Integral();
		yields[2] += hists["_faketau"]->Integral();
		
		// close file and release memory
		//f_central->Close();
		//f_jesup->Close();
		//f_jesdown->Close();

		//delete f_central;
		//delete f_jesup;
		//delete f_jesdown;
		
	} // end of sample loop

	if (channel == "Rares" or channel == "EWK") {
		cout << "Total" << setw(4) << "\t" << yields[0] << "\t" << yields[1]
			 << "\t" << yields[2] << endl;
	}
	
	return shapes;
}

vector<TH1D*> getShapesData(TString channel, vector<TString> samples)
{
	vector<TH1D*> shapes;

	vector<vector<int>> eventList; // (run, lumisection, event)

	TH1D* h = new TH1D("x_"+channel, "", 7, 0.5, 7.5);
	h->Sumw2();

	int nevents = eventList.size();
	cout << "eventList size : " << nevents << endl;
	for (const TString & sample : samples) {

		cout << "opening sample: " << sample << endl;

		// open file
		TFile* f = new TFile(dir_map.at(sample)+"output_"+channel+".root");

		// get trees
		TTree* tree = (TTree*) f->Get("ttHtaus/eventTree");

		fillHistoFromTreeData(h, tree, eventList);
		
		// close file
		//f->Close();
		//delete f;
	}
	
	// print out yields
	cout << channel << "\t" << setw(4) << h->Integral() << endl;
	
	shapes.push_back(h);

	// closure test systematics of the lepton fake rate
	if (channel == "fakes_data") {
		vector<TH1D*> hists_clos = getClosureTestShapes(h);
		shapes.insert(shapes.end(), hists_clos.begin(), hists_clos.end());
	}
	
	return shapes;
}

void fillHistoFromTreeMC(map<TString, TH1D*>& hists, TTree* tree)
{
	assert(hists.size());
	
	int run;
	int ls;
	int event;
	
	float mva_ttbar;
	float mva_ttV;

	float event_weight;
	float pu_weight;
	float mc_weight;
	float btagsf_weight;
	float leptonsf_weight;
	float tausf_weight;
	float triggersf_weight;

	float mc_weight_scale_muf0p5;
	float mc_weight_scale_muf2;
	float mc_weight_scale_mur0p5;
	float mc_weight_scale_mur2;
	
	float btagsf_weight_lfup;
	float btagsf_weight_lfdown;
	float btagsf_weight_hfup;
	float btagsf_weight_hfdown;
	float btagsf_weight_hfstats1up;
	float btagsf_weight_hfstats1down;
	float btagsf_weight_hfstats2up;
	float btagsf_weight_hfstats2down;
	float btagsf_weight_lfstats1up;
	float btagsf_weight_lfstats1down;
	float btagsf_weight_lfstats2up;
	float btagsf_weight_lfstats2down;
	float btagsf_weight_cerr1up;
	float btagsf_weight_cerr1down;
	float btagsf_weight_cerr2up;
	float btagsf_weight_cerr2down;

	int isGenMatched;
	int HiggsDecayType;

	int lepCategory;
	int btagCategory;

	int matchHLTPath;

	int ibin;

	tree->SetBranchAddress("run", &run);
	tree->SetBranchAddress("ls", &ls);
	tree->SetBranchAddress("nEvent", &event);
	tree->SetBranchAddress("MVA_2lss_ttbar", &mva_ttbar);
	tree->SetBranchAddress("MVA_2lss_ttV", &mva_ttV);
	tree->SetBranchAddress("event_weight", &event_weight);
	tree->SetBranchAddress("PU_weight", &pu_weight);
	tree->SetBranchAddress("MC_weight", &mc_weight);
	tree->SetBranchAddress("bTagSF_weight", &btagsf_weight);
	tree->SetBranchAddress("leptonSF_weight", &leptonsf_weight);
	tree->SetBranchAddress("tauSF_weight", &tausf_weight);
	tree->SetBranchAddress("triggerSF_weight", &triggersf_weight);
	tree->SetBranchAddress("MC_weight_scale_muF0p5", &mc_weight_scale_muf0p5);
	tree->SetBranchAddress("MC_weight_scale_muF2", &mc_weight_scale_muf2);
	tree->SetBranchAddress("MC_weight_scale_muR0p5", &mc_weight_scale_mur0p5);
	tree->SetBranchAddress("MC_weight_scale_muR2", &mc_weight_scale_mur2);
	tree->SetBranchAddress("btagSF_weight_LFUp", &btagsf_weight_lfup);
	tree->SetBranchAddress("btagSF_weight_LFDown", &btagsf_weight_lfdown);
	tree->SetBranchAddress("btagSF_weight_HFUp", &btagsf_weight_hfup);
	tree->SetBranchAddress("btagSF_weight_HFDown", &btagsf_weight_hfdown);
	tree->SetBranchAddress("btagSF_weight_HFStats1Up", &btagsf_weight_hfstats1up);
	tree->SetBranchAddress("btagSF_weight_HFStats1Down", &btagsf_weight_hfstats1down);
	tree->SetBranchAddress("btagSF_weight_HFStats2Up", &btagsf_weight_hfstats2up);
	tree->SetBranchAddress("btagSF_weight_HFStats2Down", &btagsf_weight_hfstats2down);
	tree->SetBranchAddress("btagSF_weight_LFStats1Up", &btagsf_weight_lfstats1up);
	tree->SetBranchAddress("btagSF_weight_LFStats1Down", &btagsf_weight_lfstats1down);
	tree->SetBranchAddress("btagSF_weight_LFStats2Up", &btagsf_weight_lfstats2up);
	tree->SetBranchAddress("btagSF_weight_LFStats2Down", &btagsf_weight_lfstats2down);
	tree->SetBranchAddress("btagSF_weight_cErr1Up", &btagsf_weight_cerr1up);
	tree->SetBranchAddress("btagSF_weight_cErr1Down", &btagsf_weight_cerr1down);
	tree->SetBranchAddress("btagSF_weight_cErr2Up", &btagsf_weight_cerr2up);
	tree->SetBranchAddress("btagSF_weight_cErr2Down", &btagsf_weight_cerr2down);
	tree->SetBranchAddress("isGenMatched", &isGenMatched);
	tree->SetBranchAddress("HiggsDecayType", &HiggsDecayType);
	tree->SetBranchAddress("lepCategory", &lepCategory);
	tree->SetBranchAddress("btagCategory", &btagCategory);
	tree->SetBranchAddress("matchHLTPath", &matchHLTPath);
	tree->SetBranchAddress("ibin", &ibin);

	// Loop over events in the TTree
	int nEntries = tree->GetEntries();

	for (int i = 0; i < nEntries; ++i) {
		tree->GetEntry(i);

		if (not matchHLTPath) continue;
		
		// Update weights here if needed

		
		// Update bin index if needed
		// ibin = partition2DBDT(mva_ttbar, mva_ttV);

		// Fill the histograms
		// Assume keys of the histogram map are in the following format:
		// (_HiggsDecayMode)+(_gentau/_faketau)+(_systSuffix)
		// empty string "" for the central inclusive one
		for (auto & hist : hists) {
			TString key = hist.first;
			// filters
			if (key.Contains("_htt") and abs(HiggsDecayType) != 15) continue;
			if (key.Contains("_hzz") and abs(HiggsDecayType) != 23) continue;
			if (key.Contains("_hww") and abs(HiggsDecayType) != 24) continue;
			if (key.Contains("_gentau") and not isGenMatched) continue;
			if (key.Contains("_faketau") and isGenMatched) continue;

			// decide which weight to use
			float w_sf = 1;
			if (key.EndsWith("_LFUp"))
				w_sf = btagsf_weight_lfup / btagsf_weight;
			else if (key.EndsWith("_LFDown"))
				w_sf = btagsf_weight_lfdown / btagsf_weight;
			else if (key.EndsWith("_HFUp"))
				w_sf = btagsf_weight_hfup / btagsf_weight;
			else if (key.EndsWith("_HFDown"))
				w_sf = btagsf_weight_hfdown / btagsf_weight;
			else if (key.EndsWith("_HFStats1Up"))
				w_sf = btagsf_weight_hfstats1up / btagsf_weight;
			else if (key.EndsWith("_HFStats1Down"))
				w_sf = btagsf_weight_hfstats1down / btagsf_weight;
			else if (key.EndsWith("_HFStats2Up"))
				w_sf = btagsf_weight_hfstats2up / btagsf_weight;
			else if (key.EndsWith("_HFStats2Down"))
				w_sf = btagsf_weight_hfstats2down / btagsf_weight;
			else if (key.EndsWith("_LFStats1Up"))
				w_sf = btagsf_weight_lfstats1up / btagsf_weight;
			else if (key.EndsWith("_LFStats1Down"))
				w_sf = btagsf_weight_lfstats1down / btagsf_weight;
			else if (key.EndsWith("_LFStats2Up"))
				w_sf = btagsf_weight_lfstats2up / btagsf_weight;
			else if (key.EndsWith("_LFStats2Down"))
				w_sf = btagsf_weight_lfstats2down / btagsf_weight;
			else if (key.EndsWith("_cErr1Up"))
				w_sf = btagsf_weight_cerr1up / btagsf_weight;
			else if (key.EndsWith("_cErr1Down"))
				w_sf = btagsf_weight_cerr1down / btagsf_weight;
			else if (key.EndsWith("_cErr2Up"))
				w_sf = btagsf_weight_cerr2up / btagsf_weight;
			else if (key.EndsWith("_cErr2Down"))
				w_sf = btagsf_weight_cerr2down / btagsf_weight;
			else if (key.EndsWith("_x1Up"))
				w_sf = mc_weight_scale_muf2 / mc_weight;
			else if (key.EndsWith("_x1Down"))
				w_sf = mc_weight_scale_muf0p5 / mc_weight;
			else if (key.EndsWith("_y1Up"))
				w_sf = mc_weight_scale_mur2 / mc_weight;
			else if (key.EndsWith("_y2Down"))
				w_sf = mc_weight_scale_mur0p5 / mc_weight;

			hist.second -> Fill(ibin, event_weight * w_sf);
		} // end of keys loop
		
	} // end of event loop

	return;
}

void fillHistoFromTreeData(TH1D* h, TTree* tree, vector<vector<int>>& eventList)
{
	int run;
	int ls;
	int event;
	// FIXME: should be unsigned long long;
	float event_weight;
	float mva_ttbar;
	float mva_ttV;
	int ibin;
	int matchHLTPath;

	tree->SetBranchAddress("nEvent", &event);
	tree->SetBranchAddress("ls", &ls);
	tree->SetBranchAddress("run", &run);
	tree->SetBranchAddress("event_weight", &event_weight);
	tree->SetBranchAddress("MVA_2lss_ttbar", &mva_ttbar);
	tree->SetBranchAddress("MVA_2lss_ttV", &mva_ttV);
	tree->SetBranchAddress("ibin", &ibin);
	tree->SetBranchAddress("matchHLTPath", &matchHLTPath);

	int nEntries = tree->GetEntries();

	for (int i = 0; i < nEntries; ++i) {
		tree->GetEntry(i);

		if (not matchHLTPath) continue;
		
		vector<int> eventid = {run, ls, event};
		//assert(run >= 0 and ls >= 0 and event >= 0);

		bool alreadyIncluded =
			find(eventList.begin(), eventList.end(), eventid) != eventList.end();

		if (not alreadyIncluded) {

			eventList.push_back(eventid);

			// update event weight here if needed

			// update bin index here if needed

			// fill histogram
			h->Fill(ibin, event_weight);
		}	
	} // end of event loop
}

void updateHistoName(pair<const TString, TH1D*>& h, TString channel)
{
	TString key = h.first;
	
	for (const TString & bsys : BTagSysts) {
		if (key.EndsWith(bsys)) {
			key.ReplaceAll(bsys, sys_coname+"btag_"+bsys);
			h.second->SetName("x_"+channel+key);
			return;
		}
	}

	for (const TString & thu : ThSysts) {
		if (key.EndsWith(thu)) {
			key.ReplaceAll(thu, sys_coname+"thu_shape_"+channel+"_"+thu);
			h.second->SetName("x_"+channel+key);
			return;
		}
	}

	if (key.EndsWith("JESUp")) {
		key.ReplaceAll("JESUp", sys_coname+"JESUp");
		h.second->SetName("x_"+channel+key);
		return;
	}
	
	if (key.EndsWith("JESDown")) {
		key.ReplaceAll("JESDown", sys_coname+"JESDown");
		h.second->SetName("x_"+channel+key);
		return;
	}

	h.second->SetName("x_"+channel+key);
	return;
}

vector<TH1D*> getClosureTestShapes(TH1D* h_nominal)
{	
	TString hname = h_nominal->GetName();
	assert(hname=="x_fakes_data");
	//h_nominal->Sumw2();
	
	vector<TH1D*> hists_clos;
	
	// open file
	TString era_clos;
	if (LUMI == 12.9 * 1000)
		era_clos = "12.9fb";
	else if (LUMI == 36.8 * 1000)
		era_clos = "36.8fb";
	else {
		cerr << "Closure test file is not available!" << endl;
		assert(0);
	}
	
	TFile* f_clos = new TFile("../data/Closure_FR_syst/Closure_FR_lepton_syst_2lSS1tau_nofaketau_MVA_2lSS_"+era_clos+".root", "READ");
	
	// get histograms
	TH1D* h_ttbar_minus_qcd_fr_ele =
		(TH1D*)f_clos->Get("x_TT_DL_FR_TT_MC_minus_FR_QCD_MC_ele");
	TH1D* h_ttbar_minus_qcd_fr_mu =
		(TH1D*)f_clos->Get("x_TT_DL_FR_TT_MC_minus_FR_QCD_MC_mu");

	//h_ttbar_minus_qcd_fr_ele->Sumw2();
	//h_ttbar_minus_qcd_fr_mu->Sumw2();
	
	TH1D* h_clos_e_shape_up =
		new TH1D("x_fakes_data_"+sys_coname+"Clos_e_shape1Up","", 7, 0.5, 7.5);
	TH1D* h_clos_e_shape_down =
		new TH1D("x_fakes_data_"+sys_coname+"Clos_e_shape1Down","", 7, 0.5, 7.5);
	TH1D* h_clos_mu_shape_up =
		new TH1D("x_fakes_data_"+sys_coname+"Clos_mu_shape1Up","", 7, 0.5, 7.5);
	TH1D* h_clos_mu_shape_down =
		new TH1D("x_fakes_data_"+sys_coname+"Clos_mu_shape1Down","", 7, 0.5, 7.5);

	h_clos_e_shape_up->Sumw2();
	h_clos_e_shape_down->Sumw2();
	h_clos_mu_shape_up->Sumw2();
	h_clos_mu_shape_down->Sumw2();
	
	h_clos_e_shape_up->Add(h_nominal, h_ttbar_minus_qcd_fr_ele, 1, 1);
	h_clos_e_shape_down->Add(h_nominal, h_ttbar_minus_qcd_fr_ele, 1, -1);
	h_clos_mu_shape_up->Add(h_nominal, h_ttbar_minus_qcd_fr_mu, 1, 1);
	h_clos_mu_shape_down->Add(h_nominal, h_ttbar_minus_qcd_fr_mu, 1, -1);

	hists_clos.push_back(h_clos_e_shape_up);
	hists_clos.push_back(h_clos_e_shape_down);
	hists_clos.push_back(h_clos_mu_shape_up);
	hists_clos.push_back(h_clos_mu_shape_down);

	return hists_clos;
}
