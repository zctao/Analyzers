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
#include "Analyzers/ttH_analyzer/interface/Types_enum.h"
#include "Analyzers/ttH_analyzer/interface/TreeAnalyzer.h"

//#include "eventSelector.h"

//const float LUMI = 12.9 * 1000.;  // 1/pb
const float LUMI = 35.867 * 1000; // 1/pb

vector<TH1D*> getShapesMC(TString, vector<TString>, bool, int);
vector<TH1D*> getShapesData(TString, vector<TString>, bool, int);
map<TString, TH1D*> setupHistoMap(map<TString, TH1D*>&, TString, bool, TString);
//void fillHistoFromTreeMC(map<TString, TH1D*>&, TTree*);
//void fillHistoFromTreeData(TH1D*, TTree*, vector<vector<unsigned long long>>&);
void updateHistoName(pair<const TString, TH1D*>&, TString);
vector<TH1D*> getClosureTestShapes(TH1D*);
void makeBinContentsPositive(TH1*, int);
double compIntegral(TH1*, bool, bool);
double square(double);

// hists map
// "", "_gentau", "_faketau", "_"+syts, "_gentau"+"_"+syst, "_faketau"+"_"+syst,
// "_JESUp", "_JESDown"; jesup jesdown on seperate loop

void makeDatacardsFromTree(TString outfile_suffix = "2016b-h", bool addSyst = true, int verbosity=0)
// 0=unverbose; 1=show passed events; 2=show failed events;
{
	using namespace std;

	TH1::AddDirectory(0);
	
	gROOT->ProcessLine(".L ../src/SFHelper.cc+");
	gROOT->ProcessLine(".L ../src/TreeAnalyzer.cc+");
	
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
			shapes = getShapesData(channel, SamplesInChannel.at(channel),
								   addSyst, verbosity);
		else {
			cout << "sample" << "\t" << "yields" << "\t" << "yields(gentau)"
				 << "\t" << "yields(faketau)" << endl;
			cout << "-----------------------------------------" << endl;
			shapes = getShapesMC(channel, SamplesInChannel.at(channel),
								 addSyst, verbosity);
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

map<TString, TH1D*> setupHistoMap(TString sample, bool addSyst, TString suffix)
{
	map<TString, TH1D*> hists;

	assert(suffix == "" or suffix == "_JESUp" or suffix == "_JESDown" or
		   suffix == "_TESUp" or suffix == "_TESDown");
	
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
		
		keys_sys.push_back("" + suffix);

		if (addSyst and suffix == "") {  // only for the central one
			if (sample=="ttH" or sample=="TTW" or sample=="TTZ") {
				for (auto & thu : ThSysts)
					keys_sys.push_back("_"+thu);
			}
			
			for (auto & btagsys : BTagSysts)
				keys_sys.push_back("_"+btagsys);

			for (auto & frjtsys : FRjtSysts) {
				if (pre.Contains("_faketau"))
					keys_sys.push_back("_"+frjtsys);
			}
		}

		for (auto & ks : keys_sys) {
			TH1D* h_s = new TH1D("x_"+sample+pre+ks, "", 7, 0.5, 7.5);
			h_s->Sumw2();
			hists[pre+ks] = h_s;
		}
	}
	
	return hists;
}

vector<TH1D*> getShapesMC( TString channel, vector<TString> samples, bool addSyst,
						   int verbosity)
{	
	Analysis_types AnaType = Analysis_types::Analyze_2lss1tau;
	Selection_types SelType = Selection_types::Signal_2lss1tau;
	
	vector<TH1D*> shapes;
    float yields[3] = {0., 0., 0.};

	bool first = true;
	
	for (const TString & sample : samples) {
		
		cout << sample << "\t";

		// set up histograms
		map<TString,TH1D*> hists = setupHistoMap(sample, addSyst, "");

		// open files
		TFile* f_central =
			new TFile(dir_map.at(sample)+"output_"+sample+".root");

		// check if file is open
		if (f_central->IsOpen()) {
			// get trees
			TTree* tree_central = (TTree*) f_central->Get("ttHtaus/eventTree");
			//fillHistoFromTreeMC(hists, tree_central);
			TreeAnalyzer tana(tree_central,AnaType,SelType,false,verbosity);
			tana.fill_Datacards_MC(hists);
			//tana.dump_Events(sample);
		}
		else
			cout << "CANNOT open file " << "output_" << sample << ".root" << endl;
		
		if (addSyst) {
			// set up histograms
			map<TString,TH1D*> hists_jesup = setupHistoMap(sample, false, "_JESUp");
			map<TString,TH1D*> hists_jesdown = setupHistoMap(sample, false, "_JESDown");
			map<TString,TH1D*> hists_tesup = setupHistoMap(sample, false, "_TESUp");
			map<TString,TH1D*> hists_tesdown = setupHistoMap(sample, false, "_TESDown");

			// open files
			TFile* f_jesup =
				new TFile(dir_map.at(sample)+"output_"+sample+"_jesup.root");
			TFile* f_jesdown =
				new TFile(dir_map.at(sample)+"output_"+sample+"_jesdown.root");
			TFile* f_tesup =
				new TFile(dir_map.at(sample)+"output_"+sample+"_tesup.root");
			TFile* f_tesdown =
				new TFile(dir_map.at(sample)+"output_"+sample+"_tesdown.root");

			// check if file is open
			if (f_jesup->IsOpen()) {
				// get trees
				TTree* tree_jesup = (TTree*) f_jesup->Get("ttHtaus/eventTree");
				//fillHistoFromTreeMC(hists_jesup, tree_jesup);
				TreeAnalyzer tana_jesup(tree_jesup,AnaType,SelType,false,verbosity);
				tana_jesup.fill_Datacards_MC(hists_jesup);
			}
			else
				cout << "CANNOT open file " << "output_" << sample << "_jesup.root" << endl;
			
			if (f_jesdown->IsOpen()) {
				// get trees
				TTree* tree_jesdown = (TTree*) f_jesdown->Get("ttHtaus/eventTree");
				//fillHistoFromTreeMC(hists_jesdown, tree_jesdown);
				TreeAnalyzer tana_jesdown(tree_jesdown,AnaType,SelType,false,verbosity);
				tana_jesdown.fill_Datacards_MC(hists_jesdown);
			}
			else
				cout << "CANNOT open file " << "output_" << sample << "_jesdown.root" << endl;
			
			if (f_tesup->IsOpen()) {
				// get trees
				TTree* tree_tesup = (TTree*) f_tesup->Get("ttHtaus/eventTree");
				//fillHistoFromTreeMC(hists_tesup, tree_tesup);
				TreeAnalyzer tana_tesup(tree_tesup,AnaType,SelType,false,verbosity);
				tana_tesup.fill_Datacards_MC(hists_tesup);
			}
			else
				cout << "CANNOT open file " << "output_" << sample << "_tesup.root" << endl;
			
			if (f_tesdown->IsOpen()) {
				// get trees
				TTree* tree_tesdown = (TTree*) f_tesdown->Get("ttHtaus/eventTree");
				//fillHistoFromTreeMC(hists_tesdown, tree_tesdown);
				TreeAnalyzer tana_tesdown(tree_tesdown,AnaType,SelType,false,verbosity);
				tana_tesdown.fill_Datacards_MC(hists_tesdown);
			}
			else
				cout << "CANNOT open file " << "output_" << sample << "_tesdown.root" << endl;
		
			// append maps
			hists.insert(hists_jesup.begin(), hists_jesup.end());
			hists.insert(hists_jesdown.begin(), hists_jesdown.end());
			hists.insert(hists_tesup.begin(), hists_tesup.end());
			hists.insert(hists_tesdown.begin(), hists_tesdown.end());
		}
		
		// scale histograms
		if (f_central->IsOpen()) {
			//TH1D* h_SumGenWeight = (TH1D*)f_central->Get("ttHtaus/h_SumGenWeight");
			TH1D* h_SumGenWeightxPU = (TH1D*)f_central->Get("ttHtaus/h_SumGenWeightxPU");
			//float nSum = h_SumGenWeight->GetBinContent(1);
			float nSum = h_SumGenWeightxPU->GetBinContent(1);
			float XS = xsection::xsection[string(sample)];
			
			int ih = 0;
			for (auto & h : hists) {
				//h.second->Sumw2();
				h.second->Scale(LUMI * XS / nSum);

				// make bin content positive
				makeBinContentsPositive(h.second,0);
				
				// update histogram names and add to the output histograms
				if (first) {
					updateHistoName(h, channel);
					shapes.push_back(h.second);
				}
				else {
					shapes.at(ih++)->Add(h.second);
				}
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
				cout << hists["_"+mode]->GetEntries() << "\t";
				cout << hists["_"+mode+"_gentau"]->GetEntries() << "\t";
				cout << hists["_"+mode+"_faketau"]->GetEntries() << endl;
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

vector<TH1D*> getShapesData(TString channel, vector<TString> samples,
							bool addSyst, int verbosity)
{	
	Analysis_types AnaType = Analysis_types::Analyze_2lss1tau;
	Selection_types SelType = Selection_types::Signal_2lss1tau;
	
	if (channel.Contains("fakes"))
		SelType = Selection_types::Control_1lfakeable;
	else if (channel.Contains("flips"))
		SelType = Selection_types::Control_2los1tau;
	
	vector<TH1D*> shapes;
	vector<vector<unsigned long long>> eventList; // (run, lumisection, event)
	vector<vector<unsigned long long>> eventList_dump;

	vector<TH1D*> hists;
	TH1D* h = new TH1D("x_"+channel, "", 7, 0.5, 7.5);
	h->Sumw2();
	hists.push_back(h);

	if (SelType == Selection_types::Control_1lfakeable and addSyst) {
		for (const auto & frlsyst : FRlSysts) {
			TH1D* h_s = new TH1D("x_"+channel+"_"+sys_coname+frlsyst, "",
								 7, 0.5, 7.5);
			h_s->Sumw2();
			hists.push_back(h_s);
		}
	}

	int nevents = eventList.size();
	cout << "eventList size : " << nevents << endl;
	for (const TString & sample : samples) {

		cout << "opening sample: " << sample << endl;

		// open file
		TFile* f = new TFile(dir_map.at(sample)+"output_"+channel+".root");

		if (f->IsOpen()) {
			// get trees
			TTree* tree = (TTree*) f->Get("ttHtaus/eventTree");
			//fillHistoFromTreeData(h, tree, eventList);
			TreeAnalyzer tana(tree,AnaType,SelType,true,verbosity);
			tana.fill_Datacards_Data(hists, eventList);
			//tana.dump_Events(channel+"_"+sample, eventList_dump);
		}
		else {
			cout << "CANNOT open file " << "output_" << channel
				 << ".root for sample " << sample << endl;
		}
		// close file
		//f->Close();
		//delete f;
	}
	
	// print out yields
	cout << channel << "\t" << setw(4) << hists.at(0)->Integral() << endl;

	// make bin content positive
	for (auto h_ : hists) {
		makeBinContentsPositive(h_,0);
		shapes.push_back(h_);
	}

	// closure test systematics of the lepton fake rate
	if (channel == "fakes_data" and addSyst) {
		vector<TH1D*> hists_clos = getClosureTestShapes(h);
		shapes.insert(shapes.end(), hists_clos.begin(), hists_clos.end());
	}
	
	return shapes;
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
			
			TString ch = channel;
			if (channel.Contains("TTW"))
				ch = "ttW";
			if (channel.Contains("TTZ"))
				ch = "ttZ";
			
			key.ReplaceAll(thu, sys_coname+"thu_shape_"+ch+"_"+thu);
			h.second->SetName("x_"+channel+key);
			return;
		}
	}

	for (const TString & frjt : FRjtSysts) {
		if (key.EndsWith(frjt)) {
			key.ReplaceAll(frjt, sys_coname+frjt);
			h.second->SetName("x_"+channel+key);
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

	if (key.EndsWith("TESUp")) {
		key.ReplaceAll("TESUp", sys_coname+"tauESUp");
		h.second->SetName("x_"+channel+key);
		return;
	}
	
	if (key.EndsWith("TESDown")) {
		key.ReplaceAll("TESDown", sys_coname+"tauESDown");
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
	if (data_lumi == "12_9fb/")
		era_clos = "12.9fb";
	else if (data_lumi == "36_8fb/")
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
		new TH1D("x_fakes_data_"+sys_coname+"Clos_e_shapeUp","", 7, 0.5, 7.5);
	TH1D* h_clos_e_shape_down =
		new TH1D("x_fakes_data_"+sys_coname+"Clos_e_shapeDown","", 7, 0.5, 7.5);
	TH1D* h_clos_mu_shape_up =
		new TH1D("x_fakes_data_"+sys_coname+"Clos_m_shapeUp","", 7, 0.5, 7.5);
	TH1D* h_clos_mu_shape_down =
		new TH1D("x_fakes_data_"+sys_coname+"Clos_m_shapeDown","", 7, 0.5, 7.5);

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

void makeBinContentsPositive(TH1* histogram, int verbosity)

{
	if ( verbosity ) {
		std::cout << "<makeBinContentsPositive>:" << std::endl;
		std::cout << " integral(" << histogram->GetName() << ") = " << histogram->Integral() << std::endl;
	}
	
	double integral_original = compIntegral(histogram, true, true);
	
	if ( integral_original < 0. ) integral_original = 0.;
	
	if ( verbosity ) {
		std::cout << " integral_original = " << integral_original << std::endl;
	}
	
	int numBins = histogram->GetNbinsX();
	
	for ( int iBin = 0; iBin <= (numBins + 1); ++iBin ) {
		double binContent_original = histogram->GetBinContent(iBin);
		double binError2_original = square(histogram->GetBinError(iBin));
		
		if ( binContent_original < 0. ) {
			double binContent_modified = 0.;
			double binError2_modified = binError2_original + square(binContent_original - binContent_modified);
			
			assert(binError2_modified >= 0.);
			
			if ( verbosity ) {
				std::cout << "bin #" << iBin << " (x =  " << histogram->GetBinCenter(iBin) << "): binContent = " << binContent_original << " +/- " << TMath::Sqrt(binError2_original) << " --> setting it to binContent = " << binContent_modified << " +/- " << TMath::Sqrt(binError2_modified) << std::endl;
			}

			histogram->SetBinContent(iBin, binContent_modified);
			histogram->SetBinError(iBin, TMath::Sqrt(binError2_modified));
		}
	}

	double integral_modified = compIntegral(histogram, true, true);
	
	if ( integral_modified < 0. ) integral_modified = 0.;
	
	if ( verbosity ) {
		std::cout << " integral_modified = " << integral_modified << std::endl;
	}
	
	if ( integral_modified > 0. ) {
		double sf = integral_original/integral_modified;
		
		if ( verbosity ) {
			std::cout << "--> scaling histogram by factor = " << sf << std::endl;
		}
		
		histogram->Scale(sf);
	} else {
		for ( int iBin = 0; iBin <= (numBins + 1); ++iBin ) {
			histogram->SetBinContent(iBin, 0.);
		}
	}

	if ( verbosity ) {
		std::cout << " integral(" << histogram->GetName() << ") = " << histogram->Integral() << std::endl;
	}
}

double compIntegral(TH1* histogram, bool includeUnderflowBin, bool includeOverflowBin)
{
	double sumBinContent = 0.;
	int numBins = histogram->GetNbinsX();
	int firstBin = ( includeUnderflowBin ) ? 0 : 1;
	int lastBin = ( includeOverflowBin  ) ? (numBins + 1) : numBins;
	
	for ( int iBin = firstBin; iBin <= lastBin; ++iBin ) {
		sumBinContent += histogram->GetBinContent(iBin);
	}
	
	return sumBinContent;
}

double square(double x)
{
	return x*x;
}
