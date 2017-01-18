#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>

using namespace std;

void makeSyncTree(TString input_file="/uscms/home/ztao/nobackup/output_sync.root",
				  TString output_file="~/nobackup/ttHTT_syncNtuple/80X/syncNtuple.root")
{
	// open input file and read tree
	cout << "Opening input file: " << input_file << endl;
	TFile* old_file = new TFile(input_file);
	TTree* old_tree = (TTree*)old_file->Get("ttHtaus/eventTree");

	// turn off extra branches
	old_tree->SetBranchStatus("event_weight",0);
	old_tree->SetBranchStatus("MC_weight_scale_muF0p5", 0);
	old_tree->SetBranchStatus("MC_weight_scale_muF2", 0);
	old_tree->SetBranchStatus("MC_weight_scale_muR0p5", 0);
	old_tree->SetBranchStatus("MC_weight_scale_muR2", 0);
	old_tree->SetBranchStatus("btagSF_weight_LFUp", 0);
	old_tree->SetBranchStatus("btagSF_weight_LFDown", 0);
	old_tree->SetBranchStatus("btagSF_weight_HFUp", 0);
	old_tree->SetBranchStatus("btagSF_weight_HFDown", 0);
	old_tree->SetBranchStatus("btagSF_weight_HFStats1Up", 0);
	old_tree->SetBranchStatus("btagSF_weight_HFStats1Down", 0);
	old_tree->SetBranchStatus("btagSF_weight_HFStats2Up", 0);
	old_tree->SetBranchStatus("btagSF_weight_HFStats2Down", 0);
	old_tree->SetBranchStatus("btagSF_weight_LFStats1Up", 0);
	old_tree->SetBranchStatus("btagSF_weight_LFStats1Down", 0);
	old_tree->SetBranchStatus("btagSF_weight_LFStats2Up", 0);
	old_tree->SetBranchStatus("btagSF_weight_LFStats2Down", 0);
	old_tree->SetBranchStatus("btagSF_weight_cErr1Up", 0);
	old_tree->SetBranchStatus("btagSF_weight_cErr1Down", 0);
	old_tree->SetBranchStatus("btagSF_weight_cErr2Up", 0);
	old_tree->SetBranchStatus("btagSF_weight_cErr2Down", 0);
	old_tree->SetBranchStatus("HiggsDecayType", 0);
	old_tree->SetBranchStatus("lepCategory", 0);
	old_tree->SetBranchStatus("btagCategory", 0);
	old_tree->SetBranchStatus("npuTrue", 0);
	old_tree->SetBranchStatus("npuInTime", 0);
	old_tree->SetBranchStatus("n_tau", 0);
	old_tree->SetBranchStatus("ibin", 0);
	old_tree->SetBranchStatus("tau0_decayMode", 0);
	old_tree->SetBranchStatus("tau1_decayMode", 0);
	old_tree->SetBranchStatus("pass_single_e", 0);
	old_tree->SetBranchStatus("pass_single_mu", 0);
	old_tree->SetBranchStatus("pass_double_e", 0);
	old_tree->SetBranchStatus("pass_double_mu", 0);
	old_tree->SetBranchStatus("pass_elemu", 0);
	old_tree->SetBranchStatus("matchHLTPath", 0);
	old_tree->SetBranchStatus("mu0_mcMatchType", 0);
	old_tree->SetBranchStatus("mu1_mcMatchType", 0);
	old_tree->SetBranchStatus("ele0_mcMatchType", 0);
	old_tree->SetBranchStatus("ele1_mcMatchType", 0);
	old_tree->SetBranchStatus("tau0_mcMatchType", 0);
	old_tree->SetBranchStatus("tau1_mcMatchType", 0);
	
	//old_tree->SetBranchStatus("PU_weight");

	// create new tree and output file
	cout << "Output file created: " << output_file << endl;
	TFile* new_file = new TFile(output_file, "recreate");
	TTree* new_tree = old_tree->CloneTree();

	delete old_tree;

	new_tree->SetName("syncTree");
	new_file->Write();

	// event count
    cout << "# events with at least 1 preselected muon:" << "\t"
		 << new_tree->GetEntries("n_presel_mu>0") << endl;
	cout << "# events with at least 1 preselected electron: " << "\t"
		 << new_tree->GetEntries("n_presel_ele>0") << endl;
	cout << "# events with at least 1 preselected tau: " << "\t"
		 << new_tree->GetEntries("n_presel_tau>0") << endl;
	cout << "# events with at least 1 preselected jet: " << "\t"
		 << new_tree->GetEntries("n_presel_jet>0") << endl;
		
	delete new_file;
	delete old_file;
}
