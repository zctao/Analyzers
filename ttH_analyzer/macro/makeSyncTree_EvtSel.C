#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>
#include <vector>

using namespace std;

void makeSyncTree_EvtSel(
						 //TString dir="~/Documents/ttH/Outputs/80X/")
						 TString dir="/uscms/home/ztao/nobackup/")
{
	// open input file and read tree
	cout << "Opening root file from directory " << dir << endl;
	TFile* old_file1 = new TFile(dir+"output_sync_event_sr.root");
	TFile* old_file2 = new TFile(dir+"output_sync_event_fake.root");
	TFile* old_file3 = new TFile(dir+"output_sync_event_flip.root");
	
	TTree* old_tree1 = (TTree*)old_file1->Get("ttHtaus/eventTree");
	TTree* old_tree2 = (TTree*)old_file2->Get("ttHtaus/eventTree");
	TTree* old_tree3 = (TTree*)old_file3->Get("ttHtaus/eventTree");

	vector<TTree*> old_trees;
	old_trees.push_back(old_tree1);
	old_trees.push_back(old_tree2);
	old_trees.push_back(old_tree3);
	
	// turn off extra branches
	for (auto tree : old_trees) {
		tree->SetBranchStatus("event_weight", 0);
		tree->SetBranchStatus("HLT_Ele27_WPTight_Gsf", 0);
		tree->SetBranchStatus("HLT_IsoMu24", 0);
		tree->SetBranchStatus("HLT_IsoTkMu24", 0);
		tree->SetBranchStatus("HLT_IsoMu22_eta2p1", 0);
		tree->SetBranchStatus("HLT_IsoTkMu22_eta2p1", 0);
		tree->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 0);
		tree->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", 0);
		tree->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", 0);
		tree->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", 0);
		tree->SetBranchStatus("MC_weight_scale_muF0p5", 0);
		tree->SetBranchStatus("MC_weight_scale_muF2", 0);
		tree->SetBranchStatus("MC_weight_scale_muR0p5", 0);
		tree->SetBranchStatus("MC_weight_scale_muR2", 0);
		tree->SetBranchStatus("btagSF_weight_LFUp", 0);
		tree->SetBranchStatus("btagSF_weight_LFDown", 0);
		tree->SetBranchStatus("btagSF_weight_HFUp", 0);
		tree->SetBranchStatus("btagSF_weight_HFDown", 0);
		tree->SetBranchStatus("btagSF_weight_HFStats1Up", 0);
		tree->SetBranchStatus("btagSF_weight_HFStats1Down", 0);
		tree->SetBranchStatus("btagSF_weight_HFStats2Up", 0);
		tree->SetBranchStatus("btagSF_weight_HFStats2Down", 0);
		tree->SetBranchStatus("btagSF_weight_LFStats1Up", 0);
		tree->SetBranchStatus("btagSF_weight_LFStats1Down", 0);
		tree->SetBranchStatus("btagSF_weight_LFStats2Up", 0);
		tree->SetBranchStatus("btagSF_weight_LFStats2Down", 0);
		tree->SetBranchStatus("btagSF_weight_cErr1Up", 0);
		tree->SetBranchStatus("btagSF_weight_cErr1Down", 0);
		tree->SetBranchStatus("btagSF_weight_cErr2Up", 0);
		tree->SetBranchStatus("btagSF_weight_cErr2Down", 0);
		tree->SetBranchStatus("HiggsDecayType", 0);
		tree->SetBranchStatus("lepCategory", 0);
		tree->SetBranchStatus("btagCategory", 0);
		tree->SetBranchStatus("npuTrue", 0);
		tree->SetBranchStatus("npuInTime", 0);
		tree->SetBranchStatus("n_tau", 0);
		tree->SetBranchStatus("ibin", 0);
		tree->SetBranchStatus("tau0_decayMode", 0);
		tree->SetBranchStatus("tau1_decayMode", 0);
		tree->SetBranchStatus("pass_single_e", 0);
		tree->SetBranchStatus("pass_single_mu", 0);
		tree->SetBranchStatus("pass_double_e", 0);
		tree->SetBranchStatus("pass_double_mu", 0);
		tree->SetBranchStatus("pass_elemu", 0);
		tree->SetBranchStatus("matchHLTPath", 0);
		tree->SetBranchStatus("mu0_mcMatchType", 0);
		tree->SetBranchStatus("mu1_mcMatchType", 0);
		tree->SetBranchStatus("ele0_mcMatchType", 0);
		tree->SetBranchStatus("ele1_mcMatchType", 0);
		tree->SetBranchStatus("tau0_mcMatchType", 0);
		tree->SetBranchStatus("tau1_mcMatchType", 0);
	}
	
	// create new tree and output file
	TString output_file = "~/nobackup/ttHTT_syncNtuple/80X/syncNtuple_event.root";
	//TString output_file = "/afs/cern.ch/work/z/ztao/public/ttHTT_syncNtuple/80X/syncNtuple_event.root";
	cout << "Output file created: " << output_file << endl;
	TFile* new_file = new TFile(output_file, "recreate");

	vector<TTree*> new_trees;
	for (auto old_tree : old_trees) {
		TTree* new_tree = old_tree->CloneTree();
		new_trees.push_back(new_tree);
	}

	// delete old tree
	for (auto old_tree : old_trees) 
		delete old_tree;

	
	new_trees[0]->SetName("syncTree_2lSS1tau_SR");
	new_trees[1]->SetName("syncTree_2lSS1tau_Fake");
	new_trees[2]->SetName("syncTree_2lSS1tau_Flip");

	new_file->Write();

	// event count
	cout << "SR : " << new_trees[0]->GetEntries() << endl;
	cout << "Fake : " << new_trees[1]->GetEntries() << endl;
	cout << "Flip : " << new_trees[2]->GetEntries() << endl;
		
	delete new_file;
	delete old_file1;
	delete old_file2;
	delete old_file3;
}
