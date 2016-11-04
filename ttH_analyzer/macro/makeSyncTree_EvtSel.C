#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>
#include <vector>

using namespace std;

void makeSyncTree_EvtSel(
						 //TString dir="~/Documents/ttH/Outputs/80X")
						 TString dir="/uscms/home/ztao/nobackup")
{
	// open input file and read tree
	TFile* old_file1 = new TFile(dir+"/output_sync_event_sr.root");
	TFile* old_file2 = new TFile(dir+"/output_sync_event_fake.root");
	TFile* old_file3 = new TFile(dir+"/output_sync_event_flip.root");
	
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
		tree->SetBranchStatus("MC_weight_scale_muF0p5", 0);
		tree->SetBranchStatus("MC_weight_scale_muF2", 0);
		tree->SetBranchStatus("MC_weight_scale_muR0p5", 0);
		tree->SetBranchStatus("MC_weight_scale_muR2", 0);
		tree->SetBranchStatus("n_tau", 0);
		tree->SetBranchStatus("tau0_decayMode", 0);
		tree->SetBranchStatus("tau1_decayMode", 0);
		tree->SetBranchStatus("pass_single_e", 0);
		tree->SetBranchStatus("pass_single_mu", 0);
		tree->SetBranchStatus("pass_double_e", 0);
		tree->SetBranchStatus("pass_double_mu", 0);
		tree->SetBranchStatus("pass_elemu", 0);

		tree->SetBranchStatus("PU_weight", 0);
	}
	
	// create new tree and output file
	//TFile* new_file = new TFile("/afs/cern.ch/work/z/ztao/public/ttHTT_syncNtuple/80X/syncNtuple_event.root", "recreate");
	TFile* new_file = new TFile("~/nobackup/ttHTT_syncNtuple/80X/syncNtuple_event.root", "recreate");

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
