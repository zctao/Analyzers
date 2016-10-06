#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>

using namespace std;

void makeSyncTree(TString input_file)
{
	// open input file and read tree
	TFile* old_file = new TFile(input_file);
	TTree* old_tree = (TTree*)old_file->Get("ttHtaus/eventTree");

	// turn off extra branches
	old_tree->SetBranchStatus("event_weight",0);
	old_tree->SetBranchStatus("genWeight", 0);
	old_tree->SetBranchStatus("n_tau", 0);
	old_tree->SetBranchStatus("tau0_decayMode", 0);
	old_tree->SetBranchStatus("tau1_decayMode", 0);

	// create new tree and output file
	TFile* new_file = new TFile("syncNtuple_ttH_80X.root", "recreate");
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
