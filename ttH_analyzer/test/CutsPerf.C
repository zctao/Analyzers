#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TPad.h"
#include "TMath.h"
#include "TPaveStats.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm> 

using namespace std;


const int nstep = 20;

void getEffArray(double ptmin, double ptmax, int nstep, double eff[4][20], TTree* tree);
void MakeROCPlot(TTree *tree_sig, int nevt_sig, TTree* tree_bkg, int nevt_bkg);
int HistFiller(TTree* tree, TString label);
void CutHistDrawer(TString histfile1, TString histfile2);
void NtupleHistDrawer(TH1F* h_cutflow_sig, TH1F* h_cutflow_TTJets,
				   TH1F* h_njets_sig, TH1F* h_njets_TTJets,
				   TH1F* h_nbtags_sig, TH1F* h_nbtags_TTJets,
				   TH1F* h_ntauID_sig, TH1F* h_ntauID_TTJets);

int nsig = -99;
int nTTJets = -99;
int nsample_sig = -999;
int nsample_TTJets = -999;

//scale factor
double sf_sig = -99.9;
double sf_TTJets = -99.9;

const int sig_scale = 10;

void CutsPerf(
			  //const TString sig_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_sig.root",
			  //const TString bkg_file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_TTJets.root"
			  const TString sig_file = "/eos/uscms/store/user/ztao/ttHToTT_M125_13TeV_powheg_pythia8/ttHToTauTau_Ntuple_signal/151119_151045/0000/CU_ttH_EDA_output_tmp.root",
			  const TString bkg_file = "/eos/uscms/store/user/ztao/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ttHToTauTau_Ntuple_TTJets/151118_223202/0000/CU_ttH_EDA_output_tmp.root"
							)
{
	// Read ntuples
	TFile* f_sig = new TFile(sig_file);
	TFile* f_TTJets = new TFile(bkg_file);
	
	TTree* tree_sig = (TTree*) f_sig->Get("ttHtautau/EventTree");
	TTree* tree_TTJets = (TTree*) f_TTJets->Get("ttHtautau/EventTree");

	if (tree_sig->GetEntries() == 0 or tree_TTJets->GetEntries() == 0) {
		cout << "File doesn't exist or is empty, returning..." << endl;
		cout << endl;
		return;
	}
	else {
		cout << "Number of entries (signal):" << tree_sig->GetEntries() << endl;
		cout << "Number of entries (TTJets) :" << tree_TTJets->GetEntries() << endl;
	}

	// Total number of events
	// (before basic selection, not the number of events in the tree)
	// Get the total number of events from first bin of cut flow histogram

	TH1F* h_cutflow_sig =
		(TH1F*) f_sig->Get("ttHtautau/h_tth_syncex_dileptauh");
	TH1F* h_cutflow_TTJets =
		(TH1F*) f_TTJets->Get("ttHtautau/h_tth_syncex_dileptauh");
	
	nsample_sig = h_cutflow_sig -> GetBinContent(1);
	nsample_TTJets = h_cutflow_TTJets -> GetBinContent(1);

	cout << "Total number of signal sample processed :" << nsample_sig << endl;
	cout << "Total number of TTJets sample processed :" << nsample_TTJets << endl;

	// Calculate scale factors
	// assume 100 fb^-1
	sf_sig = 0.5085 * 0.0637 * 100000 / nsample_sig;
	sf_TTJets = 832 * 100000 / nsample_TTJets;
	
	// Get jet multiplicity histograms and tauID histograms
	TH1F* h_njets_sig =
		(TH1F*) f_sig->Get("ttHtautau/h_njets");
	TH1F* h_nbtags_sig =
		(TH1F*) f_sig->Get("ttHtautau/h_nbtags");
	TH1F* h_ntauID_sig =
		(TH1F*) f_sig->Get("ttHtautau/h_ntauID");
	TH1F* h_njets_TTJets =
		(TH1F*) f_TTJets->Get("ttHtautau/h_njets");
	TH1F* h_nbtags_TTJets =
		(TH1F*) f_TTJets->Get("ttHtautau/h_nbtags");
	TH1F* h_ntauID_TTJets =
		(TH1F*) f_TTJets->Get("ttHtautau/h_ntauID");

	//MakeROCPlot(tree_sig, nevt_sig, tree_TTJets, nevt_TTJets);

	nsig = HistFiller(tree_sig, "sig");
	nTTJets = HistFiller(tree_TTJets, "TTJets");

	std::cout << "number of events after cuts (signal) :" << nsig << endl;
	std::cout <<  "number of events after cuts (TTJets) :" << nTTJets << endl;
	
	// Draw Histograms
	CutHistDrawer("/uscms/home/ztao/work/CU_ttH_WD/Outputs/cuts_hist_sig.root", "/uscms/home/ztao/work/CU_ttH_WD/Outputs/cuts_hist_TTJets.root");
	
	NtupleHistDrawer(h_cutflow_sig, h_cutflow_TTJets,
				  h_njets_sig,   h_njets_TTJets,
				  h_nbtags_sig,  h_nbtags_TTJets,
				  h_ntauID_sig,  h_ntauID_TTJets);
	
}


// ------------------------------------------------
	
void getEffArray(double ptmin /*=20*/, double ptmax, int nstep, double eff[4][20], TTree* tree)
// ** NEED TO BE UPDATED **
{
	// initialize arrays
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 20; ++j)
			eff[i][j] = 0;
	}
	
	const int nEntries = tree->GetEntries();
	
	// Define leafs and branches

	std::vector<float>* loose_tau_pt;
	std::vector<float>* loose_tau_eta;
	std::vector<float>* medium_tau_pt;
	std::vector<float>* medium_tau_eta;
	std::vector<float>* tight_tau_pt;
	std::vector<float>* tight_tau_eta;
	
	TBranch* b_loose_tau_pt;
	TBranch* b_loose_tau_eta;
	TBranch* b_medium_tau_pt;
	TBranch* b_medium_tau_eta;
	TBranch* b_tight_tau_pt;
	TBranch* b_tight_tau_eta;

	loose_tau_pt = 0;
	loose_tau_eta = 0;
	medium_tau_pt = 0;
	medium_tau_eta = 0;
	tight_tau_pt = 0;
	tight_tau_eta = 0;

	tree->SetBranchAddress("loose_tau_pt", &loose_tau_pt, &b_loose_tau_pt);
	tree->SetBranchAddress("loose_tau_eta", &loose_tau_eta, &b_loose_tau_eta);
	tree->SetBranchAddress("medium_tau_pt", &medium_tau_pt, &b_medium_tau_pt);
	tree->SetBranchAddress("medium_tau_eta", &medium_tau_eta, &b_medium_tau_eta);
	tree->SetBranchAddress("tight_tau_pt", &tight_tau_pt, &b_tight_tau_pt);
	tree->SetBranchAddress("tight_tau_eta", &tight_tau_eta, &b_tight_tau_eta);

	// ----------------------------------------------------------------------------------------------------------------
	//        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
	// ----------------------------------------------------------------------------------------------------------------
	
	// event loop
	for (int ievt=0; ievt<nEntries; ++ievt) {
		tree -> GetEntry(ievt);

		for (int i = 0; i<nstep; ++i) {
			double ptcut = ptmin+(ptmax-ptmin)/nstep*i;

			// tight
			for (auto & taupt : *tight_tau_pt) {
				if (taupt < ptcut)
					continue;
				
				eff[3][i]++;
				break;  // find at least one tau that passes pt cut
				
			}

			// medium
			for (auto & taupt : *medium_tau_pt) {
				if (taupt < ptcut)
					continue;
				
				eff[2][i]++;
				break;  // find at least one tau that passes pt cut
				
			}

			// loose
			for (auto & taupt : *loose_tau_pt) {
				if (taupt < ptcut)
					continue;
				
				eff[1][i]++;
				break;  // find at least one tau that passes pt cut
				
			}
			
		} // end of pt_cut loop

		
	} // end of event loop
	
	delete tree;
	
}

void MakeROCPlot(TTree *tree_sig, int nevt_sig, TTree* tree_bkg, int nevt_bkg)
{ // FIX ME
	
	// Define efficiency arrays
	//const int nstep = 20;
	double eff_sig[4][20];
	double eff_bkg[4][20];
	
	getEffArray(20, 120, nstep, eff_sig, tree_sig);
	getEffArray(20, 120, nstep, eff_bkg, tree_bkg);
	
	
	
	for (int i = 1; i < 4; ++i) {
		for (int j = 0; j < nstep; ++j){
			eff_sig[i][j]/=nevt_sig;
		}
	}

	for (int i = 1; i < 4; ++i) {
		for (int j = 0; j < nstep; ++j){
			eff_bkg[i][j]/=nevt_bkg;
		}
	}

	TCanvas c;
	gPad->SetGrid();
	gPad->SetLogy();

	//TGraphErrors* gr_noniso =
	//new TGraphErrors(nstep, eff_sig[0], eff_bkg[0], 0, 0);
	TGraphErrors* gr_loose =
		new TGraphErrors(nstep, eff_sig[1], eff_bkg[1], 0, 0);
	TGraphErrors* gr_medium =
		new TGraphErrors(nstep, eff_sig[2], eff_bkg[2], 0, 0);
	TGraphErrors* gr_tight =
		new TGraphErrors(nstep, eff_sig[3], eff_bkg[3], 0, 0);

	//gr_noniso->SetLineColor(4);
	//gr_noniso->SetMarkerColor(4);
	gr_loose->SetLineColor(3);
    gr_loose->SetMarkerColor(3);
	gr_medium->SetLineColor(2);
    gr_medium->SetMarkerColor(2);
	gr_tight->SetLineColor(1);
    gr_tight->SetMarkerColor(1);
	
	gr_loose->SetTitle("Tau ID WP");
	gr_loose->GetXaxis()->SetTitle("Signal efficiency");
	gr_loose->GetYaxis()->SetTitle("Background efficiency");
	gr_loose->SetMarkerStyle(21);
	gr_medium->SetMarkerStyle(20);
	gr_tight->SetMarkerStyle(22);
	//gr_noniso->SetMarkerStyle(7);

	gr_loose->Draw("APZ0");
	gr_medium->Draw("same PZ0");
	gr_tight->Draw("same PZ0");
	//gr_noniso->Draw("same PZ0");
	
	TLegend* leg=new TLegend(0.23,0.7,0.50,0.9);
	//leg->AddEntry(gr_noniso,"nonIso","p");
	leg->AddEntry(gr_loose,"loose","p");
	leg->AddEntry(gr_medium,"medium","p");
	leg->AddEntry(gr_tight,"tight","p");
	leg->Draw("same");

	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/tauWPROC.pdf");
}

int HistFiller(TTree* tree, TString label) 
{
	int nevtpass = 0;

	int bug_cnt = 0;
	
	const int nEntries = tree->GetEntries();
	cout << "nEntries :" << nEntries << endl;
	
	// Define leafs and branches
	int n_electrons = -99;
	int n_muons = -99;
	int n_loose_taus = -99;
	int n_medium_taus = -99;
	int n_tight_taus = -99;
	int n_jets = -99;
	int n_btags = -99;
	float pv_x = -99;
	float pv_y = -99;
	float pv_z = -99;
	std::vector<float>* e_pt = 0;
	std::vector<float>* e_eta = 0;
	std::vector<float>* e_phi = 0;
	std::vector<float>* e_mass = 0;
	std::vector<int>* e_charges = 0;
	std::vector<float>* e_vtx_dz = 0;
	std::vector<float>* e_vtx_dxy = 0;
	std::vector<bool>* e_isGsfCtfScPixChargeConsistent = 0;
	std::vector<float>* mu_pt = 0;
	std::vector<float>* mu_eta = 0;
	std::vector<float>* mu_phi = 0;
	std::vector<float>* mu_mass = 0;
	std::vector<int>* mu_charges = 0;
	std::vector<float>* mu_vtx_dz = 0;
	std::vector<float>* mu_vtx_dxy = 0;
	std::vector<float>* mu_relTrkPtError = 0;
	std::vector<float>* loose_tau_pt = 0;
	std::vector<float>* loose_tau_eta = 0;
	std::vector<float>* loose_tau_phi = 0;
	std::vector<float>* loose_tau_mass = 0;
	std::vector<int>* loose_tau_charges = 0;
	std::vector<size_t>* loose_tau_decaymode = 0;
	std::vector<float>* jet_pt = 0;
	std::vector<float>* jet_eta = 0;
	std::vector<float>* jet_phi = 0;
	std::vector<float>* jet_mass = 0;
	std::vector<float>* jet_charges = 0;
	std::vector<float>* bjet_pt = 0;
	std::vector<float>* bjet_eta = 0;
	std::vector<float>* bjet_phi = 0;
	std::vector<float>* bjet_mass = 0;
	std::vector<float>* bjet_charges = 0;
	float met_x = -9999.9;
	float met_y = -9999.9;

	TBranch* b_n_electrons;
	TBranch* b_n_muons;
	TBranch* b_n_loose_taus;
	TBranch* b_n_medium_taus;
	TBranch* b_n_tight_taus;
	TBranch* b_n_jets;
	TBranch* b_n_btags;
	TBranch* b_pv_x;
	TBranch* b_pv_y;
	TBranch* b_pv_z;
	TBranch* b_e_pt;
	TBranch* b_e_eta;
	TBranch* b_e_phi;
	TBranch* b_e_mass;
	TBranch* b_e_charges;
	TBranch* b_e_vtx_dz;
	TBranch* b_e_vtx_dxy;
	TBranch* b_e_isGsfCtfScPixChargeConsistent;
	TBranch* b_mu_pt;
	TBranch* b_mu_eta;
	TBranch* b_mu_phi;
	TBranch* b_mu_mass;
	TBranch* b_mu_charges;
	TBranch* b_mu_vtx_dz;
	TBranch* b_mu_vtx_dxy;
	TBranch* b_mu_relTrkPtError;
	TBranch* b_loose_tau_pt;
	TBranch* b_loose_tau_eta;
	TBranch* b_loose_tau_phi;
	TBranch* b_loose_tau_mass;
	TBranch* b_loose_tau_charges;
	TBranch* b_loose_tau_decaymode;
	TBranch* b_jet_pt;
	TBranch* b_jet_eta;
	TBranch* b_jet_phi;
	TBranch* b_jet_mass;
	TBranch* b_jet_charges;
	TBranch* b_bjet_pt;
	TBranch* b_bjet_eta;
	TBranch* b_bjet_phi;
	TBranch* b_bjet_mass;
	TBranch* b_bjet_charges;
	TBranch* b_met_x;
	TBranch* b_met_y;

	tree->SetBranchAddress("n_electrons", &n_electrons, &b_n_electrons);
	tree->SetBranchAddress("n_loose_taus", &n_loose_taus, &b_n_loose_taus);
	tree->SetBranchAddress("n_medium_taus", &n_medium_taus, &b_n_medium_taus);
	tree->SetBranchAddress("n_tight_taus", &n_tight_taus, &b_n_tight_taus);
	tree->SetBranchAddress("n_jets", &n_jets, &b_n_jets);
	tree->SetBranchAddress("n_btags", &n_btags, &b_n_btags);
	tree->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
	tree->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
	tree->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
	tree->SetBranchAddress("e_pt", &e_pt, &b_e_pt);
	tree->SetBranchAddress("e_eta", &e_eta, &b_e_eta);
	tree->SetBranchAddress("e_phi", &e_phi, &b_e_phi);
	tree->SetBranchAddress("e_mass", &e_mass, &b_e_mass);
	tree->SetBranchAddress("e_charges", &e_charges, &b_e_charges);
	tree->SetBranchAddress("e_vtx_dz", &e_vtx_dz, &b_e_vtx_dz);
	tree->SetBranchAddress("e_vtx_dxy", &e_vtx_dxy, &b_e_vtx_dxy);
	tree->SetBranchAddress("e_isGsfCtfScPixChargeConsistent", &e_isGsfCtfScPixChargeConsistent, &b_e_isGsfCtfScPixChargeConsistent);
	tree->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
	tree->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
	tree->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
	tree->SetBranchAddress("mu_mass", &mu_mass, &b_mu_mass);
	tree->SetBranchAddress("mu_charges", &mu_charges, &b_mu_charges);
	tree->SetBranchAddress("mu_vtx_dz", &mu_vtx_dz, &b_mu_vtx_dz);
	tree->SetBranchAddress("mu_vtx_dxy", &mu_vtx_dxy, &b_mu_vtx_dxy);
	tree->SetBranchAddress("mu_relTrkPtError", &mu_relTrkPtError, &b_mu_relTrkPtError);
	tree->SetBranchAddress("loose_tau_pt", &loose_tau_pt, &b_loose_tau_pt);
	tree->SetBranchAddress("loose_tau_eta", &loose_tau_eta, &b_loose_tau_eta);
	tree->SetBranchAddress("loose_tau_phi", &loose_tau_phi, &b_loose_tau_phi);
	tree->SetBranchAddress("loose_tau_mass", &loose_tau_mass, &b_loose_tau_mass);
	tree->SetBranchAddress("loose_tau_charges", &loose_tau_charges, &b_loose_tau_charges);
	tree->SetBranchAddress("loose_tau_decaymode", &loose_tau_decaymode, &b_loose_tau_decaymode);
	tree->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
	tree->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
	tree->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
	tree->SetBranchAddress("jet_mass", &jet_mass, &b_jet_mass);
	tree->SetBranchAddress("jet_charges", &jet_charges, &b_jet_charges);
	tree->SetBranchAddress("bjet_pt", &bjet_pt, &b_bjet_pt);
	tree->SetBranchAddress("bjet_eta", &bjet_eta, &b_bjet_eta);
	tree->SetBranchAddress("bjet_phi", &bjet_phi, &b_bjet_phi);
	tree->SetBranchAddress("bjet_mass", &bjet_mass, &b_bjet_mass);
	tree->SetBranchAddress("bjet_charges", &bjet_charges, &b_bjet_charges);
	tree->SetBranchAddress("MET_x", &met_x, &b_met_x);
	tree->SetBranchAddress("MET_y", &met_y, &b_met_y);

	// Histograms
	TH1F* h_lep1pt = new TH1F("h_lep1pt", "", 30, 0., 300.);
	TH1F* h_lep1eta = new TH1F("h_lep1eta", "", 25, -2.5, 2.5);
	TH1F* h_lep2pt = new TH1F("h_lep2pt", "", 30, 0., 300.);
	TH1F* h_lep2eta = new TH1F("h_lep2eta", "", 25, -2.5, 2.5);
	TH1F* h_taupt = new TH1F("h_taupt", "", 30, 0., 300.);
	TH1F* h_taueta = new TH1F("h_taueta", "", 25, -2.5, 2.5);
	TH1F* h_bjet1pt = new TH1F("h_bjet1pt", "", 50, 0., 500.);
	TH1F* h_bjet1eta = new TH1F("h_bjet1eta", "", 25, -2.5, 2.5);
	TH1F* h_dRlep1tau = new TH1F("h_dRlep1tau", "", 25, 0., 5.);
	TH1F* h_dRlep2tau = new TH1F("h_dRlep2tau", "", 25, 0., 5.);
	TH1F* h_lep1dxy = new TH1F("h_lep1dxy", "", 20, 0, 0.02);
	TH1F* h_lep1dz = new TH1F("h_lep1dz", "", 20, 0., 0.05);
	TH1F* h_lep2dxy = new TH1F("h_lep2dxy", "", 20, 0, 0.02);
	TH1F* h_lep2dz = new TH1F("h_lep2dz", "", 20, 0., 0.05);
	//TH1F* h_tau1dz = new TH1F("h_tau1dz", "", 20, 0., 0.005);
	//TH1F* h_tau1dxy = new TH1F("h_tau1dxy", "", 20, 0., 0.005);
	TH1F* h_lepmaxdxy = new TH1F("h_lepmaxdxy", "", 20, 0, 0.05);
	TH1F* h_visM_taulep1 = new TH1F("h_visM_taulep1", "", 30, 0, 300);
	TH1F* h_visM_taulep2 = new TH1F("h_visM_taulep2", "", 30, 0, 300);
	TH1F* h_visM_taulepmaxdxy = new TH1F("h_visM_taulepmaxdxy", "", 30, 0, 300);
	TH1F* h_visM_taulepmindR = new TH1F("h_visM_taulepmindR", "", 30, 0, 300);
	
	
	// ----------------------------------------------------------------------------------------------------------------
	//        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
	// ----------------------------------------------------------------------------------------------------------------
	
	// event loop
	for (int ievt=0; ievt<nEntries; ++ievt) {

		tree -> GetEntry(ievt);
			
		// combine lepton vectors
		std::vector<float> *lep_pt, *lep_eta, *lep_phi, *lep_mass;
		std::vector<float> *lep_vtx_dz, *lep_vtx_dxy;
		std::vector<int> *lep_charges;

		/*
		////////////////////////////////////
		// debug
		int lepcomb = 10*e_pt->size()+mu_pt->size();
		size_t esize = e_pt->size();
		size_t musize = mu_pt->size();
		size_t lepsize = e_pt->size()+mu_pt->size();

		// check how many events have >=2 leptons with |dz|<cuts(0.01)
		int vzcutpass = 0;
		for (size_t ive=0; ive<e_vz->size(); ++ive) {
			if (abs(e_vz->at(ive)) <= 0.01)
				vzcutpass++;
		}
		for (size_t ivm=0; ivm<mu_vz->size(); ++ivm) {
			if (abs(mu_vz->at(ivm)) <= 0.01)
				vzcutpass++;
		}

		if (vzcutpass < 2)
			continue;
		
		////////////////////////////////////
		*/
		
		lep_pt = e_pt;
		lep_eta = e_eta;
		lep_phi = e_phi;
		lep_mass = e_mass;
		lep_charges = e_charges;
		lep_vtx_dz = e_vtx_dz;
		lep_vtx_dxy = e_vtx_dxy;
		lep_pt->reserve(e_pt->size()+mu_pt->size());
		lep_eta->reserve(e_eta->size()+mu_eta->size());
		lep_phi->reserve(e_phi->size()+mu_phi->size());
		lep_mass->reserve(e_mass->size()+mu_mass->size());
		lep_charges->reserve(e_charges->size()+mu_charges->size());
		lep_vtx_dz->reserve(e_vtx_dz->size()+mu_vtx_dz->size());
		lep_vtx_dxy->reserve(e_vtx_dxy->size()+mu_vtx_dxy->size());
		lep_pt->insert(lep_pt->end(), mu_pt->begin(), mu_pt->end());
		lep_eta->insert(lep_eta->end(), mu_eta->begin(), mu_eta->end());
		lep_phi->insert(lep_phi->end(), mu_phi->begin(), mu_phi->end());
		lep_mass->insert(lep_mass->end(), mu_mass->begin(), mu_mass->end());
		lep_charges->insert(lep_charges->end(), mu_charges->begin(), mu_charges->end());
		lep_vtx_dz->insert(lep_vtx_dz->end(), mu_vtx_dz->begin(), mu_vtx_dz->end());
		lep_vtx_dxy->insert(lep_vtx_dxy->end(), mu_vtx_dxy->begin(), mu_vtx_dxy->end());
				
		
		// cut on dz
		/* !!!NotRight
		for (size_t ilep=0; ilep<lep_vz->size(); ++ilep) {
			if (fabs(lep_vz->at(ilep))>0.01) {
				lep_pt->erase(lep_pt->begin()+ilep);
				lep_eta->erase(lep_eta->begin()+ilep);
				lep_phi->erase(lep_phi->begin()+ilep);
				lep_mass->erase(lep_mass->begin()+ilep);
				lep_charges->erase(lep_charges->begin()+ilep);
				lep_vx->erase(lep_vx->begin()+ilep);
				lep_vy->erase(lep_vy->begin()+ilep);
				lep_vz->erase(lep_vz->begin()+ilep);
			}
		}
		
		if (lep_vz->size()<2)
			continue;
		
		else if (lep_pt->size()==2) {
			if (lep_charges->at(0) * lep_charges->at(1) < 0)
				continue;
		}
		*/
		/*
		// debug
		// check how many events have >=2 leptons with |dz|<cuts(0.01)
		int vzcutpass = 0;
		for (size_t ilep = 0; ilep<lep_vz->size(); ++ilep) {
			if (abs(lep_vz->at(ilep)) <= 0.01)
				vzcutpass++;
		}
		if (vzcutpass < 2)
			continue;
		
		////////////////////////////////////
		// debug
		int netlepq = 0;
		for (size_t il=0; il<lep_charges->size(); ++il) {
			netlepq += lep_charges->at(il);
		}
		int nettauq = 0;
		for (size_t itau=0; itau<loose_tau_charges->size(); ++itau) {
			nettauq += loose_tau_charges->at(itau);
		}
		if (netlepq * nettauq >=0) {
			++bug_cnt;
			std::cout << "lep_pt size :" << lep_pt->size() << endl;
			std::cout << "e+mu size :" << lepsize << endl;
			std::cout << "leps net charge :" << netlepq << endl;
			std::cout << "lep charges :" << lep_charges->at(0) << lep_charges->at(1) << endl;
			std::cout << "lepcomb :" << lepcomb << endl;
			std::cout << "e size :" << esize << "\t" << "mu size :" << musize << endl;
			std::cout << "tau net charge :" << nettauq << endl;
			std::cout << "loose_tau size :" << loose_tau_charges->size() << endl;
			std::cout << "bug_cnt :" << bug_cnt << endl;
			std::cout << "lep vertices :" << endl;
			for (size_t ilep=0; ilep<lep_vz->size(); ++ilep) {
				cout << "dz :" << lep_vz->at(ilep) - pv_z << endl;
			}
		}
		////////////////////////////////////
		*/

		// Look for leading and sub-leading leptons
		// require both lepton charges are of opposite sign to tau charge

		float ptmax = 0, ptsubmax = 0;
		size_t ilep1=-99, ilep2=-99;
		//size_t ilep1=0, ilep2=0;
		
		for (size_t itr=0; itr < lep_charges->size(); ++itr) {
			if (lep_charges->at(itr) * loose_tau_charges->at(0) > 0)
				continue;
			if (lep_pt->at(itr)>ptmax) {
				ptsubmax = ptmax;
				ptmax = lep_pt->at(itr);
				ilep2 = ilep1;
				ilep1 = itr;
			}
			else if (lep_pt->at(itr)>ptsubmax) {
				ptsubmax = lep_pt->at(itr);
				ilep2 = itr;
			}
		}
		
		
		/*
		if (lep_charges->at(0) * tau_charge > 0 and lep_pt->size() < 3) {
			//std::cout << "initialize ilep1" << endl;
		
		}
		
		for (size_t i=0; i<lep_pt->size(); ++i) {
			if (lep_charges->at(i) * tau_charge  < 0) {
				ilep1 = i;
				break;
			}
		}
		//std::cout << "ilep1 :" << ilep1 << endl;
		if (ilep1 >= lep_pt->size()) {
			std::cout << "Oooooops" << endl;
			break;
		}
		
		for (size_t i=0; i<lep_pt->size(); ++i) {
			if (lep_charges->at(i)* tau_charge > 0 )
				continue;
			if (lep_pt->at(i) < lep_pt->at(ilep1))
				continue;
			ilep1 = i;			
		}
		
		//std::cout << "initialize ilep2" << endl;
		for (size_t i=0; i<lep_pt->size(); ++i) {
			if (i != ilep1 and lep_charges->at(i) * tau_charge < 0) {
				ilep2 = i;
				break;
			}
				
		}
		if (ilep2 >= lep_pt->size()) {
			std::cout << "Oooooops2" << endl;
			break;
		}
		
		for (size_t i=0; i<lep_pt->size(); ++i) {
			if (lep_charges->at(i) * tau_charge > 0 ) continue;
			if (i==ilep1) continue;
			if (lep_pt->at(i) < lep_pt->at(ilep2)) continue;
			ilep2 = i;
		}
		
		////////////////////////////////////
		// debug
		bool sort_bug = false;
		if (lep_pt->at(ilep1)<lep_pt->at(ilep2)) sort_bug = true;
		for (size_t i = 0; i< lep_pt->size(); ++i) {
			if (lep_charges->at(i) * tau_charge > 0) continue;
			if (i == ilep1 or i == ilep2) continue;
			if (lep_pt->at(i)>lep_pt->at(ilep2)) sort_bug = true;
		}
		if (sort_bug) std::cout << "sorting oops" << endl;
		////////////////////////////////////
		*/
		
		float dz_lep1 = lep_vtx_dz->at(ilep1);
		float dz_lep2 = lep_vtx_dz->at(ilep2);
		float dxy_lep1 = lep_vtx_dxy->at(ilep1);
		float dxy_lep2 = lep_vtx_dxy->at(ilep2);
		//float dxy_tau1 = TMath::Sqrt((loose_tau_vx->at(itau1)-pv_x)*(loose_tau_vx->at(itau1)-pv_x)+(loose_tau_vy->at(itau1)-pv_y)*(loose_tau_vy->at(itau1)-pv_y));

		
		++nevtpass;

		//std::cout << "fill histograms" << endl;
		/// Fill histograms
		h_lep1pt -> Fill(lep_pt->at(ilep1));
		h_lep1eta -> Fill(lep_eta->at(ilep1));
		h_lep2pt -> Fill(lep_pt->at(ilep2));
		h_lep2eta -> Fill(lep_eta->at(ilep2));

		h_taupt -> Fill(loose_tau_pt->at(0));
		h_taueta -> Fill(loose_tau_eta->at(0));

		h_bjet1pt -> Fill(bjet_pt->at(0));
		h_bjet1eta -> Fill(bjet_eta->at(0));

		float deta_l1T = lep_eta->at(ilep1)-loose_tau_eta->at(0);
		float dphi_l1T = TVector2::Phi_mpi_pi(lep_phi->at(ilep1)-loose_tau_phi->at(0));
		float dRlep1tau = TMath::Sqrt(deta_l1T*deta_l1T+dphi_l1T*dphi_l1T);
		
		float deta_l2T = lep_eta->at(ilep2)-loose_tau_eta->at(0);
		float dphi_l2T = TVector2::Phi_mpi_pi(lep_phi->at(ilep2)-loose_tau_phi->at(0));
		float dRlep2tau = TMath::Sqrt(deta_l2T*deta_l2T+dphi_l2T*dphi_l2T);

		h_dRlep1tau -> Fill(dRlep1tau);
		h_dRlep2tau -> Fill(dRlep2tau);

		// vertex
		h_lep1dxy -> Fill(lep_vtx_dxy->at(ilep1));
		h_lep1dz -> Fill(lep_vtx_dz->at(ilep1));
		h_lep2dxy -> Fill(lep_vtx_dxy->at(ilep2));
		h_lep2dz -> Fill(lep_vtx_dz->at(ilep2));
		//h_tau1dz -> Fill(dz_tau1);
		//h_tau1dxy -> Fill(dxy_tau1);
		h_lepmaxdxy -> Fill(std::max(lep_vtx_dxy->at(ilep1),lep_vtx_dxy->at(ilep2)));

		// visible mass
		TLorentzVector tau_p4, lep1_p4, lep2_p4, lep_maxdxy_p4,lep_mindR_p4;
		tau_p4.SetPtEtaPhiM(loose_tau_pt->at(0),loose_tau_eta->at(0),
							loose_tau_phi->at(0),loose_tau_mass->at(0));
		lep1_p4.SetPtEtaPhiM(lep_pt->at(ilep1), lep_eta->at(ilep1),
							 lep_phi->at(ilep1), lep_mass->at(ilep1));
		lep2_p4.SetPtEtaPhiM(lep_pt->at(ilep2), lep_eta->at(ilep2),
							 lep_phi->at(ilep2), lep_mass->at(ilep2));

		size_t ilep_maxdxy, ilep_mindR;
		if (lep_vtx_dxy->at(ilep1) > lep_vtx_dxy->at(ilep2))
			ilep_maxdxy = ilep1;
		else
			ilep_maxdxy = ilep2;
		
		if (dRlep1tau < dRlep2tau)
			ilep_mindR = ilep1;
		else
			ilep_mindR = ilep2;

		lep_maxdxy_p4.SetPtEtaPhiM(lep_pt->at(ilep_maxdxy),lep_eta->at(ilep_maxdxy),
							 lep_phi->at(ilep_maxdxy), lep_mass->at(ilep_maxdxy));
		lep_mindR_p4.SetPtEtaPhiM(lep_pt->at(ilep_mindR),lep_eta->at(ilep_mindR),
								  lep_phi->at(ilep_mindR),lep_mass->at(ilep_mindR));
		h_visM_taulep1 -> Fill((tau_p4+lep1_p4).M());
		h_visM_taulep2 -> Fill((tau_p4+lep2_p4).M());
		h_visM_taulepmaxdxy -> Fill((tau_p4+lep_maxdxy_p4).M());
		h_visM_taulepmindR -> Fill((tau_p4+lep_mindR_p4).M());
		
	} // end of event loop

	TFile *outputfile = new TFile(
				     "/uscms/home/ztao/work/CU_ttH_WD/Outputs/cuts_hist_"+label+".root",
					 "RECREATE");
	
	h_lep1pt -> Write();
	h_lep1eta -> Write();
	h_lep2pt -> Write();
	h_lep2eta -> Write();
	h_taupt -> Write();
	h_taueta -> Write();
	h_bjet1pt -> Write();
	h_bjet1eta -> Write();
	h_dRlep1tau -> Write();
	h_dRlep2tau -> Write();
	h_lep1dxy -> Write();
	h_lep1dz -> Write();
	h_lep2dxy -> Write();
	h_lep2dz -> Write();
	//h_tau1dz -> Write();
	//h_tau1dxy -> Write();
	h_lepmaxdxy -> Write();
	h_visM_taulep1 -> Write();
	h_visM_taulep2 -> Write();
	h_visM_taulepmaxdxy -> Write();
	h_visM_taulepmindR -> Write();
	
	delete outputfile;

	return nevtpass;
}


void CutHistDrawer(TString histfile1, TString histfile2)
{

	TFile* f_sig = new TFile(histfile1);
	TFile* f_TTJets = new TFile(histfile2);

	TH1F* h_lep1pt_sig = (TH1F*)f_sig->Get("h_lep1pt");
	TH1F* h_lep1eta_sig = (TH1F*)f_sig->Get("h_lep1eta");
	TH1F* h_lep2pt_sig = (TH1F*)f_sig->Get("h_lep2pt");
	TH1F* h_lep2eta_sig = (TH1F*)f_sig->Get("h_lep2eta");
	TH1F* h_taupt_sig = (TH1F*)f_sig->Get("h_taupt");
	TH1F* h_taueta_sig = (TH1F*)f_sig->Get("h_taueta");
	TH1F* h_bjet1pt_sig = (TH1F*)f_sig->Get("h_bjet1pt");
	TH1F* h_bjet1eta_sig = (TH1F*)f_sig->Get("h_bjet1eta");
	TH1F* h_dRlep1tau_sig = (TH1F*)f_sig->Get("h_dRlep1tau");
	TH1F* h_dRlep2tau_sig = (TH1F*)f_sig->Get("h_dRlep2tau");
	TH1F* h_lep1dxy_sig = (TH1F*)f_sig->Get("h_lep1dxy");
	TH1F* h_lep2dxy_sig = (TH1F*)f_sig->Get("h_lep2dxy");
	TH1F* h_lep1dz_sig = (TH1F*)f_sig->Get("h_lep1dz");
	TH1F* h_lep2dz_sig = (TH1F*)f_sig->Get("h_lep2dz");
	//TH1F* h_tau1dz_sig = (TH1F*)f_sig->Get("h_tau1dz");
	TH1F* h_lepmaxdxy_sig = (TH1F*)f_sig->Get("h_lepmaxdxy");
	TH1F* h_visM_taulep1_sig = (TH1F*)f_sig->Get("h_visM_taulep1");
	TH1F* h_visM_taulep2_sig = (TH1F*)f_sig->Get("h_visM_taulep2");
	TH1F* h_visM_taulepmaxdxy_sig = (TH1F*)f_sig->Get("h_visM_taulepmaxdxy");
	TH1F* h_visM_taulepmindR_sig = (TH1F*)f_sig->Get("h_visM_taulepmindR");
	
	TH1F* h_lep1pt_TTJets = (TH1F*)f_TTJets->Get("h_lep1pt");
	TH1F* h_lep1eta_TTJets = (TH1F*)f_TTJets->Get("h_lep1eta");
	TH1F* h_lep2pt_TTJets = (TH1F*)f_TTJets->Get("h_lep2pt");
	TH1F* h_lep2eta_TTJets = (TH1F*)f_TTJets->Get("h_lep2eta");
	TH1F* h_taupt_TTJets = (TH1F*)f_TTJets->Get("h_taupt");
	TH1F* h_taueta_TTJets = (TH1F*)f_TTJets->Get("h_taueta");
	TH1F* h_bjet1pt_TTJets = (TH1F*)f_TTJets->Get("h_bjet1pt");
	TH1F* h_bjet1eta_TTJets = (TH1F*)f_TTJets->Get("h_bjet1eta");
	TH1F* h_dRlep1tau_TTJets = (TH1F*)f_TTJets->Get("h_dRlep1tau");
	TH1F* h_dRlep2tau_TTJets = (TH1F*)f_TTJets->Get("h_dRlep2tau");
	TH1F* h_lep1dxy_TTJets = (TH1F*)f_TTJets->Get("h_lep1dxy");
	TH1F* h_lep2dxy_TTJets = (TH1F*)f_TTJets->Get("h_lep2dxy");
	TH1F* h_lep1dz_TTJets = (TH1F*)f_TTJets->Get("h_lep1dz");
	TH1F* h_lep2dz_TTJets = (TH1F*)f_TTJets->Get("h_lep2dz");
	//TH1F* h_tau1dz_TTJets = (TH1F*)f_TTJets->Get("h_tau1dz");
	TH1F* h_lepmaxdxy_TTJets = (TH1F*)f_TTJets->Get("h_lepmaxdxy");
	TH1F* h_visM_taulep1_TTJets = (TH1F*)f_TTJets->Get("h_visM_taulep1");
	TH1F* h_visM_taulep2_TTJets = (TH1F*)f_TTJets->Get("h_visM_taulep2");
	TH1F* h_visM_taulepmaxdxy_TTJets = (TH1F*)f_TTJets->Get("h_visM_taulepmaxdxy");
	TH1F* h_visM_taulepmindR_TTJets = (TH1F*)f_TTJets->Get("h_visM_taulepmindR");
	
	/*
	h_lep1pt_sig -> Sumw2();
	h_lep1eta_sig -> Sumw2();
	h_lep2pt_sig -> Sumw2();
	h_lep2eta_sig -> Sumw2();
	h_taupt_sig -> Sumw2();
	h_taueta_sig -> Sumw2();
	h_bjet1pt_sig -> Sumw2();
	h_bjet1eta_sig -> Sumw2();
	h_dRlep1tau_sig -> Sumw2();
	h_dRlep2tau_sig -> Sumw2();

	h_lep1pt_TTJets -> Sumw2();
	h_lep1eta_TTJets -> Sumw2();
	h_lep2pt_TTJets -> Sumw2();
	h_lep2eta_TTJets -> Sumw2();
	h_taupt_TTJets -> Sumw2();
	h_taueta_TTJets -> Sumw2();
	h_bjet1pt_TTJets -> Sumw2();
	h_bjet1eta_TTJets -> Sumw2();
	h_dRlep1tau_TTJets -> Sumw2();
	h_dRlep2tau_TTJets -> Sumw2();
	*/
	
	h_lep1pt_sig -> Scale(1.0/nsig);
	h_lep1eta_sig -> Scale(1.0/nsig);
	h_lep2pt_sig -> Scale(1.0/nsig);
	h_lep2eta_sig -> Scale(1.0/nsig);
	h_taupt_sig -> Scale(1.0/nsig);
	h_taueta_sig -> Scale(1.0/nsig);
	h_bjet1pt_sig -> Scale(1.0/nsig);
	h_bjet1eta_sig -> Scale(1.0/nsig);
	h_dRlep1tau_sig -> Scale(1.0/nsig);
	h_dRlep2tau_sig -> Scale(1.0/nsig);
	h_lep1dxy_sig -> Scale(1.0/nsig);
	h_lep1dz_sig -> Scale(1.0/nsig);
	h_lep2dxy_sig -> Scale(1.0/nsig);
	h_lep2dz_sig -> Scale(1.0/nsig);
	//h_tau1dz_sig -> Scale(1.0/nsig);
	h_lepmaxdxy_sig -> Scale(1.0/nsig);
	h_visM_taulep1_sig -> Scale(1.0/nsig);
	h_visM_taulep2_sig -> Scale(1.0/nsig);
	h_visM_taulepmaxdxy_sig -> Scale(1.0/nsig);
	h_visM_taulepmindR_sig -> Scale(1.0/nsig);
	
	h_lep1pt_TTJets -> Scale(1.0/nTTJets);
	h_lep1eta_TTJets -> Scale(1.0/nTTJets);
	h_lep2pt_TTJets -> Scale(1.0/nTTJets);
	h_lep2eta_TTJets -> Scale(1.0/nTTJets);
	h_taupt_TTJets -> Scale(1.0/nTTJets);
	h_taueta_TTJets -> Scale(1.0/nTTJets);
	h_bjet1pt_TTJets -> Scale(1.0/nTTJets);
	h_bjet1eta_TTJets -> Scale(1.0/nTTJets);
	h_dRlep1tau_TTJets -> Scale(1.0/nTTJets);
	h_dRlep2tau_TTJets -> Scale(1.0/nTTJets);
	h_lep1dxy_TTJets -> Scale(1.0/nTTJets);
	h_lep1dz_TTJets -> Scale(1.0/nTTJets);
	h_lep2dxy_TTJets -> Scale(1.0/nTTJets);
	h_lep2dz_TTJets -> Scale(1.0/nTTJets);
	//h_tau1dz_TTJets -> Scale(1.0/nTTJets);
	h_lepmaxdxy_TTJets -> Scale(1.0/nTTJets);
	h_visM_taulep1_TTJets -> Scale(1.0/nTTJets);
	h_visM_taulep2_TTJets -> Scale(1.0/nTTJets);
	h_visM_taulepmaxdxy_TTJets -> Scale(1.0/nTTJets);
	h_visM_taulepmindR_TTJets -> Scale(1.0/nTTJets);
	
	TCanvas c;

	//// h_lep1pt
	h_lep1pt_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lep1pt_TTJets->GetXaxis()->SetTitle("Leading Lepton p_{T} [GeV]");
	h_lep1pt_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_lep1pt_TTJets->SetLineColor(kBlue);
	float lep1pt_max =
		max(h_lep1pt_TTJets->GetMaximum(), h_lep1pt_sig->GetMaximum());
	h_lep1pt_TTJets->SetMaximum(lep1pt_max*1.1);
	h_lep1pt_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_lep1pt_bkg = (TPaveStats*) h_lep1pt_TTJets->GetListOfFunctions()->FindObject("stats");
	st_lep1pt_bkg->SetTextColor(kBlue);
	st_lep1pt_bkg->SetOptStat(1110);
	h_lep1pt_sig->SetLineColor(kRed);
	h_lep1pt_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_lep1pt_sig = (TPaveStats*) h_lep1pt_sig->GetListOfFunctions()->FindObject("stats");
	st_lep1pt_sig->SetTextColor(kRed);
	st_lep1pt_sig->SetOptStat(1110);
	st_lep1pt_bkg->SetY1NDC(2*(st_lep1pt_sig->GetY1NDC())-st_lep1pt_sig->GetY2NDC());
	st_lep1pt_bkg->SetY2NDC(st_lep1pt_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_lep1pt = new TLegend(0.72,0.25,0.85,0.38);
	leg_lep1pt -> AddEntry(h_lep1pt_sig, "signal", "L");
	leg_lep1pt -> AddEntry(h_lep1pt_TTJets, "tt+Jets", "L");
	leg_lep1pt -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lep1pt.pdf");

	//// h_lep1eta
	h_lep1eta_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lep1eta_TTJets->GetXaxis()->SetTitle("Leading Lepton #eta");
	h_lep1eta_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_lep1eta_TTJets->SetLineColor(kBlue);
	float lep1eta_max =
		max(h_lep1eta_TTJets->GetMaximum(), h_lep1eta_sig->GetMaximum());
	h_lep1eta_TTJets->SetMaximum(lep1eta_max*1.1);
	h_lep1eta_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_lep1eta_bkg = (TPaveStats*) h_lep1eta_TTJets->GetListOfFunctions()->FindObject("stats");
	st_lep1eta_bkg->SetTextColor(kBlue);
	st_lep1eta_bkg->SetOptStat(1110);
	h_lep1eta_sig->SetLineColor(kRed);
	h_lep1eta_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_lep1eta_sig = (TPaveStats*) h_lep1eta_sig->GetListOfFunctions()->FindObject("stats");
	st_lep1eta_sig->SetTextColor(kRed);
	st_lep1eta_sig->SetOptStat(1110);
	st_lep1eta_bkg->SetY1NDC(2*(st_lep1eta_sig->GetY1NDC())-st_lep1eta_sig->GetY2NDC());
	st_lep1eta_bkg->SetY2NDC(st_lep1eta_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_lep1eta = new TLegend(0.72,0.25,0.85,0.38);
	leg_lep1eta -> AddEntry(h_lep1eta_sig, "signal", "L");
	leg_lep1eta -> AddEntry(h_lep1eta_TTJets, "tt+Jets", "L");
	leg_lep1eta -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lep1eta.pdf");

	//// h_lep2pt
	h_lep2pt_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lep2pt_TTJets->GetXaxis()->SetTitle("Subleading Lepton p_{T} [GeV]");
	h_lep2pt_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_lep2pt_TTJets->SetLineColor(kBlue);
	float lep2pt_max =
		max(h_lep2pt_TTJets->GetMaximum(), h_lep2pt_sig->GetMaximum());
	h_lep2pt_TTJets->SetMaximum(lep2pt_max*1.1);
	h_lep2pt_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_lep2pt_bkg = (TPaveStats*) h_lep2pt_TTJets->GetListOfFunctions()->FindObject("stats");
	st_lep2pt_bkg->SetTextColor(kBlue);
	st_lep2pt_bkg->SetOptStat(1110);
	h_lep2pt_sig->SetLineColor(kRed);
	h_lep2pt_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_lep2pt_sig = (TPaveStats*) h_lep2pt_sig->GetListOfFunctions()->FindObject("stats");
	st_lep2pt_sig->SetTextColor(kRed);
	st_lep2pt_sig->SetOptStat(1110);
	st_lep2pt_bkg->SetY1NDC(2*(st_lep2pt_sig->GetY1NDC())-st_lep2pt_sig->GetY2NDC());
	st_lep2pt_bkg->SetY2NDC(st_lep2pt_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_lep2pt = new TLegend(0.72,0.25,0.85,0.38);
	leg_lep2pt -> AddEntry(h_lep2pt_sig, "signal", "L");
	leg_lep2pt -> AddEntry(h_lep2pt_TTJets, "tt+Jets", "L");
	leg_lep2pt -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lep2pt.pdf");

	//// h_lep2eta
	h_lep2eta_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lep2eta_TTJets->GetXaxis()->SetTitle("Subleading Lepton #eta");
	h_lep2eta_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_lep2eta_TTJets->SetLineColor(kBlue);
	float lep2eta_max =
		max(h_lep2eta_TTJets->GetMaximum(), h_lep2eta_sig->GetMaximum());
	h_lep2eta_TTJets->SetMaximum(lep2eta_max*1.1);
	h_lep2eta_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_lep2eta_bkg = (TPaveStats*) h_lep2eta_TTJets->GetListOfFunctions()->FindObject("stats");
	st_lep2eta_bkg->SetTextColor(kBlue);
	st_lep2eta_bkg->SetOptStat(1110);
	h_lep2eta_sig->SetLineColor(kRed);
	h_lep2eta_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_lep2eta_sig = (TPaveStats*) h_lep2eta_sig->GetListOfFunctions()->FindObject("stats");
	st_lep2eta_sig->SetTextColor(kRed);
	st_lep2eta_sig->SetOptStat(1110);
	st_lep2eta_bkg->SetY1NDC(2*(st_lep2eta_sig->GetY1NDC())-st_lep2eta_sig->GetY2NDC());
	st_lep2eta_bkg->SetY2NDC(st_lep2eta_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_lep2eta = new TLegend(0.72,0.25,0.85,0.38);
	leg_lep2eta -> AddEntry(h_lep2eta_sig, "signal", "L");
	leg_lep2eta -> AddEntry(h_lep2eta_TTJets, "tt+Jets", "L");
	leg_lep2eta -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lep2eta.pdf");

	//// h_taupt
	h_taupt_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_taupt_TTJets->GetXaxis()->SetTitle("#tau p_{T} [GeV]");
	h_taupt_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_taupt_TTJets->SetLineColor(kBlue);
	float taupt_max =
		max(h_taupt_TTJets->GetMaximum(), h_taupt_sig->GetMaximum());
	h_taupt_TTJets->SetMaximum(taupt_max*1.1);
	h_taupt_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_taupt_bkg = (TPaveStats*) h_taupt_TTJets->GetListOfFunctions()->FindObject("stats");
	st_taupt_bkg->SetTextColor(kBlue);
	st_taupt_bkg->SetOptStat(1110);
	h_taupt_sig->SetLineColor(kRed);
	h_taupt_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_taupt_sig = (TPaveStats*) h_taupt_sig->GetListOfFunctions()->FindObject("stats");
	st_taupt_sig->SetTextColor(kRed);
	st_taupt_sig->SetOptStat(1110);
	st_taupt_bkg->SetY1NDC(2*(st_taupt_sig->GetY1NDC())-st_taupt_sig->GetY2NDC());
	st_taupt_bkg->SetY2NDC(st_taupt_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_taupt = new TLegend(0.72,0.25,0.85,0.38);
	leg_taupt -> AddEntry(h_taupt_sig, "signal", "L");
	leg_taupt -> AddEntry(h_taupt_TTJets, "tt+Jets", "L");
	leg_taupt -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/taupt.pdf");

	//// h_taueta
	h_taueta_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_taueta_TTJets->GetXaxis()->SetTitle("#tau #eta");
	h_taueta_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_taueta_TTJets->SetLineColor(kBlue);
	float taueta_max =
		max(h_taueta_TTJets->GetMaximum(), h_taueta_sig->GetMaximum());
	h_taueta_TTJets->SetMaximum(taueta_max*1.1);
	h_taueta_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_taueta_bkg = (TPaveStats*) h_taueta_TTJets->GetListOfFunctions()->FindObject("stats");
	st_taueta_bkg->SetTextColor(kBlue);
	st_taueta_bkg->SetOptStat(1110);
	h_taueta_sig->SetLineColor(kRed);
	h_taueta_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_taueta_sig = (TPaveStats*) h_taueta_sig->GetListOfFunctions()->FindObject("stats");
	st_taueta_sig->SetTextColor(kRed);
	st_taueta_sig->SetOptStat(1110);
	st_taueta_bkg->SetY1NDC(2*(st_taueta_sig->GetY1NDC())-st_taueta_sig->GetY2NDC());
	st_taueta_bkg->SetY2NDC(st_taueta_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_taueta = new TLegend(0.72,0.25,0.85,0.38);
	leg_taueta -> AddEntry(h_taueta_sig, "signal", "L");
	leg_taueta -> AddEntry(h_taueta_TTJets, "tt+Jets", "L");
	leg_taueta -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/taueta.pdf");

	//// h_bjet1pt
	h_bjet1pt_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_bjet1pt_TTJets->GetXaxis()->SetTitle("Leading b-jet p_{T} [GeV]");
	h_bjet1pt_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_bjet1pt_TTJets->SetLineColor(kBlue);
	float bjet1pt_max =
		max(h_bjet1pt_TTJets->GetMaximum(), h_bjet1pt_sig->GetMaximum());
	h_bjet1pt_TTJets->SetMaximum(bjet1pt_max*1.1);
	h_bjet1pt_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_bjet1pt_bkg = (TPaveStats*) h_bjet1pt_TTJets->GetListOfFunctions()->FindObject("stats");
	st_bjet1pt_bkg->SetTextColor(kBlue);
	st_bjet1pt_bkg->SetOptStat(1110);
	h_bjet1pt_sig->SetLineColor(kRed);
	h_bjet1pt_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_bjet1pt_sig = (TPaveStats*) h_bjet1pt_sig->GetListOfFunctions()->FindObject("stats");
	st_bjet1pt_sig->SetTextColor(kRed);
	st_bjet1pt_sig->SetOptStat(1110);
	st_bjet1pt_bkg->SetY1NDC(2*(st_bjet1pt_sig->GetY1NDC())-st_bjet1pt_sig->GetY2NDC());
	st_bjet1pt_bkg->SetY2NDC(st_bjet1pt_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_bjet1pt = new TLegend(0.72,0.25,0.85,0.38);
	leg_bjet1pt -> AddEntry(h_bjet1pt_sig, "signal", "L");
	leg_bjet1pt -> AddEntry(h_bjet1pt_TTJets, "tt+Jets", "L");
	leg_bjet1pt -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/bjet1pt.pdf");

	//// h_bjet1eta
	h_bjet1eta_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_bjet1eta_TTJets->GetXaxis()->SetTitle("Leading b-jet #eta");
	h_bjet1eta_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_bjet1eta_TTJets->SetLineColor(kBlue);
	float bjet1eta_max =
		max(h_bjet1eta_TTJets->GetMaximum(), h_bjet1eta_sig->GetMaximum());
	h_bjet1eta_TTJets->SetMaximum(bjet1eta_max*1.1);
	h_bjet1eta_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_bjet1eta_bkg = (TPaveStats*) h_bjet1eta_TTJets->GetListOfFunctions()->FindObject("stats");
	st_bjet1eta_bkg->SetTextColor(kBlue);
	st_bjet1eta_bkg->SetOptStat(1110);
	h_bjet1eta_sig->SetLineColor(kRed);
	h_bjet1eta_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_bjet1eta_sig = (TPaveStats*) h_bjet1eta_sig->GetListOfFunctions()->FindObject("stats");
	st_bjet1eta_sig->SetTextColor(kRed);
	st_bjet1eta_sig->SetOptStat(1110);
	st_bjet1eta_bkg->SetY1NDC(2*(st_bjet1eta_sig->GetY1NDC())-st_bjet1eta_sig->GetY2NDC());
	st_bjet1eta_bkg->SetY2NDC(st_bjet1eta_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_bjet1eta = new TLegend(0.72,0.25,0.85,0.38);
	leg_bjet1eta -> AddEntry(h_bjet1eta_sig, "signal", "L");
	leg_bjet1eta -> AddEntry(h_bjet1eta_TTJets, "tt+Jets", "L");
	leg_bjet1eta -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/bjet1eta.pdf");

	//// h_dRlep1tau
	h_dRlep1tau_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_dRlep1tau_TTJets->GetXaxis()->SetTitle("#DeltaR(l_{ldg},#tau_{h})");
	h_dRlep1tau_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_dRlep1tau_TTJets->SetLineColor(kBlue);
	float dRlep1tau_max =
		max(h_dRlep1tau_TTJets->GetMaximum(), h_dRlep1tau_sig->GetMaximum());
	h_dRlep1tau_TTJets->SetMaximum(dRlep1tau_max*1.1);
	h_dRlep1tau_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_dRlep1tau_bkg = (TPaveStats*) h_dRlep1tau_TTJets->GetListOfFunctions()->FindObject("stats");
	st_dRlep1tau_bkg->SetTextColor(kBlue);
	st_dRlep1tau_bkg->SetOptStat(1110);
	h_dRlep1tau_sig->SetLineColor(kRed);
	h_dRlep1tau_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_dRlep1tau_sig = (TPaveStats*) h_dRlep1tau_sig->GetListOfFunctions()->FindObject("stats");
	st_dRlep1tau_sig->SetTextColor(kRed);
	st_dRlep1tau_sig->SetOptStat(1110);
	st_dRlep1tau_bkg->SetY1NDC(2*(st_dRlep1tau_sig->GetY1NDC())-st_dRlep1tau_sig->GetY2NDC());
	st_dRlep1tau_bkg->SetY2NDC(st_dRlep1tau_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_dRlep1tau = new TLegend(0.72,0.25,0.85,0.38);
	leg_dRlep1tau -> AddEntry(h_dRlep1tau_sig, "signal", "L");
	leg_dRlep1tau -> AddEntry(h_dRlep1tau_TTJets, "tt+Jets", "L");
	leg_dRlep1tau -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dRlep1tau.pdf");

	//// h_dRlep2tau
	h_dRlep2tau_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_dRlep2tau_TTJets->GetXaxis()->SetTitle("#DeltaR(l_{sub-ldg},#tau_{h})");
	h_dRlep2tau_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_dRlep2tau_TTJets->SetLineColor(kBlue);
	float dRlep2tau_max =
		max(h_dRlep2tau_TTJets->GetMaximum(), h_dRlep2tau_sig->GetMaximum());
	h_dRlep2tau_TTJets->SetMaximum(dRlep2tau_max*1.1);
	h_dRlep2tau_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_dRlep2tau_bkg = (TPaveStats*) h_dRlep2tau_TTJets->GetListOfFunctions()->FindObject("stats");
	st_dRlep2tau_bkg->SetTextColor(kBlue);
	st_dRlep2tau_bkg->SetOptStat(1110);
	h_dRlep2tau_sig->SetLineColor(kRed);
	h_dRlep2tau_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_dRlep2tau_sig = (TPaveStats*) h_dRlep2tau_sig->GetListOfFunctions()->FindObject("stats");
	st_dRlep2tau_sig->SetTextColor(kRed);
	st_dRlep2tau_sig->SetOptStat(1110);
	st_dRlep2tau_bkg->SetY1NDC(2*(st_dRlep2tau_sig->GetY1NDC())-st_dRlep2tau_sig->GetY2NDC());
	st_dRlep2tau_bkg->SetY2NDC(st_dRlep2tau_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_dRlep2tau = new TLegend(0.72,0.25,0.85,0.38);
	leg_dRlep2tau -> AddEntry(h_dRlep2tau_sig, "signal", "L");
	leg_dRlep2tau -> AddEntry(h_dRlep2tau_TTJets, "tt+Jets", "L");
	leg_dRlep2tau -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dRlep2tau.pdf");

	//// h_lep1dxy
	h_lep1dxy_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lep1dxy_TTJets->GetXaxis()->SetTitle("dxy_{ldg lep}");
	h_lep1dxy_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_lep1dxy_TTJets->SetLineColor(kBlue);
	float lep1dxy_max =
		max(h_lep1dxy_TTJets->GetMaximum(), h_lep1dxy_sig->GetMaximum());
	h_lep1dxy_TTJets->SetMaximum(lep1dxy_max*1.1);
	h_lep1dxy_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_lep1dxy_bkg = (TPaveStats*) h_lep1dxy_TTJets->GetListOfFunctions()->FindObject("stats");
	st_lep1dxy_bkg->SetTextColor(kBlue);
	st_lep1dxy_bkg->SetOptStat(1110);
	h_lep1dxy_sig->SetLineColor(kRed);
	h_lep1dxy_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_lep1dxy_sig = (TPaveStats*) h_lep1dxy_sig->GetListOfFunctions()->FindObject("stats");
	st_lep1dxy_sig->SetTextColor(kRed);
	st_lep1dxy_sig->SetOptStat(1110);
	st_lep1dxy_bkg->SetY1NDC(2*(st_lep1dxy_sig->GetY1NDC())-st_lep1dxy_sig->GetY2NDC());
	st_lep1dxy_bkg->SetY2NDC(st_lep1dxy_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_lep1dxy = new TLegend(0.72,0.25,0.85,0.38);
	leg_lep1dxy -> AddEntry(h_lep1dxy_sig, "signal", "L");
	leg_lep1dxy -> AddEntry(h_lep1dxy_TTJets, "tt+Jets", "L");
	leg_lep1dxy -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lep1dxy.pdf");

	//// h_lep2dxy
	h_lep2dxy_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lep2dxy_TTJets->GetXaxis()->SetTitle("dxy_{sub-ldg lep}");
	h_lep2dxy_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_lep2dxy_TTJets->SetLineColor(kBlue);
	float lep2dxy_max =
		max(h_lep2dxy_TTJets->GetMaximum(), h_lep2dxy_sig->GetMaximum());
	h_lep2dxy_TTJets->SetMaximum(lep2dxy_max*1.1);
	h_lep2dxy_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_lep2dxy_bkg = (TPaveStats*) h_lep2dxy_TTJets->GetListOfFunctions()->FindObject("stats");
	st_lep2dxy_bkg->SetTextColor(kBlue);
	st_lep2dxy_bkg->SetOptStat(1110);
	h_lep2dxy_sig->SetLineColor(kRed);
	h_lep2dxy_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_lep2dxy_sig = (TPaveStats*) h_lep2dxy_sig->GetListOfFunctions()->FindObject("stats");
	st_lep2dxy_sig->SetTextColor(kRed);
	st_lep2dxy_sig->SetOptStat(1110);
	st_lep2dxy_bkg->SetY1NDC(2*(st_lep2dxy_sig->GetY1NDC())-st_lep2dxy_sig->GetY2NDC());
	st_lep2dxy_bkg->SetY2NDC(st_lep2dxy_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_lep2dxy = new TLegend(0.72,0.25,0.85,0.38);
	leg_lep2dxy -> AddEntry(h_lep2dxy_sig, "signal", "L");
	leg_lep2dxy -> AddEntry(h_lep2dxy_TTJets, "tt+Jets", "L");
	leg_lep2dxy -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lep2dxy.pdf");

	//// h_lep1dz
	h_lep1dz_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lep1dz_TTJets->GetXaxis()->SetTitle("dz_{ldg lep}");
	h_lep1dz_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_lep1dz_TTJets->SetLineColor(kBlue);
	float lep1dz_max =
		max(h_lep1dz_TTJets->GetMaximum(), h_lep1dz_sig->GetMaximum());
	h_lep1dz_TTJets->SetMaximum(lep1dz_max*1.1);
	h_lep1dz_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_lep1dz_bkg = (TPaveStats*) h_lep1dz_TTJets->GetListOfFunctions()->FindObject("stats");
	st_lep1dz_bkg->SetTextColor(kBlue);
	st_lep1dz_bkg->SetOptStat(1110);
	h_lep1dz_sig->SetLineColor(kRed);
	h_lep1dz_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_lep1dz_sig = (TPaveStats*) h_lep1dz_sig->GetListOfFunctions()->FindObject("stats");
	st_lep1dz_sig->SetTextColor(kRed);
	st_lep1dz_sig->SetOptStat(1110);
	st_lep1dz_bkg->SetY1NDC(2*(st_lep1dz_sig->GetY1NDC())-st_lep1dz_sig->GetY2NDC());
	st_lep1dz_bkg->SetY2NDC(st_lep1dz_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_lep1dz = new TLegend(0.72,0.25,0.85,0.38);
	leg_lep1dz -> AddEntry(h_lep1dz_sig, "signal", "L");
	leg_lep1dz -> AddEntry(h_lep1dz_TTJets, "tt+Jets", "L");
	leg_lep1dz -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lep1dz.pdf");

	//// h_lep2dz
	h_lep2dz_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lep2dz_TTJets->GetXaxis()->SetTitle("dz_{sub-ldg lep}");
	h_lep2dz_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_lep2dz_TTJets->SetLineColor(kBlue);
	float lep2dz_max =
		max(h_lep2dz_TTJets->GetMaximum(), h_lep2dz_sig->GetMaximum());
	h_lep2dz_TTJets->SetMaximum(lep2dz_max*1.1);
	h_lep2dz_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_lep2dz_bkg = (TPaveStats*) h_lep2dz_TTJets->GetListOfFunctions()->FindObject("stats");
	st_lep2dz_bkg->SetTextColor(kBlue);
	st_lep2dz_bkg->SetOptStat(1110);
	h_lep2dz_sig->SetLineColor(kRed);
	h_lep2dz_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_lep2dz_sig = (TPaveStats*) h_lep2dz_sig->GetListOfFunctions()->FindObject("stats");
	st_lep2dz_sig->SetTextColor(kRed);
	st_lep2dz_sig->SetOptStat(1110);
	st_lep2dz_bkg->SetY1NDC(2*(st_lep2dz_sig->GetY1NDC())-st_lep2dz_sig->GetY2NDC());
	st_lep2dz_bkg->SetY2NDC(st_lep2dz_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_lep2dz = new TLegend(0.72,0.25,0.85,0.38);
	leg_lep2dz -> AddEntry(h_lep2dz_sig, "signal", "L");
	leg_lep2dz -> AddEntry(h_lep2dz_TTJets, "tt+Jets", "L");
	leg_lep2dz -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lep2dz.pdf");

	/*
	h_tau1dz_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_tau1dz_TTJets->GetXaxis()->SetTitle("dz_{#tau_{h}}");
	h_tau1dz_TTJets->GetYaxis()->SetTitle("Normalized");
	h_tau1dz_TTJets->SetLineColor(kBlue);
	h_tau1dz_TTJets->Draw("");
	h_tau1dz_sig->SetLineColor(kRed);
	h_tau1dz_sig->Draw("same");
	TLegend *leg_tau1dz = new TLegend(0.72,0.25,0.85,0.38);
	leg_tau1dz -> AddEntry(h_tau1dz_sig, "signal", "L");
	leg_tau1dz -> AddEntry(h_tau1dz_TTJets, "tt+Jets", "L");
	leg_tau1dz -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/tau1dz.pdf");
	*/

	//// h_lepmaxdxy
	h_lepmaxdxy_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_lepmaxdxy_TTJets->GetXaxis()->SetTitle("max(dxy_{ldg lep},dxy_{sub-ldg lep})");
	h_lepmaxdxy_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_lepmaxdxy_TTJets->SetLineColor(kBlue);
	float lepmaxdxy_max =
		max(h_lepmaxdxy_TTJets->GetMaximum(), h_lepmaxdxy_sig->GetMaximum());
	h_lepmaxdxy_TTJets->SetMaximum(lepmaxdxy_max*1.1);
	h_lepmaxdxy_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_lepmaxdxy_bkg = (TPaveStats*) h_lepmaxdxy_TTJets->GetListOfFunctions()->FindObject("stats");
	st_lepmaxdxy_bkg->SetTextColor(kBlue);
	st_lepmaxdxy_bkg->SetOptStat(1110);
	h_lepmaxdxy_sig->SetLineColor(kRed);
	h_lepmaxdxy_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_lepmaxdxy_sig = (TPaveStats*) h_lepmaxdxy_sig->GetListOfFunctions()->FindObject("stats");
	st_lepmaxdxy_sig->SetTextColor(kRed);
	st_lepmaxdxy_sig->SetOptStat(1110);
	st_lepmaxdxy_bkg->SetY1NDC(2*(st_lepmaxdxy_sig->GetY1NDC())-st_lepmaxdxy_sig->GetY2NDC());
	st_lepmaxdxy_bkg->SetY2NDC(st_lepmaxdxy_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_lepmaxdxy = new TLegend(0.72,0.25,0.85,0.38);
	leg_lepmaxdxy -> AddEntry(h_lepmaxdxy_sig, "signal", "L");
	leg_lepmaxdxy -> AddEntry(h_lepmaxdxy_TTJets, "tt+Jets", "L");
	leg_lepmaxdxy -> Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/lepmaxdxy.pdf");

	//// h_visM_taulep1
	h_visM_taulep1_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_visM_taulep1_TTJets->GetXaxis()->SetTitle("visMass(#tau_{h},l_{ldg})");
	h_visM_taulep1_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_visM_taulep1_TTJets->SetLineColor(kBlue);
	float visM_taulep1_max =
		max(h_visM_taulep1_TTJets->GetMaximum(),h_visM_taulep1_sig->GetMaximum());
	h_visM_taulep1_TTJets->SetMaximum(visM_taulep1_max*1.1);
	h_visM_taulep1_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_visM_taulep1_bkg = (TPaveStats*) h_visM_taulep1_TTJets->GetListOfFunctions()->FindObject("stats");
	st_visM_taulep1_bkg->SetTextColor(kBlue);
	st_visM_taulep1_bkg->SetOptStat(1110);
	h_visM_taulep1_sig->SetLineColor(kRed);
	h_visM_taulep1_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_visM_taulep1_sig = (TPaveStats*) h_visM_taulep1_sig->GetListOfFunctions()->FindObject("stats");
	st_visM_taulep1_sig->SetTextColor(kRed);
	st_visM_taulep1_sig->SetOptStat(1110);
	st_visM_taulep1_bkg->SetY1NDC(2*(st_visM_taulep1_sig->GetY1NDC())-st_visM_taulep1_sig->GetY2NDC());
	st_visM_taulep1_bkg->SetY2NDC(st_visM_taulep1_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_visM_taulep1 = new TLegend(0.72,0.25,0.85,0.38);
	leg_visM_taulep1 -> AddEntry(h_visM_taulep1_sig, "signal", "L");
	leg_visM_taulep1 -> AddEntry(h_visM_taulep1_TTJets, "tt+Jets", "L");
	leg_visM_taulep1->Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/visM_taulep1.pdf");

	//// h_visM_taulep2
	h_visM_taulep2_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_visM_taulep2_TTJets->GetXaxis()->SetTitle("visMass(#tau_{h},l_{subldg})");
	h_visM_taulep2_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_visM_taulep2_TTJets->SetLineColor(kBlue);
	float visM_taulep2_max =
		max(h_visM_taulep2_TTJets->GetMaximum(),h_visM_taulep2_sig->GetMaximum());
	h_visM_taulep2_TTJets->SetMaximum(visM_taulep2_max*1.1);
	h_visM_taulep2_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_visM_taulep2_bkg = (TPaveStats*) h_visM_taulep2_TTJets->GetListOfFunctions()->FindObject("stats");
	st_visM_taulep2_bkg->SetTextColor(kBlue);
	st_visM_taulep2_bkg->SetOptStat(1110);
	h_visM_taulep2_sig->SetLineColor(kRed);
	h_visM_taulep2_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_visM_taulep2_sig = (TPaveStats*) h_visM_taulep2_sig->GetListOfFunctions()->FindObject("stats");
	st_visM_taulep2_sig->SetTextColor(kRed);
	st_visM_taulep2_sig->SetOptStat(1110);
	st_visM_taulep2_bkg->SetY1NDC(2*(st_visM_taulep2_sig->GetY1NDC())-st_visM_taulep2_sig->GetY2NDC());
	st_visM_taulep2_bkg->SetY2NDC(st_visM_taulep2_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_visM_taulep2 = new TLegend(0.72,0.25,0.85,0.38);
	leg_visM_taulep2 -> AddEntry(h_visM_taulep2_sig, "signal", "L");
	leg_visM_taulep2 -> AddEntry(h_visM_taulep2_TTJets, "tt+Jets", "L");
	leg_visM_taulep2->Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/visM_taulep2.pdf");

	//// h_visM_taulepmaxdxy
	h_visM_taulepmaxdxy_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_visM_taulepmaxdxy_TTJets->GetXaxis()->SetTitle("visMass(#tau_{h},l_{max(dxy)})");
	h_visM_taulepmaxdxy_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_visM_taulepmaxdxy_TTJets->SetLineColor(kBlue);
	float visM_taulepmaxdxy_max =
		max(h_visM_taulepmaxdxy_TTJets->GetMaximum(),h_visM_taulepmaxdxy_sig->GetMaximum());
	h_visM_taulepmaxdxy_TTJets->SetMaximum(visM_taulepmaxdxy_max*1.1);
	h_visM_taulepmaxdxy_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_visM_taulepmaxdxy_bkg = (TPaveStats*) h_visM_taulepmaxdxy_TTJets->GetListOfFunctions()->FindObject("stats");
	st_visM_taulepmaxdxy_bkg->SetTextColor(kBlue);
	st_visM_taulepmaxdxy_bkg->SetOptStat(1110);
	h_visM_taulepmaxdxy_sig->SetLineColor(kRed);
	h_visM_taulepmaxdxy_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_visM_taulepmaxdxy_sig = (TPaveStats*) h_visM_taulepmaxdxy_sig->GetListOfFunctions()->FindObject("stats");
	st_visM_taulepmaxdxy_sig->SetTextColor(kRed);
	st_visM_taulepmaxdxy_sig->SetOptStat(1110);
	st_visM_taulepmaxdxy_bkg->SetY1NDC(2*(st_visM_taulepmaxdxy_sig->GetY1NDC())-st_visM_taulepmaxdxy_sig->GetY2NDC());
	st_visM_taulepmaxdxy_bkg->SetY2NDC(st_visM_taulepmaxdxy_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_visM_taulepmaxdxy = new TLegend(0.72,0.25,0.85,0.38);
	leg_visM_taulepmaxdxy -> AddEntry(h_visM_taulepmaxdxy_sig, "signal", "L");
	leg_visM_taulepmaxdxy -> AddEntry(h_visM_taulepmaxdxy_TTJets, "tt+Jets", "L");
	leg_visM_taulepmaxdxy->Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/visM_taulepmaxdxy.pdf");

	//// h_visM_taulepmindR
	h_visM_taulepmindR_TTJets->SetTitle("same sign leptons + #tau_{h}");
	h_visM_taulepmindR_TTJets->GetXaxis()->SetTitle("visMass(#tau_{h},l_{min(dR)})");
	h_visM_taulepmindR_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_visM_taulepmindR_TTJets->SetLineColor(kBlue);
	float visM_taulepmindR_max =
		max(h_visM_taulepmindR_TTJets->GetMaximum(),h_visM_taulepmindR_sig->GetMaximum());
	h_visM_taulepmindR_TTJets->SetMaximum(visM_taulepmindR_max*1.1);
	h_visM_taulepmindR_TTJets->Draw("");
	gPad->Update();
	TPaveStats *st_visM_taulepmindR_bkg = (TPaveStats*) h_visM_taulepmindR_TTJets->GetListOfFunctions()->FindObject("stats");
	st_visM_taulepmindR_bkg->SetTextColor(kBlue);
	st_visM_taulepmindR_bkg->SetOptStat(1110);
	h_visM_taulepmindR_sig->SetLineColor(kRed);
	h_visM_taulepmindR_sig->Draw("sames");
	gPad->Update();
	TPaveStats *st_visM_taulepmindR_sig = (TPaveStats*) h_visM_taulepmindR_sig->GetListOfFunctions()->FindObject("stats");
	st_visM_taulepmindR_sig->SetTextColor(kRed);
	st_visM_taulepmindR_sig->SetOptStat(1110);
	st_visM_taulepmindR_bkg->SetY1NDC(2*(st_visM_taulepmindR_sig->GetY1NDC())-st_visM_taulepmindR_sig->GetY2NDC());
	st_visM_taulepmindR_bkg->SetY2NDC(st_visM_taulepmindR_sig->GetY1NDC());
	gPad->Modified();
	
	TLegend *leg_visM_taulepmindR = new TLegend(0.72,0.25,0.85,0.38);
	leg_visM_taulepmindR -> AddEntry(h_visM_taulepmindR_sig, "signal", "L");
	leg_visM_taulepmindR -> AddEntry(h_visM_taulepmindR_TTJets, "tt+Jets", "L");
	leg_visM_taulepmindR->Draw("same");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/visM_taulepmindR.pdf");
	
}

void NtupleHistDrawer(TH1F* h_cutflow_sig, TH1F* h_cutflow_TTJets,
					  TH1F* h_njets_sig, TH1F* h_njets_TTJets,
					  TH1F* h_nbtags_sig, TH1F* h_nbtags_TTJets,
					  TH1F* h_ntauID_sig, TH1F* h_ntauID_TTJets)
{

	/// Cutflow
	//int nsample_sig = h_cutflow_sig -> GetBinContent(1);
	//int nsample_TTJets = h_cutflow_TTJets -> GetBinContent(1);	
	h_cutflow_sig->Sumw2();
	h_cutflow_TTJets->Sumw2();
	h_cutflow_sig->Scale(sf_sig);
	h_cutflow_TTJets->Scale(sf_TTJets);
	
	TCanvas *c0 = new TCanvas("c0", "", 800, 800);
	TPad *pad_1 = new TPad("pad_1", "", 0, 0.3, 1, 1.0);
	pad_1->SetBottomMargin(2);
	pad_1->SetGridx();
	pad_1->SetLogy();
	pad_1->Draw();
	pad_1->cd();
	h_cutflow_TTJets->SetStats(0);
	h_cutflow_TTJets->SetLineColor(kBlue);
	h_cutflow_TTJets->SetMarkerStyle(20);
	h_cutflow_TTJets->SetMarkerColor(kBlue);
	h_cutflow_TTJets->GetXaxis()->SetTitle("Cuts");
	h_cutflow_TTJets->GetYaxis()->SetTitle("");
	h_cutflow_TTJets->GetYaxis()->SetTitleSize(0.03);
	h_cutflow_TTJets->GetYaxis()->SetTitleOffset(0.8);
	h_cutflow_TTJets->GetYaxis()->SetLabelSize(0.02);
	h_cutflow_TTJets->SetTitle("Cut Flow (100 fb^{-1})");
	h_cutflow_TTJets->GetXaxis()->SetBinLabel(7, ">= 2 jets");
	h_cutflow_TTJets->GetXaxis()->SetBinLabel(8, ">= 1 b-tag");
	h_cutflow_TTJets->Draw();
	h_cutflow_sig->Scale(sig_scale);
	h_cutflow_sig->SetLineColor(kRed);
	h_cutflow_sig->SetMarkerStyle(20);
	h_cutflow_sig->SetMarkerColor(kRed);
	h_cutflow_sig->GetXaxis()->SetLabelSize(0.);
	h_cutflow_sig->Draw("same");

	TLegend *leg_cf = new TLegend(0.72, 0.65, 0.85, 0.78);
	string siglabel_cf = "signal (x"+to_string(sig_scale)+")";
	leg_cf->AddEntry(h_cutflow_sig, siglabel_cf.c_str(), "p");
	leg_cf->AddEntry(h_cutflow_TTJets, "tt+jets", "p");
	leg_cf->Draw("same");

	c0->cd();
	TPad *pad_2 = new TPad("pad_2", "", 0, 0.05, 1, 0.3);
	pad_2->SetTopMargin(0);
	pad_2->SetBottomMargin(0.01);
	pad_2->SetGridx();
	pad_2->SetLogy();
	pad_2->Draw();
	pad_2->cd();

	TH1F* h_cutflow_ratio = (TH1F*) h_cutflow_sig->Clone("h_cutflow_ratio");
	h_cutflow_ratio->Scale(1.0/sig_scale);
	h_cutflow_ratio->SetLineColor(kBlack);
	h_cutflow_ratio->SetStats(0);
	h_cutflow_ratio->Divide(h_cutflow_TTJets);
	h_cutflow_ratio->SetMarkerStyle(21);
	h_cutflow_ratio->SetMarkerColor(kBlack);
	h_cutflow_ratio->SetTitle("");
	h_cutflow_ratio->GetXaxis()->SetLabelSize(0.);
	h_cutflow_ratio->GetYaxis()->SetTitle("sig/bkg");
	h_cutflow_ratio->GetYaxis()->SetTitleSize(0.1);
	h_cutflow_ratio->GetYaxis()->SetLabelSize(0.08);
	h_cutflow_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_cutflow_ratio->GetYaxis()->SetNdivisions(505);
	h_cutflow_ratio->Draw("p");

	c0->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/cutflow.pdf");

	
	/// Jet multiplicity
	h_njets_sig->Sumw2();
	h_njets_TTJets->Sumw2();
	h_nbtags_sig->Sumw2();
	h_nbtags_TTJets->Sumw2();
	h_njets_sig->Scale(1.0/nsample_sig);
	h_njets_TTJets->Scale(1.0/nsample_TTJets);
	h_nbtags_sig->Scale(1.0/nsample_sig);
	h_nbtags_TTJets->Scale(1.0/nsample_TTJets);
	
	TCanvas *c1 = new TCanvas("c1", "", 800, 800);
	TPad *pad1 = new TPad("pad1", "", 0, 0.3, 1, 1.0);
	//pad1->SetTopMargin(0.06);
	pad1->SetBottomMargin(2);
	pad1->SetGridx();
	pad1->Draw();
	pad1->cd();
	h_njets_TTJets->SetStats(0);
	h_njets_TTJets->SetLineColor(kBlue);
	h_njets_TTJets->SetMarkerStyle(20);
	h_njets_TTJets->SetMarkerColor(kBlue);
	h_njets_TTJets->GetXaxis()->SetTitle("number of jets");
	//h_njets_TTJets->SetTitleOffset(0.1);
	h_njets_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_njets_TTJets->GetYaxis()->SetTitleSize(0.04);
	h_njets_TTJets->GetYaxis()->SetTitleOffset(0.8);
	h_njets_TTJets->GetYaxis()->SetLabelSize(0.03);
	h_njets_TTJets->SetTitle("SS leptons + #tau_{h}");
	h_njets_TTJets->Draw();
	h_njets_sig->SetLineColor(kRed);
	h_njets_sig->SetMarkerStyle(20);
	h_njets_sig->SetMarkerColor(kRed);;
	h_njets_sig->Draw("same");
	
	TLegend *leg_njets = new TLegend(0.72,0.25,0.85,0.38);
	leg_njets->AddEntry(h_njets_sig, "signal","p");
	leg_njets->AddEntry(h_njets_TTJets, "tt+Jets","p");
	leg_njets->Draw("same");

	c1->cd();
	TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.01);
	pad2->SetGridx();
	pad2->Draw();
	pad2->cd();

	TH1F* h_njets_ratio = (TH1F*) h_njets_sig->Clone("h_njets_ratio");
	h_njets_ratio->SetLineColor(kBlack);
	h_njets_ratio->Sumw2();
	h_njets_ratio->SetStats(0);
	h_njets_ratio->Divide(h_njets_TTJets);
	h_njets_ratio->SetMarkerStyle(21);
	h_njets_ratio->SetMarkerColor(kBlack);
	h_njets_ratio->SetTitle("");
	h_njets_ratio->GetXaxis()->SetLabelSize(0.);
	h_njets_ratio->GetYaxis()->SetTitle("sig/bkg");
	h_njets_ratio->GetYaxis()->SetTitleSize(0.1);
	h_njets_ratio->GetYaxis()->SetLabelSize(0.08);
	h_njets_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_njets_ratio->GetYaxis()->SetNdivisions(505);
	h_njets_ratio->Draw("ep");

	c1->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/njets.pdf");
	
	// btags
	TCanvas *c2 = new TCanvas("c2", "", 800, 800);
	TPad *pad3 = new TPad("pad3", "", 0, 0.3, 1, 1.0);
	//pad3->SetTopMargin(0.06);
	pad3->SetBottomMargin(2);
	pad3->SetGridx();
	pad3->Draw();
	pad3->cd();
	h_nbtags_TTJets->SetStats(0);
	h_nbtags_TTJets->SetLineColor(kBlue);
	h_nbtags_TTJets->SetMarkerStyle(20);
	h_nbtags_TTJets->SetMarkerColor(kBlue);
	h_nbtags_TTJets->GetXaxis()->SetTitle("number of btags");
	//h_nbtags_TTJets->SetTitleOffset(0.1);
	h_nbtags_TTJets->GetYaxis()->SetTitle("Normalized to 1");
	h_nbtags_TTJets->GetYaxis()->SetTitleSize(0.04);
	h_nbtags_TTJets->GetYaxis()->SetTitleOffset(0.8);
	h_nbtags_TTJets->GetYaxis()->SetLabelSize(0.03);
	h_nbtags_TTJets->SetTitle("SS leptons + #tau_{h}");
	h_nbtags_TTJets->Draw();
	h_nbtags_sig->SetLineColor(kRed);
	h_nbtags_sig->SetMarkerStyle(20);
	h_nbtags_sig->SetMarkerColor(kRed);
	h_nbtags_sig->Draw("same");

	TLegend *leg_nbtags = new TLegend(0.72,0.25,0.85,0.38);
	leg_nbtags->AddEntry(h_nbtags_sig, "signal","p");
	leg_nbtags->AddEntry(h_nbtags_TTJets, "tt+Jets","p");
	leg_nbtags->Draw("same");

	c2->cd();
	TPad *pad4 = new TPad("pad4", "", 0, 0.05, 1, 0.3);
	pad4->SetTopMargin(0);
	pad4->SetBottomMargin(0.01);
	pad4->SetGridx();
	pad4->Draw();
	pad4->cd();

	TH1F* h_nbtags_ratio = (TH1F*) h_nbtags_sig->Clone("h_nbtags_ratio");
	h_nbtags_ratio->SetLineColor(kBlack);
	h_nbtags_ratio->Sumw2();
	h_nbtags_ratio->SetStats(0);
	h_nbtags_ratio->Divide(h_nbtags_TTJets);
	h_nbtags_ratio->SetMarkerStyle(21);
	h_nbtags_ratio->SetMarkerColor(kBlack);
	h_nbtags_ratio->SetTitle("");
	h_nbtags_ratio->GetXaxis()->SetLabelSize(0.);
	h_nbtags_ratio->GetYaxis()->SetTitle("sig/bkg");
	h_nbtags_ratio->GetYaxis()->SetTitleSize(0.1);
	h_nbtags_ratio->GetYaxis()->SetLabelSize(0.08);
	h_nbtags_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_nbtags_ratio->GetYaxis()->SetNdivisions(505);
	h_nbtags_ratio->Draw("ep");

	c2->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/nbtags.pdf");

	/// Tau ID
	h_ntauID_sig->Sumw2();
	h_ntauID_TTJets->Sumw2();
	h_ntauID_sig->Scale(1.0/nsample_sig);
	h_ntauID_TTJets->Scale(1.0/nsample_TTJets);
	
	TCanvas *c3 = new TCanvas("c3", "", 800, 800);
	TPad *pad5 = new TPad("pad5", "", 0, 0.3, 1, 1.0);
	//pad5->SetTopMargin(0.06);
	pad5->SetBottomMargin(2);
	pad5->SetGridx();
	pad5->Draw();
	pad5->cd();
	h_ntauID_TTJets->SetStats(0);
	h_ntauID_TTJets->SetLineColor(kRed);
	h_ntauID_TTJets->SetMarkerStyle(20);
	h_ntauID_TTJets->SetMarkerColor(kRed);
	h_ntauID_TTJets->GetXaxis()->SetTitle("TauID WP");
	//h_ntauID_TTJets->GetXaxis()->SetBinLabel(1,"NonIso");
	//h_ntauID_TTJets->GetXaxis()->SetBinLabel(2,"Loose");
	//h_ntauID_TTJets->GetXaxis()->SetBinLabel(3,"Medium");
	//h_ntauID_TTJets->GetXaxis()->SetBinLabel(4,"Tight");
	//h_ntauID_TTJets->SetTitleOffset(0.1);
	h_ntauID_TTJets->GetYaxis()->SetTitle("Normalized");
	h_ntauID_TTJets->GetYaxis()->SetTitleSize(0.04);
	h_ntauID_TTJets->GetYaxis()->SetTitleOffset(0.8);
	h_ntauID_TTJets->GetYaxis()->SetLabelSize(0.03);
	h_ntauID_TTJets->SetTitle("SS leptons + #tau_{h} (>=1 tau)");
	h_ntauID_TTJets->SetMinimum(0);
	h_ntauID_TTJets->SetMaximum(1.1);
	h_ntauID_TTJets->Draw();
	h_ntauID_sig->SetLineColor(kBlue);
	h_ntauID_sig->SetMarkerStyle(20);
	h_ntauID_sig->SetMarkerColor(kBlue);
	h_ntauID_sig->Draw("same");

	TLegend *leg_ntauID = new TLegend(0.72,0.27,0.85,0.40);
	leg_ntauID->AddEntry(h_ntauID_sig, "signal","p");
	leg_ntauID->AddEntry(h_ntauID_TTJets, "tt+Jets","p");
	leg_ntauID->Draw("same");

	c3->cd();
	TPad *pad6 = new TPad("pad6", "", 0, 0.05, 1, 0.3);
	pad6->SetTopMargin(0);
	pad6->SetBottomMargin(0.01);
	pad6->SetGridx();
	pad6->Draw();
	pad6->cd();

	TH1F* h_ntauID_ratio = (TH1F*) h_ntauID_sig->Clone("h_ntauID_ratio");
	h_ntauID_ratio->SetLineColor(kBlack);
	h_ntauID_ratio->Sumw2();
	h_ntauID_ratio->SetStats(0);
	h_ntauID_ratio->Divide(h_ntauID_TTJets);
	h_ntauID_ratio->SetMarkerStyle(21);
	h_ntauID_ratio->SetMarkerColor(kBlack);
	h_ntauID_ratio->SetTitle("");
	h_ntauID_ratio->GetXaxis()->SetLabelSize(0.);
	h_ntauID_ratio->GetYaxis()->SetTitle("sig/bkg");
	h_ntauID_ratio->GetYaxis()->SetTitleSize(0.1);
	h_ntauID_ratio->GetYaxis()->SetLabelSize(0.08);
	h_ntauID_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_ntauID_ratio->GetYaxis()->SetNdivisions(505);
	h_ntauID_ratio->Draw("ep");

	c3->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/ntauID.pdf");
	
	/*
	/// cut on number of jets
	TH1F* h_njetscut = new TH1F("h_njetscut", "", 9, -0.5, 8.5);
	TH1F* h_nbtagscut = new TH1F("h_nbtagscut", "", 5, -0.5, 4.5);

	TCanvas *c4 = new TCanvas("c4", "", 800, 800);
	TPad *pad7 = new TPad("pad7", "", 0, 0.3, 1, 1.0);
	//pad7->SetTopMargin(0.06);
	pad7->SetBottomMargin(2);
	pad7->SetGridx();
	pad7->Draw();
	pad7->cd();
	h_njetscut_sig->SetStats(0);
	h_njetscut_sig->SetLineColor(kRed);
	h_njetscut_sig->SetMarkerStyle(20);
	h_njetscut_sig->SetMarkerColor(kRed);
	h_njetscut_sig->GetXaxis()->SetTitle(">= n_jets");
	//h_njetscut_sig->SetTitleOffset(0.1);
	h_njetscut_sig->GetYaxis()->SetTitle("(Normalized)");
	h_njetscut_sig->GetYaxis()->SetTitleSize(0.04);
	h_njetscut_sig->GetYaxis()->SetTitleOffset(0.8);
	h_njetscut_sig->GetYaxis()->SetLabelSize(0.03);
	h_njetscut_sig->SetTitle("Cut on number of jets");
	h_njetscut_sig->Draw();
	h_njetscut_TTJets->SetLineColor(kBlue);
	h_njetscut_TTJets->SetMarkerStyle(20);
	h_njetscut_TTJets->SetMarkerColor(kBlue);
	h_njetscut_TTJets->Draw("same");

	TLegend *leg_njetscut = new TLegend(0.72,0.25,0.85,0.38);
	leg_njetscut->AddEntry(h_njetscut_sig, "signal","p");
	leg_njetscut->AddEntry(h_njetscut_TTJets, "tt+Jets","p");
	leg_njetscut->Draw("same");

	c4->cd();
	TPad *pad8 = new TPad("pad8", "", 0, 0.05, 1, 0.3);
	pad8->SetTopMargin(0);
	pad8->SetBottomMargin(0.01);
	pad8->SetGridx();
	pad8->Draw();
	pad8->cd();

	TH1F* h_njetscut_ratio = (TH1F*) h_njetscut_TTJets->Clone("h_njetscut_ratio");
	h_njetscut_ratio->SetLineColor(kBlack);
	h_njetscut_ratio->Sumw2();
	h_njetscut_ratio->SetStats(0);
	h_njetscut_ratio->Divide(h_njetscut_sig);
	h_njetscut_ratio->SetMarkerStyle(21);
	h_njetscut_ratio->SetMarkerColor(kBlack);
	h_njetscut_ratio->SetTitle("");
	h_njetscut_ratio->GetXaxis()->SetLabelSize(0.);
	h_njetscut_ratio->GetYaxis()->SetTitle("TTJets/Signal");
	h_njetscut_ratio->GetYaxis()->SetTitleSize(0.1);
	h_njetscut_ratio->GetYaxis()->SetLabelSize(0.08);
	h_njetscut_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_njetscut_ratio->GetYaxis()->SetNdivisions(505);
	h_njetscut_ratio->Draw("ep");

	c4->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/njetscut.pdf");
	
	/// Cut on number of b-jets
	TCanvas *c5 = new TCanvas("c5", "", 800, 800);
	TPad *pad9 = new TPad("pad9", "", 0, 0.3, 1, 1.0);
	//pad9->SetTopMargin(0.06);
	pad9->SetBottomMargin(2);
	pad9->SetGridx();
	pad9->Draw();
	pad9->cd();
	h_nbtagscut_sig->SetStats(0);
	h_nbtagscut_sig->SetLineColor(kRed);
	h_nbtagscut_sig->SetMarkerStyle(20);
	h_nbtagscut_sig->SetMarkerColor(kRed);
	h_nbtagscut_sig->GetXaxis()->SetTitle(">= n_btags");
	//h_nbtagscut_sig->SetTitleOffset(0.1);
	h_nbtagscut_sig->GetYaxis()->SetTitle("(Normalized)");
	h_nbtagscut_sig->GetYaxis()->SetTitleSize(0.04);
	h_nbtagscut_sig->GetYaxis()->SetTitleOffset(0.8);
	h_nbtagscut_sig->GetYaxis()->SetLabelSize(0.03);
	h_nbtagscut_sig->SetTitle("Cut on number of b-tagged jets");
	h_nbtagscut_sig->Draw();
	h_nbtagscut_TTJets->SetLineColor(kBlue);
	h_nbtagscut_TTJets->SetMarkerStyle(20);
	h_nbtagscut_TTJets->SetMarkerColor(kBlue);
	h_nbtagscut_TTJets->Draw("same");

	TLegend *leg_nbtagscut = new TLegend(0.72,0.25,0.85,0.38);
	leg_nbtagscut->AddEntry(h_nbtagscut_sig, "signal","p");
	leg_nbtagscut->AddEntry(h_nbtagscut_TTJets, "tt+Jets","p");
	leg_nbtagscut->Draw("same");

	c5->cd();
	TPad *pad10 = new TPad("pad10", "", 0, 0.05, 1, 0.3);
	pad10->SetTopMargin(0);
	pad10->SetBottomMargin(0.01);
	pad10->SetGridx();
	pad10->Draw();
	pad10->cd();

	TH1F* h_nbtagscut_ratio = (TH1F*) h_nbtagscut_TTJets->Clone("h_nbtagscut_ratio");
	h_nbtagscut_ratio->SetLineColor(kBlack);
	h_nbtagscut_ratio->Sumw2();
	h_nbtagscut_ratio->SetStats(0);
	h_nbtagscut_ratio->Divide(h_nbtagscut_sig);
	h_nbtagscut_ratio->SetMarkerStyle(21);
	h_nbtagscut_ratio->SetMarkerColor(kBlack);
	h_nbtagscut_ratio->SetTitle("");
	h_nbtagscut_ratio->GetXaxis()->SetLabelSize(0.);
	h_nbtagscut_ratio->GetYaxis()->SetTitle("TTJets/Signal");
	h_nbtagscut_ratio->GetYaxis()->SetTitleSize(0.1);
	h_nbtagscut_ratio->GetYaxis()->SetLabelSize(0.08);
	h_nbtagscut_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_nbtagscut_ratio->GetYaxis()->SetNdivisions(505);
	h_nbtagscut_ratio->Draw("ep");

	c5->SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/nbtagscut.pdf");
	*/
}
