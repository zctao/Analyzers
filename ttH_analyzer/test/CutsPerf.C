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
void Draw_Histogram_W2Stat(TH1F* h_sig, TH1F* h_bkg, TString pname, TString xtitle, TString ytitle, TString title, bool setlogy);
void Draw_Histogram_WRatio(TH1F* h_sig, TH1F* h_bkg, TString xtitle, TString ytitle, TString title);

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
			  const TString sig_file = "/eos/uscms/store/user/ztao/ttHToTT_M125_13TeV_powheg_pythia8/ttHToTauTau_Ntuple_signal/151209_202153/0000/CU_ttH_EDA_output_tmp.root",
			  const TString bkg_file = "/eos/uscms/store/user/ztao/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ttHToTauTau_Ntuple_TTJets/151209_202324/0000/CU_ttH_EDA_output_tmp.root"
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
	float metsig = -999.9;
	float metcov00 = -9999.9;
	float metcov01 = -9999.9;  // metcov10 = metcov01
	float metcov11 = -9999.9;

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
	TBranch* b_metsig;
	TBranch* b_metcov00;
	TBranch* b_metcov01;
	TBranch* b_metcov11;

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
	tree->SetBranchAddress("METSignificance", &metsig, &b_metsig);
	tree->SetBranchAddress("METCovariance00", &metcov00, &b_metcov00);
	tree->SetBranchAddress("METCovariance01", &metcov01, &b_metcov01);
	tree->SetBranchAddress("METCovariance11", &metcov11, &b_metcov11);

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

	TH1F* h_MET = new TH1F("h_MET", "", 50, 0 ,500);
	TH1F* h_MHT = new TH1F("h_MHT", "", 50, 0, 500);
	TH1F* h_METSig = new TH1F("h_METSig", "", 40, 0, 200);
	TH1F* h_HT = new TH1F("h_HT", "", 100, 0, 1000);
	TH1F* h_metLD = new TH1F("h_metLD", "", 50, 0, 5);
	TH1F* h_visMtot = new TH1F("h_visMtot", "", 50, 0, 2000);
	TH1F* h_dilepmass = new TH1F("h_dilepmass", "", 20, 0, 100);

	TH1F* h_dRlep1jet1 = new TH1F("h_dRlep1jet1", "", 25, 0., 5.);
	TH1F* h_dRlep1jet2 = new TH1F("h_dRlep1jet2", "", 25, 0., 5.);
	TH1F* h_dRlep2jet1 = new TH1F("h_dRlep2jet1", "", 25, 0., 5.);
	TH1F* h_dRlep2jet2 = new TH1F("h_dRlep2jet2", "", 25, 0., 5.);
	
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

		// MET
		float MET = TMath::Sqrt(met_x*met_x+met_y*met_y);
		h_MET -> Fill(MET);
		// MET Significance
		h_METSig -> Fill(metsig);
		
		/// Calculate variables that are not yet included in the ntuple
		std::vector<TLorentzVector> leptons;
		std::vector<TLorentzVector> jets;
		std::vector<TLorentzVector> tau_soft; // pT < jet pT cut (30)
		std::vector<TLorentzVector> b_soft; // pT < jet pT cut (30)

		for (size_t itr=0; itr < lep_pt->size(); ++itr) {
			TLorentzVector lepp4;
			lepp4.SetPtEtaPhiM(lep_pt->at(itr),lep_eta->at(itr),lep_phi->at(itr),lep_mass->at(itr));
			leptons.push_back(lepp4);
		}
		for (size_t itr=0; itr < jet_pt->size(); ++itr) {
			TLorentzVector jetp4;
			jetp4.SetPtEtaPhiM(jet_pt->at(itr),jet_eta->at(itr),jet_phi->at(itr),jet_mass->at(itr));
			jets.push_back(jetp4);
		}

		// bjet and tau pT cuts are smaller than jet pT cut, include tau and bjet with pT < jet pT cut(30)?
		// pT cuts for tau and bjet are 20 GeV
		TLorentzVector bp4, taup4;
		for (size_t ib=0; ib < bjet_pt->size(); ++ib) {
			if (bjet_pt->at(ib) < 30) {
				bp4.SetPtEtaPhiM(bjet_pt->at(ib),bjet_eta->at(ib),
								 bjet_phi->at(ib),bjet_mass->at(ib));
				b_soft.push_back(bp4);
			}
		}
		for (size_t itau=0; itau < loose_tau_pt->size(); ++itau) {
			if (loose_tau_pt->at(itau) < 30) {
				taup4.SetPtEtaPhiM(loose_tau_pt->at(itau),loose_tau_eta->at(itau),
								 loose_tau_phi->at(itau),loose_tau_mass->at(itau));
				tau_soft.push_back(taup4);
			}
		}

		// upper bound on number of jets?
		
		// HT
		float HT = 0;
		for (auto & pt : *lep_pt)
			HT += pt;
		for (auto & pt : *jet_pt)
			HT += pt;
		for (auto & pt : *loose_tau_pt) {
			if (pt < 30)
				HT += pt;
		}
		for (auto & pt : *bjet_pt) {
			if (pt < 30)
				HT += pt;
		}
		h_HT -> Fill(HT);
		
		// MHT
		float MHT_x = 0;
		float MHT_y = 0;
		for (auto & l : leptons) {
			MHT_x -= l.Px();
			MHT_y -= l.Py();
		}
		for (auto & j : jets) {
			MHT_x -= j.Px();
			MHT_y -= j.Py();
		}
		for (auto & tau: tau_soft) {
			MHT_x -= tau.Px();
			MHT_y -= tau.Py();
		}
		for (auto & b : b_soft) {
			MHT_x -= b.Px();
			MHT_y -= b.Py();
		}
		float MHT = TMath::Sqrt(MHT_x*MHT_x+MHT_y*MHT_y);
		h_MHT -> Fill(MHT);

		// metLD
		float metLD = 0.00397*MET+ 0.00265*MHT;
		h_metLD -> Fill(metLD);
		
		// total visible mass
		TLorentzVector p4tot;
		for (auto & l : leptons)
			p4tot += l;
		for (auto & j : jets)
			p4tot += j;
		for (auto & tau : tau_soft)
			p4tot += tau;
		for (auto & b : b_soft)
			p4tot += b;
		
		h_visMtot -> Fill(p4tot.M());
		
		/// Additional cuts if applicale
		
		++nevtpass;
		
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
	
		
		float dz_lep1 = lep_vtx_dz->at(ilep1);
		float dz_lep2 = lep_vtx_dz->at(ilep2);
		float dxy_lep1 = lep_vtx_dxy->at(ilep1);
		float dxy_lep2 = lep_vtx_dxy->at(ilep2);

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

		// Di-lepton mass
		if ((lep1_p4+lep2_p4).M() < 100)
			h_dilepmass -> Fill((lep1_p4+lep2_p4).M());

		// dR between leptons and jets
		h_dRlep1jet1 -> Fill(lep1_p4.DeltaR(jets[0]));
		h_dRlep1jet2 -> Fill(lep1_p4.DeltaR(jets[1]));
		h_dRlep2jet1 -> Fill(lep2_p4.DeltaR(jets[0]));
		h_dRlep2jet2 -> Fill(lep2_p4.DeltaR(jets[1]));
			
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

	h_MET -> Write();
	h_METSig -> Write();
	h_MHT -> Write();
	h_HT -> Write();
	h_metLD -> Write();
	h_visMtot -> Write();
	h_dilepmass -> Write();

	h_dRlep1jet1 -> Write();
	h_dRlep1jet2 -> Write();
	h_dRlep2jet1 -> Write();
	h_dRlep2jet2 -> Write();
	
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
	TH1F* h_MET_sig = (TH1F*)f_sig->Get("h_MET");
	TH1F* h_METSig_sig = (TH1F*)f_sig->Get("h_METSig");
	TH1F* h_MHT_sig = (TH1F*)f_sig->Get("h_MHT");
	TH1F* h_HT_sig = (TH1F*)f_sig->Get("h_HT");
	TH1F* h_metLD_sig = (TH1F*)f_sig->Get("h_metLD");
	TH1F* h_visMtot_sig = (TH1F*)f_sig->Get("h_visMtot");
	TH1F* h_dilepmass_sig = (TH1F*)f_sig->Get("h_dilepmass");
	TH1F* h_dRlep1jet1_sig = (TH1F*)f_sig->Get("h_dRlep1jet1");
	TH1F* h_dRlep1jet2_sig = (TH1F*)f_sig->Get("h_dRlep1jet2");
	TH1F* h_dRlep2jet1_sig = (TH1F*)f_sig->Get("h_dRlep2jet1");
	TH1F* h_dRlep2jet2_sig = (TH1F*)f_sig->Get("h_dRlep2jet2");
	
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
	TH1F* h_MET_TTJets = (TH1F*)f_TTJets->Get("h_MET");
	TH1F* h_METSig_TTJets = (TH1F*)f_TTJets->Get("h_METSig");
	TH1F* h_MHT_TTJets = (TH1F*)f_TTJets->Get("h_MHT");
	TH1F* h_HT_TTJets = (TH1F*)f_TTJets->Get("h_HT");
	TH1F* h_metLD_TTJets = (TH1F*)f_TTJets->Get("h_metLD");
	TH1F* h_visMtot_TTJets = (TH1F*)f_TTJets->Get("h_visMtot");
	TH1F* h_dilepmass_TTJets = (TH1F*)f_TTJets->Get("h_dilepmass");
	TH1F* h_dRlep1jet1_TTJets = (TH1F*)f_TTJets->Get("h_dRlep1jet1");
	TH1F* h_dRlep1jet2_TTJets = (TH1F*)f_TTJets->Get("h_dRlep1jet2");
	TH1F* h_dRlep2jet1_TTJets = (TH1F*)f_TTJets->Get("h_dRlep2jet1");
	TH1F* h_dRlep2jet2_TTJets = (TH1F*)f_TTJets->Get("h_dRlep2jet2");
	
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
	h_MET_sig -> Scale(1.0/nsig);
	h_METSig_sig -> Scale(1.0/nsig);
	h_MHT_sig -> Scale(1.0/nsig);
	h_HT_sig -> Scale(1.0/nsig);
	h_metLD_sig -> Scale(1.0/nsig);
	h_visMtot_sig -> Scale(1.0/nsig);
	h_dilepmass_sig -> Scale(1.0/nsig);
	h_dRlep1jet1_sig -> Scale(1.0/nsig);
	h_dRlep1jet2_sig -> Scale(1.0/nsig);
	h_dRlep2jet1_sig -> Scale(1.0/nsig);
	h_dRlep2jet2_sig -> Scale(1.0/nsig);
	
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
	h_MET_TTJets -> Scale(1.0/nTTJets);
	h_METSig_TTJets -> Scale(1.0/nTTJets);
	h_MHT_TTJets -> Scale(1.0/nTTJets);
	h_HT_TTJets -> Scale(1.0/nTTJets);
	h_metLD_TTJets -> Scale(1.0/nTTJets);
	h_visMtot_TTJets -> Scale(1.0/nTTJets);
	h_dilepmass_TTJets -> Scale(1.0/nTTJets);
	h_dRlep1jet1_TTJets -> Scale(1.0/nTTJets);
	h_dRlep1jet2_TTJets -> Scale(1.0/nTTJets);
	h_dRlep2jet1_TTJets -> Scale(1.0/nTTJets);
	h_dRlep2jet2_TTJets -> Scale(1.0/nTTJets);

	// Plotting
	//// h_lep1pt
	Draw_Histogram_W2Stat(h_lep1pt_sig, h_lep1pt_TTJets, "lep1pt",
						  "Leading Lepton p_{T} [GeV]",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_lep1eta
	Draw_Histogram_W2Stat(h_lep1eta_sig, h_lep1eta_TTJets, "lep1eta",
						  "Leading Lepton #eta",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_lep2pt
	Draw_Histogram_W2Stat(h_lep2pt_sig, h_lep2pt_TTJets, "lep2pt",
						  "Subleading Lepton p_{T} [GeV]",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_lep2eta
	Draw_Histogram_W2Stat(h_lep2eta_sig, h_lep2eta_TTJets, "lep2eta",
						  "Subleading Lepton #eta",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_taupt
	Draw_Histogram_W2Stat(h_taupt_sig, h_taupt_TTJets, "taupt",
						  "#tau p_{T} [GeV]",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_taueta
	Draw_Histogram_W2Stat(h_taueta_sig, h_taueta_TTJets, "taueta",
						  "#tau #eta",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_bjet1pt
	Draw_Histogram_W2Stat(h_bjet1pt_sig, h_bjet1pt_TTJets, "bjet1pt",
						  "Leading b-jet p_{T} [GeV]",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_bjet1eta
	Draw_Histogram_W2Stat(h_bjet1eta_sig, h_bjet1eta_TTJets, "bjet1eta",
						  "Leading b-jet #eta",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_dRlep1tau
	Draw_Histogram_W2Stat(h_dRlep1tau_sig, h_dRlep1tau_TTJets, "dRlep1tau",
						  "#DeltaR(l_{ldg},#tau_{h})",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_dRlep2tau
	Draw_Histogram_W2Stat(h_dRlep2tau_sig, h_dRlep2tau_TTJets, "dRlep2tau",
						  "#DeltaR(l_{sub-ldg},#tau_{h})",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_lep1dxy
	Draw_Histogram_W2Stat(h_lep1dxy_sig, h_lep1dxy_TTJets, "lep1dxy",
						  "dxy_{ldg lep}",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_lep2dxy
	Draw_Histogram_W2Stat(h_lep2dxy_sig, h_lep2dxy_TTJets, "lep2dxy",
						  "dxy_{sub-ldg lep}",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_lep1dz
	Draw_Histogram_W2Stat(h_lep1dz_sig, h_lep1dz_TTJets, "lep1dz",
						  "dz_{ldg lep}",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_lep2dz
	Draw_Histogram_W2Stat(h_lep2dz_sig, h_lep2dz_TTJets, "lep2dz",
						  "dz_{sub-ldg lep}",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_lepmaxdxy
	Draw_Histogram_W2Stat(h_lepmaxdxy_sig, h_lepmaxdxy_TTJets, "lepmaxdxy",
						  "max(dxy_{ldg lep},dxy_{sub-ldg lep})",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_visM_taulep1
	Draw_Histogram_W2Stat(h_visM_taulep1_sig, h_visM_taulep1_TTJets, "visM_taulep1",
						  "visMass(#tau_{h},l_{ldg})",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_visM_taulep2
	Draw_Histogram_W2Stat(h_visM_taulep2_sig, h_visM_taulep2_TTJets, "visM_taulep2",
						  "visMass(#tau_{h},l_{subldg})",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_visM_taulepmaxdxy
	Draw_Histogram_W2Stat(h_visM_taulepmaxdxy_sig, h_visM_taulepmaxdxy_TTJets, "visM_taulepmaxdxy",
						  "visMass(#tau_{h},l_{max(dxy)})",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_visM_taulepmindR
	Draw_Histogram_W2Stat(h_visM_taulepmindR_sig, h_visM_taulepmindR_TTJets, "visM_taulepmindR",
						  "visMass(#tau_{h},l_{min(dR)})",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_MET
	Draw_Histogram_W2Stat(h_MET_sig, h_MET_TTJets, "MET",
						  "MET [GeV]",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	// h_METSig
	Draw_Histogram_W2Stat(h_METSig_sig, h_METSig_TTJets, "METSig_log",
						  "MET Significance [GeV]",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", true);
	Draw_Histogram_W2Stat(h_METSig_sig, h_METSig_TTJets, "METSig",
						  "MET Significance [GeV]",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_MHT
	Draw_Histogram_W2Stat(h_MHT_sig, h_MHT_TTJets, "MHT",
						  "MHT [GeV]",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_HT
	Draw_Histogram_W2Stat(h_HT_sig, h_HT_TTJets, "HT",
						  "HT [GeV]",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_metLD
	Draw_Histogram_W2Stat(h_metLD_sig, h_metLD_TTJets, "metLD",
						  "metLD [GeV]",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_visMtot
	Draw_Histogram_W2Stat(h_visMtot_sig, h_visMtot_TTJets, "visMtot",
						  "visMtot [GeV]",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	//// h_dilepmass
	Draw_Histogram_W2Stat(h_dilepmass_sig, h_dilepmass_TTJets, "dilepmass",
						  "dilepmass [GeV] (< 100 GeV)",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	// dR between lepton and jet
	Draw_Histogram_W2Stat(h_dRlep1jet1_sig, h_dRlep1jet1_TTJets, "dRlep1jet1",
						  "#DeltaR(l_{ldg},jet_{ldg})",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	Draw_Histogram_W2Stat(h_dRlep1jet2_sig, h_dRlep1jet2_TTJets, "dRlep1jet2",
						  "#DeltaR(l_{ldg},jet_{subldg})",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	Draw_Histogram_W2Stat(h_dRlep2jet1_sig, h_dRlep2jet1_TTJets, "dRlep2jet1",
						  "#DeltaR(l_{subldg},jet_{ldg})",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	Draw_Histogram_W2Stat(h_dRlep2jet2_sig, h_dRlep2jet2_TTJets, "dRlep2jet2",
						  "#DeltaR(l_{subldg},jet_{subldg})",
						  "Normalized to 1",
						  "same sign leptons + #tau_{h}", false);
	
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
	h_cutflow_TTJets->GetXaxis()->SetBinLabel(2, "HLT");;
	h_cutflow_TTJets->GetXaxis()->SetBinLabel(7, ">= 2 jets");
	h_cutflow_TTJets->GetXaxis()->SetBinLabel(8, "n_ btags");
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
	h_nbtags_TTJets->GetXaxis()->SetTitle("number of btags (loose)");
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

void Draw_Histogram_W2Stat(TH1F* h_sig, TH1F* h_bkg, TString pname, TString xtitle, TString ytitle, TString title, bool setlogy)
{
	TCanvas c;

	if (setlogy)
		gPad->SetLogy(1);
	else
		gPad->SetLogy(0);
	
	h_sig->SetTitle(title);
	h_sig->GetXaxis()->SetTitle(xtitle);
	h_sig->GetYaxis()->SetTitle(ytitle);
	h_sig->SetLineColor(kRed);
	float ymax = max(h_sig->GetMaximum(), h_bkg->GetMaximum());
	h_sig->SetMaximum(ymax*1.1);
	h_sig->Draw("");
	gPad->Update();
	TPaveStats *st_sig = (TPaveStats*) h_sig->GetListOfFunctions()->FindObject("stats");
	st_sig->SetTextColor(kRed);
	st_sig->SetOptStat(1110);

	h_bkg->SetLineColor(kBlue);
	h_bkg->Draw("sames");
	gPad->Update();
	TPaveStats *st_bkg = (TPaveStats*) h_bkg->GetListOfFunctions()->FindObject("stats");
	st_bkg->SetTextColor(kBlue);
	st_bkg->SetOptStat(1110);
	st_bkg->SetY1NDC(2*(st_sig->GetY1NDC())-st_sig->GetY2NDC());
	st_bkg->SetY2NDC(st_sig->GetY1NDC());
	gPad->Modified();

	TLegend *leg = new TLegend(0.72,0.25,0.85,0.38);
	leg -> AddEntry(h_sig, "signal", "L");
	leg -> AddEntry(h_bkg, "tt_jets", "L");
	leg -> Draw("same");

	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/"+pname+".pdf");
}
