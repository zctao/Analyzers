#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1.h"

#include <iostream>
#include <vector>
#include <algorithm>

#include "../macro/MyPlottingUtils.h"

using namespace std;

void fillHistofromTree(TTree*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,
					   TH1D*,TH1D*,TH1D*,TH1D*, TH1D*);

void treeAnalyzer
(
 //vector<TString> samples = {"ttH_htt","ttH_hww","ttH_hzz", "TTW", "TTZ", "TTJets"},
 vector<TString> samples = {"ttH", "TTW", "TTZ", "TTJets"},
 TString directory = "/Users/ztao/Documents/ttH/Outputs/80X/"
 )
{
	
	vector<TFile*> files;
	vector<TTree*> trees;

	// define histograms
	vector<TH1D*> hists_MVA_2lss_ttV;
	vector<TH1D*> hists_MVA_2lss_ttbar;
	vector<TH1D*> hists_MT_met_lep0;
	vector<TH1D*> hists_mindr_lep0_jet;
	vector<TH1D*> hists_mindr_lep1_jet;
	vector<TH1D*> hists_lep0_conept;
	vector<TH1D*> hists_lep1_conept;
	vector<TH1D*> hists_avg_dr_jet;

	vector<TH1D*> hists_tau_decaymode;
	vector<TH1D*> hists_dr_lep0_tau;
	vector<TH1D*> hists_dr_lep1_tau;
	vector<TH1D*> hists_mass_lep0_tau;
	vector<TH1D*> hists_mass_lep1_tau;

	for (auto& sn : samples) {
		// Read ntuples
		TFile* f = new TFile(directory+"output_"+ sn +".root");
		TTree* tree = (TTree*) f->Get("ttHtaus/eventTree");
		files.push_back(f);
		trees.push_back(tree);

		// Setup histograms
		TH1D* h_MVA_2lss_ttV = new TH1D("h_MVA_2lss_ttV", "", 10, -1.0, 1.0);
		TH1D* h_MVA_2lss_ttbar = new TH1D("h_MVA_2lss_ttbar", "", 10, -1.0, 1.0);
		TH1D* h_MT_met_lep0 = new TH1D("h_MT_met_lep0", "", 10, 0, 400);
		TH1D* h_mindr_lep0_jet = new TH1D("h_mindr_lep0_jet", "", 10, 0, 4);
		TH1D* h_mindr_lep1_jet = new TH1D("h_mindr_lep1_jet", "", 10, 0, 4);
		TH1D* h_lep0_conept = new TH1D("h_lep0_conept", "", 10, 0, 250);
		TH1D* h_lep1_conept = new TH1D("h_lep1_conept", "", 10, 0, 150);
		TH1D* h_avg_dr_jet = new TH1D("h_avg_dr_jet", "", 10, 0, 4);
		TH1D* h_tau_decaymode = new TH1D("h_tau_decaymode", "", 18, 0, 18);
		TH1D* h_dr_lep0_tau = new TH1D("h_dr_lep0_tau", "", 10, 0., 4.);
		TH1D* h_dr_lep1_tau = new TH1D("h_dr_lep1_tau", "", 10, 0., 4.);
		TH1D* h_mass_lep0_tau = new TH1D("h_mass_lep0_tau", "", 10, 0., 400.);
		TH1D* h_mass_lep1_tau = new TH1D("h_mass_lep1_tau", "", 10, 0., 400.);

		hists_MVA_2lss_ttV.push_back(h_MVA_2lss_ttV);
		hists_MVA_2lss_ttbar.push_back(h_MVA_2lss_ttbar);
		hists_MT_met_lep0.push_back(h_MT_met_lep0);
		hists_mindr_lep0_jet.push_back(h_mindr_lep0_jet);
		hists_mindr_lep1_jet.push_back(h_mindr_lep1_jet);
		hists_lep0_conept.push_back(h_lep0_conept);
		hists_lep1_conept.push_back(h_lep1_conept);
		hists_avg_dr_jet.push_back(h_avg_dr_jet);
		hists_tau_decaymode.push_back(h_tau_decaymode);
		hists_dr_lep0_tau.push_back(h_dr_lep0_tau);
		hists_dr_lep1_tau.push_back(h_dr_lep1_tau);
		hists_mass_lep0_tau.push_back(h_mass_lep0_tau);
		hists_mass_lep1_tau.push_back(h_mass_lep1_tau);
				
	}


	// Fill histograms
	int nchannels = samples.size();
	
	for (int i = 0; i < nchannels; i++) {
		fillHistofromTree(trees[i],
						  hists_MVA_2lss_ttV[i],
						  hists_MVA_2lss_ttbar[i],
						  hists_MT_met_lep0[i],
						  hists_avg_dr_jet[i],
						  hists_mindr_lep0_jet[i],
						  hists_mindr_lep1_jet[i],
						  hists_lep0_conept[i],
						  hists_lep1_conept[i],
						  hists_tau_decaymode[i],
						  hists_dr_lep0_tau[i],
						  hists_dr_lep1_tau[i],
						  hists_mass_lep0_tau[i],
						  hists_mass_lep1_tau[i]);
	}

	// Add signal histograms
	
	// Draw histograms
	drawHistograms<TH1D*>("MVA_2lss_ttV", hists_MVA_2lss_ttV, samples);
	drawHistograms<TH1D*>("MVA_2lss_ttbar", hists_MVA_2lss_ttbar, samples);
	drawHistograms<TH1D*>("MT_met_lep0", hists_MT_met_lep0, samples);
	drawHistograms<TH1D*>("mindr_lep0_jet", hists_lep0_conept, samples);
	drawHistograms<TH1D*>("mindr_lep1_jet", hists_lep1_conept, samples);
	drawHistograms<TH1D*>("lep0_conept", hists_lep0_conept, samples);
	drawHistograms<TH1D*>("lep1_conept", hists_lep1_conept, samples);
	drawHistograms<TH1D*>("avg_dr_jet", hists_avg_dr_jet, samples);
	drawHistograms<TH1D*>("tau_decaymode", hists_tau_decaymode, samples);
	drawHistograms<TH1D*>("dr_lep0_tau", hists_dr_lep0_tau, samples);
	drawHistograms<TH1D*>("dr_lep1_tau", hists_dr_lep1_tau, samples);
	drawHistograms<TH1D*>("mass_lep0_tau", hists_mass_lep0_tau, samples);
	drawHistograms<TH1D*>("mass_lep1_tau", hists_mass_lep1_tau, samples);
	
	// free memories
	for (auto& h : hists_MVA_2lss_ttV) delete h;	
	for (auto& h : hists_MVA_2lss_ttbar) delete h;
	for (auto& h : hists_MT_met_lep0) delete h;
	for (auto& h : hists_mindr_lep0_jet) delete h;
	for (auto& h : hists_mindr_lep1_jet) delete h;
	for (auto& h : hists_lep0_conept) delete h;
	for (auto& h : hists_lep1_conept) delete h;
	for (auto& h : hists_avg_dr_jet) delete h;
	for (auto& h : hists_tau_decaymode) delete h;
	for (auto& h : hists_dr_lep0_tau) delete h;
	for (auto& h : hists_dr_lep1_tau) delete h;
	for (auto& h : hists_mass_lep0_tau) delete h;
	for (auto& h : hists_mass_lep1_tau) delete h;

	for (auto& tree : trees) delete tree;
	for (auto& file : files) delete file;
}


void fillHistofromTree(TTree* tree,
					   TH1D* h_MVA_2lss_ttV,
					   TH1D* h_MVA_2lss_ttbar,
					   TH1D* h_MT_met_lep0,
					   TH1D* h_avg_dr_jet,
					   TH1D* h_mindr_lep0_jet,
					   TH1D* h_mindr_lep1_jet,
					   TH1D* h_lep0_conept,
					   TH1D* h_lep1_conept,
					   TH1D* h_tau_decaymode,
					   TH1D* h_dr_lep0_tau,
					   TH1D* h_dr_lep1_tau,
					   TH1D* h_mass_lep0_tau,
					   TH1D* h_mass_lep1_tau
					   )
{
	int nEntries = tree->GetEntries();

	double evtWeight;
	double MVA_2lss_ttV;
	double MVA_2lss_ttbar;
	double MT_met_lep0;
	double mindr_lep0_jet;
	double mindr_lep1_jet;
	double lep0_conept;
	double lep1_conept;
	double avg_dr_jet;

	int n_mu;
	double mu0_pt;
	double mu0_eta;
	double mu0_phi;
	double mu0_E;
	int mu0_charge;
	double mu0_dxy;
	double mu0_dz;
	double mu1_pt;
	double mu1_eta;
	double mu1_phi;
	double mu1_E;
	int mu1_charge;
	double mu1_dxy;
	double mu1_dz;	
	int n_ele;
	double ele0_pt;
	double ele0_eta;
	double ele0_phi;
	double ele0_E;
	int ele0_charge;
	double ele0_dxy;
	double ele0_dz;
	double ele1_pt;
	double ele1_eta;
	double ele1_phi;
	double ele1_E;
	int ele1_charge;
	double ele1_dxy;
	double ele1_dz;
	
	double tau_pt;
	double tau_eta;
	double tau_phi;
	double tau_E;
	double tau_dxy;
	double tau_dz;
	int tau_charge;
	int tau_decayMode;
	
	tree->SetBranchAddress("evtWeight", &evtWeight);
	tree->SetBranchAddress("MVA_2lss_ttV", &MVA_2lss_ttV);
	tree->SetBranchAddress("MVA_2lss_ttbar", &MVA_2lss_ttbar);
	tree->SetBranchAddress("MT_met_lep0", &MT_met_lep0);
	tree->SetBranchAddress("mindr_lep0_jet", &mindr_lep0_jet);
	tree->SetBranchAddress("mindr_lep1_jet", &mindr_lep1_jet);
	tree->SetBranchAddress("lep0_conept", &lep0_conept);
	tree->SetBranchAddress("lep1_conept", &lep1_conept);
	tree->SetBranchAddress("avg_dr_jet", &avg_dr_jet);

	tree->SetBranchAddress("n_fakeablesel_mu", &n_mu);
	tree->SetBranchAddress("mu0_pt", &mu0_pt);
	tree->SetBranchAddress("mu0_eta", &mu0_eta);
	tree->SetBranchAddress("mu0_phi", &mu0_phi);
	tree->SetBranchAddress("mu0_E", &mu0_E);
	tree->SetBranchAddress("mu0_charge", &mu0_charge);
	tree->SetBranchAddress("mu0_dxy", &mu0_dxy);
	tree->SetBranchAddress("mu0_dz", &mu0_dz);
	tree->SetBranchAddress("mu1_pt", &mu1_pt);
	tree->SetBranchAddress("mu1_eta", &mu1_eta);
	tree->SetBranchAddress("mu1_phi", &mu1_phi);
	tree->SetBranchAddress("mu1_E", &mu1_E);
	tree->SetBranchAddress("mu1_charge", &mu1_charge);
	tree->SetBranchAddress("mu1_dxy", &mu1_dxy);
	tree->SetBranchAddress("mu1_dz", &mu1_dz);
	tree->SetBranchAddress("n_fakeablesel_ele", &n_ele);
	tree->SetBranchAddress("ele0_pt", &ele0_pt);
	tree->SetBranchAddress("ele0_eta", &ele0_eta);
	tree->SetBranchAddress("ele0_phi", &ele0_phi);
	tree->SetBranchAddress("ele0_E", &ele0_E);
	tree->SetBranchAddress("ele0_charge", &ele0_charge);
	tree->SetBranchAddress("ele0_dxy", &ele0_dxy);
	tree->SetBranchAddress("ele0_dz", &ele0_dz);
	tree->SetBranchAddress("ele1_pt", &ele1_pt);
	tree->SetBranchAddress("ele1_eta", &ele1_eta);
	tree->SetBranchAddress("ele1_phi", &ele1_phi);
	tree->SetBranchAddress("ele1_E", &ele1_E);
	tree->SetBranchAddress("ele1_charge", &ele1_charge);
	tree->SetBranchAddress("ele1_dxy", &ele1_dxy);
	tree->SetBranchAddress("ele1_dz", &ele1_dz);

	tree->SetBranchAddress("tau0_pt", &tau_pt);
	tree->SetBranchAddress("tau0_eta", &tau_eta);
	tree->SetBranchAddress("tau0_phi", &tau_phi);
	tree->SetBranchAddress("tau0_E", &tau_E);
	tree->SetBranchAddress("tau0_dxy", &tau_dxy);
	tree->SetBranchAddress("tau0_dz", &tau_dz);
	tree->SetBranchAddress("tau0_charge", &tau_charge);
	tree->SetBranchAddress("tau0_decayMode", &tau_decayMode);

	for (int ievt = 0; ievt < nEntries; ++ievt) {
		tree->GetEntry(ievt);

		h_MVA_2lss_ttV -> Fill(MVA_2lss_ttV);
		h_MVA_2lss_ttbar -> Fill(MVA_2lss_ttbar);
		h_MT_met_lep0 -> Fill(MT_met_lep0);
		h_mindr_lep0_jet -> Fill(mindr_lep0_jet);
		h_mindr_lep1_jet -> Fill(mindr_lep1_jet);
		h_lep0_conept -> Fill(lep0_conept);
		h_lep1_conept -> Fill(lep1_conept);
		h_avg_dr_jet -> Fill(avg_dr_jet);

		assert( n_mu+n_ele==2 );
		
		TLorentzVector lep0;
		TLorentzVector lep1;
		
		if (n_mu == 2) {
			lep0.SetPtEtaPhiE(mu0_pt,mu0_eta,mu0_phi,mu0_E);
			lep1.SetPtEtaPhiE(mu1_pt,mu1_eta,mu1_phi,mu1_E);
		}
		else if (n_ele == 2) {
			lep0.SetPtEtaPhiE(ele0_pt,ele0_eta,ele0_phi,ele0_E);
			lep1.SetPtEtaPhiE(ele1_pt,ele1_eta,ele1_phi,ele1_E);
		}
		else {
			if (mu0_pt >= ele0_pt) {
				lep0.SetPtEtaPhiE(mu0_pt,mu0_eta,mu0_phi,mu0_E);
				lep1.SetPtEtaPhiE(ele0_pt,ele0_eta,ele0_phi,ele0_E);
			}
			else {
				lep0.SetPtEtaPhiE(ele0_pt,ele0_eta,ele0_phi,ele0_E);
				lep1.SetPtEtaPhiE(mu0_pt,mu0_eta,mu0_phi,mu0_E);
			}
		}

		TLorentzVector tau;
		tau.SetPtEtaPhiE(tau_pt,tau_eta,tau_phi,tau_E);

		double dR_lep0_tau = lep0.DeltaR(tau);
		double dR_lep1_tau = lep1.DeltaR(tau);
		double mass_lep0_tau = (lep0+tau).M();
		double mass_lep1_tau = (lep1+tau).M();

		h_tau_decaymode -> Fill(tau_decayMode);
		h_dr_lep0_tau -> Fill(dR_lep0_tau);
		h_dr_lep1_tau -> Fill(dR_lep1_tau);
		h_mass_lep0_tau -> Fill(mass_lep0_tau);
		h_mass_lep1_tau -> Fill(mass_lep1_tau);
		
	}
}
