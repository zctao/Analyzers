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
#include "../macro/Cross_Sections.h"

using namespace std;

int fillHistofromTree(TTree*,
					   TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,
					   TH1D*,TH1D*,TH1D*,TH1D*,TH1D*);
vector<TH1D*> combineHistofromTrees(vector<TTree*>, vector<int>, vector<TString>);
TString getCommonPrefix(TString, TString, char);
void scale_histograms(TH1D*,int,float,float);

float LUMI = 12.9 * 1000;  // 1/pb

void treeAnalyzer
(
 //vector<TString> samples = {"ttH_htt","ttH_hww","ttH_hzz", "TTW", "TTZ", "TTJets"},
 //vector<TString> samples = {"ttH", "TTW", "TTZ", "TTJets"},
 vector<vector<TString>> samples = {{"ttH"}, {"TTW"}, {"TTZ"}, {"TTJets_ll", "TTJets_lt","TTJets_ltbar"}},
 TString directory = "/Users/ztao/Documents/ttH/Outputs/80X/"
 )
{

	int nchannels = samples.size();
	vector<TString> channels;

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

	// reserve space
	hists_MVA_2lss_ttV.resize(nchannels);
	hists_MVA_2lss_ttbar.resize(nchannels);
	hists_MT_met_lep0.resize(nchannels);
	hists_mindr_lep0_jet.resize(nchannels);
	hists_mindr_lep1_jet.resize(nchannels);
	hists_lep0_conept.resize(nchannels);
	hists_lep1_conept.resize(nchannels);
	hists_avg_dr_jet.resize(nchannels);
	hists_tau_decaymode.resize(nchannels);
	hists_dr_lep0_tau.resize(nchannels);
	hists_dr_lep1_tau.resize(nchannels);
	hists_mass_lep0_tau.resize(nchannels);
	hists_mass_lep1_tau.resize(nchannels);
	
	// Fill histograms
	
	for (int i = 0; i < nchannels; i++) {

		vector<TString> snames = samples[i];

		// Get channel name
		assert(snames.size()>0);
		if (snames.size()>1) {
			TString prefix = getCommonPrefix(snames[0],snames[1], '_');
			channels.push_back(prefix);
		}
		else
			channels.push_back(snames[0]);

		cout << "channel : " << channels[i] << endl;
		
		vector<TTree*> trees;
		vector<int> nProcessed;
		
		for (auto& sn : snames) {
			// Read file
			TFile* f = new TFile(directory+"output_"+ sn +".root");
			TTree* tree = (TTree*) f->Get("ttHtaus/eventTree");
			TH1I* h_nProcessed = (TH1I*)f->Get("ttHtaus/h_nProcessed");
			int nentries = h_nProcessed->GetEntries();

			trees.push_back(tree);
			nProcessed.push_back(nentries);
		}

		vector<TH1D*> histo_set = combineHistofromTrees(trees, nProcessed, snames);
		
		// delete trees
		for (auto& t : trees)
			delete t;
		
		float yields = histo_set[0]->Integral();
		cout << "yields : " << yields << endl;

		hists_MVA_2lss_ttV[i] = histo_set[0];
		hists_MVA_2lss_ttbar[i] = histo_set[1];
		hists_MT_met_lep0[i] = histo_set[2];
		hists_avg_dr_jet[i] = histo_set[3];
		hists_mindr_lep0_jet[i] = histo_set[4];
		hists_mindr_lep1_jet[i] = histo_set[5];
		hists_lep0_conept[i] = histo_set[6];
		hists_lep1_conept[i] = histo_set[7];
		hists_tau_decaymode[i] = histo_set[8];
		hists_dr_lep0_tau[i] = histo_set[9];
		hists_dr_lep1_tau[i] = histo_set[10];
		hists_mass_lep0_tau[i] = histo_set[11];
		hists_mass_lep1_tau[i] = histo_set[12];
			
	}
	
	// Draw histograms
	drawHistograms<TH1D*>("MVA_2lss_ttV", hists_MVA_2lss_ttV, channels);
	drawHistograms<TH1D*>("MVA_2lss_ttbar", hists_MVA_2lss_ttbar, channels);
	drawHistograms<TH1D*>("MT_met_lep0", hists_MT_met_lep0, channels);
	drawHistograms<TH1D*>("mindr_lep0_jet", hists_lep0_conept, channels);
	drawHistograms<TH1D*>("mindr_lep1_jet", hists_lep1_conept, channels);
	drawHistograms<TH1D*>("lep0_conept", hists_lep0_conept, channels);
	drawHistograms<TH1D*>("lep1_conept", hists_lep1_conept, channels);
	drawHistograms<TH1D*>("avg_dr_jet", hists_avg_dr_jet, channels);
	drawHistograms<TH1D*>("tau_decaymode", hists_tau_decaymode, channels);
	drawHistograms<TH1D*>("dr_lep0_tau", hists_dr_lep0_tau, channels);
	drawHistograms<TH1D*>("dr_lep1_tau", hists_dr_lep1_tau, channels);
	drawHistograms<TH1D*>("mass_lep0_tau", hists_mass_lep0_tau, channels);
	drawHistograms<TH1D*>("mass_lep1_tau", hists_mass_lep1_tau, channels);
	
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

}


int fillHistofromTree(TTree* tree,
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
	
	return nEntries;
}

TString getCommonPrefix(TString s1, TString s2, char stopper)
{
	int nsize = min(s1.Sizeof(),s2.Sizeof());

	TString os = "";

	for (int i = 0; i < nsize; i++) {
		if (s1[i]==stopper)
			return os;	
		if (s1[i]==s2[i])
			os.Append(s1[i]);
	}

	return os;
}

vector<TH1D*> combineHistofromTrees(vector<TTree*> trees,
							vector<int> nProcessed,
							vector<TString> samples
						   )
{
	assert(trees.size()==nProcessed.size());
	assert(trees.size()==samples.size());

	// number of samples to combine
	int nsample = trees.size();

	vector<TH1D*> vec_h_MVA_2lss_ttV;
	vector<TH1D*> vec_h_MVA_2lss_ttbar;
	vector<TH1D*> vec_h_MT_met_lep0;
	vector<TH1D*> vec_h_avg_dr_jet;
	vector<TH1D*> vec_h_mindr_lep0_jet;
	vector<TH1D*> vec_h_mindr_lep1_jet;
	vector<TH1D*> vec_h_lep0_conept;
	vector<TH1D*> vec_h_lep1_conept;
	vector<TH1D*> vec_h_tau_decaymode;
	vector<TH1D*> vec_h_dr_lep0_tau;
	vector<TH1D*> vec_h_dr_lep1_tau;
	vector<TH1D*> vec_h_mass_lep0_tau;
	vector<TH1D*> vec_h_mass_lep1_tau;

    vec_h_MVA_2lss_ttV.resize(nsample+1);
	vec_h_MVA_2lss_ttbar.resize(nsample+1);
    vec_h_MT_met_lep0.resize(nsample+1);
	vec_h_avg_dr_jet.resize(nsample+1);
	vec_h_mindr_lep0_jet.resize(nsample+1);
    vec_h_mindr_lep1_jet.resize(nsample+1);
	vec_h_lep0_conept.resize(nsample+1);
    vec_h_lep1_conept.resize(nsample+1);
	vec_h_tau_decaymode.resize(nsample+1);
	vec_h_dr_lep0_tau.resize(nsample+1);
	vec_h_dr_lep1_tau.resize(nsample+1);
	vec_h_mass_lep0_tau.resize(nsample+1);
	vec_h_mass_lep1_tau.resize(nsample+1);

	for (int i = 0; i < nsample+1; i++) {
		
		// setup histograms
		vec_h_MVA_2lss_ttV[i] =
			new TH1D(Form("h_MVA_2lss_ttV_%d",i), "", 10, -1.0, 1.0);
	    vec_h_MVA_2lss_ttbar[i] =
			new TH1D(Form("h_MVA_2lss_ttbar_%d",i), "", 10, -1.0, 1.0);
		vec_h_MT_met_lep0[i] =
			new TH1D(Form("h_MT_met_lep0_%d",i), "", 10, 0, 400);
		vec_h_avg_dr_jet[i] =
			new TH1D(Form("h_avg_dr_jet_%d",i), "", 10, 0, 4);
		vec_h_mindr_lep0_jet[i] =
			new TH1D(Form("h_mindr_lep0_jet_%d",i), "", 10, 0, 4);
		vec_h_mindr_lep1_jet[i] =
			new TH1D(Form("h_mindr_lep1_jet_%d",i), "", 10, 0, 4);
		vec_h_lep0_conept[i] =
			new TH1D(Form("h_lep0_conept_%d",i), "", 10, 0, 250);
		vec_h_lep1_conept[i] =
			new TH1D(Form("h_lep1_conept_%d",i), "", 10, 0, 150);
		vec_h_tau_decaymode[i] =
			new TH1D(Form("h_tau_decaymode_%d",i), "", 18, 0, 18);
		vec_h_dr_lep0_tau[i] =
			new TH1D(Form("h_dr_lep0_tau_%d",i), "", 10, 0., 4.);
		vec_h_dr_lep1_tau[i] =
			new TH1D(Form("h_dr_lep1_tau_%d",i), "", 10, 0., 4.);
		vec_h_mass_lep0_tau[i] =
			new TH1D(Form("h_mass_lep0_tau_%d",i), "", 10, 0., 400.);
		vec_h_mass_lep1_tau[i] =
			new TH1D(Form("h_mass_lep1_tau_%d",i), "", 10, 0., 400.);
	}


	cout << "sample tree nEntries :";
	
	for (int i = 0; i< nsample; i++){
		
		// fill histograms
		int nentries = fillHistofromTree(trees[i],
						  vec_h_MVA_2lss_ttV[i],
						  vec_h_MVA_2lss_ttbar[i],
						  vec_h_MT_met_lep0[i],
						  vec_h_avg_dr_jet[i],
						  vec_h_mindr_lep0_jet[i],
						  vec_h_mindr_lep1_jet[i],
						  vec_h_lep0_conept[i],
						  vec_h_lep1_conept[i],
						  vec_h_tau_decaymode[i],
						  vec_h_dr_lep0_tau[i],
						  vec_h_dr_lep1_tau[i],
						  vec_h_mass_lep0_tau[i],
						  vec_h_mass_lep1_tau[i]
						  );
		
		cout << " " << nentries;
		
		// scale histograms
		float xs = xsection::xsection[string(samples[i])];

		scale_histograms(vec_h_MVA_2lss_ttV[i],nProcessed[i], xs, LUMI);
		scale_histograms(vec_h_MVA_2lss_ttbar[i],nProcessed[i], xs, LUMI);
		scale_histograms(vec_h_MT_met_lep0[i],nProcessed[i], xs, LUMI);
		scale_histograms(vec_h_avg_dr_jet[i],nProcessed[i], xs, LUMI);
		scale_histograms(vec_h_mindr_lep0_jet[i],nProcessed[i], xs, LUMI);
		scale_histograms(vec_h_mindr_lep1_jet[i],nProcessed[i], xs, LUMI);
		scale_histograms(vec_h_lep0_conept[i],nProcessed[i], xs, LUMI);
		scale_histograms(vec_h_lep1_conept[i],nProcessed[i], xs, LUMI);
		scale_histograms(vec_h_tau_decaymode[i],nProcessed[i], xs, LUMI);
		scale_histograms(vec_h_dr_lep0_tau[i],nProcessed[i], xs, LUMI);
		scale_histograms(vec_h_dr_lep1_tau[i],nProcessed[i], xs, LUMI);
		scale_histograms(vec_h_mass_lep0_tau[i],nProcessed[i], xs, LUMI);
		scale_histograms(vec_h_mass_lep1_tau[i],nProcessed[i], xs, LUMI);
		
	}

	cout << endl;
	
	// add up histograms
	for (int i = 0; i < nsample; i++) {
		vec_h_MVA_2lss_ttV[nsample]->Add(vec_h_MVA_2lss_ttV[i]);
		vec_h_MVA_2lss_ttbar[nsample]->Add(vec_h_MVA_2lss_ttbar[i]);
		vec_h_MT_met_lep0[nsample]->Add(vec_h_MT_met_lep0[i]);
		vec_h_avg_dr_jet[nsample]->Add(vec_h_avg_dr_jet[i]);
		vec_h_mindr_lep0_jet[nsample]->Add(vec_h_mindr_lep0_jet[i]);
		vec_h_mindr_lep1_jet[nsample]->Add(vec_h_mindr_lep1_jet[i]);
		vec_h_lep0_conept[nsample]->Add(vec_h_lep0_conept[i]);
		vec_h_lep1_conept[nsample]->Add(vec_h_lep1_conept[i]);
		vec_h_tau_decaymode[nsample]->Add(vec_h_tau_decaymode[i]);
		vec_h_dr_lep0_tau[nsample]->Add(vec_h_dr_lep0_tau[i]);
		vec_h_dr_lep1_tau[nsample]->Add(vec_h_dr_lep1_tau[i]);
		vec_h_mass_lep0_tau[nsample]->Add(vec_h_mass_lep0_tau[i]);
		vec_h_mass_lep1_tau[nsample]->Add(vec_h_mass_lep1_tau[i]);
	}
	
	// free memory
	for (int i = 0; i < nsample; i++) {
		delete vec_h_MVA_2lss_ttV[i];
		delete vec_h_MVA_2lss_ttbar[i];
		delete vec_h_MT_met_lep0[i];
		delete vec_h_avg_dr_jet[i];
		delete vec_h_mindr_lep0_jet[i];
		delete vec_h_mindr_lep1_jet[i];
		delete vec_h_lep0_conept[i];
		delete vec_h_lep1_conept[i];
		delete vec_h_tau_decaymode[i];
		delete vec_h_dr_lep0_tau[i];
		delete vec_h_dr_lep1_tau[i];
		delete vec_h_mass_lep0_tau[i];
		delete vec_h_mass_lep1_tau[i];
	}
	
	vector<TH1D*> histograms;
	
	histograms.push_back(vec_h_MVA_2lss_ttV[nsample]);
	histograms.push_back(vec_h_MVA_2lss_ttbar[nsample]);
	histograms.push_back(vec_h_MT_met_lep0[nsample]);
	histograms.push_back(vec_h_avg_dr_jet[nsample]);
	histograms.push_back(vec_h_mindr_lep0_jet[nsample]);
	histograms.push_back(vec_h_mindr_lep1_jet[nsample]);
	histograms.push_back(vec_h_lep0_conept[nsample]);
	histograms.push_back(vec_h_lep1_conept[nsample]);
	histograms.push_back(vec_h_tau_decaymode[nsample]);
	histograms.push_back(vec_h_dr_lep0_tau[nsample]);
	histograms.push_back(vec_h_dr_lep1_tau[nsample]);
	histograms.push_back(vec_h_mass_lep0_tau[nsample]);
	histograms.push_back(vec_h_mass_lep1_tau[nsample]);
	
	return histograms;
}

void scale_histograms(TH1D* h,
					 int nProcessed, float xsection, float lumi)
{
	h -> Sumw2();
	h -> Scale(lumi * xsection / nProcessed);
}
