#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

void MVATreeMaker(TString label)
{
	// Read ntuple
    TString ntuple;
	if (label == "signal") 
	    ntuple = "/eos/uscms/store/user/ztao/ttHToTT_M125_13TeV_powheg_pythia8/ttHToTauTau_Ntuple_signal/151209_202153/0000/CU_ttH_EDA_output.root";
	else if (label == "TTJets")
		ntuple = "/eos/uscms/store/user/ztao/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ttHToTauTau_Ntuple_TTJets/151209_202324/0000/CU_ttH_EDA_output.root";
	
	TFile* input = new TFile(ntuple);
	TTree* tree = (TTree*) input -> Get("ttHtautau/EventTree");

	const int nEntries = tree->GetEntries();
	
	if (nEntries == 0) {
		cout << "File doesn't exist or is empty, returning..." << endl;
		cout << endl;
		return;
	}
	else {
		cout << "Number of entries :" << nEntries << endl;
	}

	// Define source tree leafs and branches
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

	// Create a target tree and output file
	TFile *output = new TFile("./mvaTree_"+label+".root", "RECREATE");
	TTree* tree_mva = new TTree("mvaTree", "mvaTree");
	// Define target tree leafs and branches
	float dRlep1tau;
	float dRlep2tau;
	float lep1dxy;
	float lep1dz;
	float lep2dxy;
	float lep2dz;
	float lepmaxdxy;
	float visM_taulep1;
	float visM_taulep2;
	float visM_taulepmaxdxy;
	float visM_taulepmindR;
	float MET;
	float METSig;
	float MHT;
	float HT;
	float metLD;
	float visMtot;
	float dilepmass;
	float dRlep1jet1;
	float dRlep1jet2;
	float dRlep2jet1;
	float dRlep2jet2;
	// Setup branches
	tree_mva -> Branch("dRlep1tau", &dRlep1tau);
	tree_mva -> Branch("dRlep2tau", &dRlep2tau);
	tree_mva -> Branch("lep1dxy", &lep1dxy);
	tree_mva -> Branch("lep1dz", &lep1dz);
	tree_mva -> Branch("lep2dxy", &lep2dxy);
	tree_mva -> Branch("lep2dz", &lep2dz);
	tree_mva -> Branch("lepmaxdxy", &lepmaxdxy);
	tree_mva -> Branch("visM_taulep1", &visM_taulep1);
	tree_mva -> Branch("visM_taulep2", &visM_taulep2);
	tree_mva -> Branch("visM_taulepmaxdxy", &visM_taulepmaxdxy);
	tree_mva -> Branch("visM_taulepmindR", &visM_taulepmindR);
	tree_mva -> Branch("MET", &MET);
	tree_mva -> Branch("METSig", &METSig);
	tree_mva -> Branch("MHT", &MHT);
	tree_mva -> Branch("HT", &HT);
	tree_mva -> Branch("metLD", &metLD);
	tree_mva -> Branch("visMtot", &visMtot);
	tree_mva -> Branch("dilepmass", &dilepmass);
	tree_mva -> Branch("dRlep1jet1", &dRlep1jet1);
	tree_mva -> Branch("dRlep1jet2", &dRlep1jet2);
	tree_mva -> Branch("dRlep2jet1", &dRlep2jet1);
	tree_mva -> Branch("dRlep2jet2", &dRlep2jet2);
	
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
		MET = TMath::Sqrt(met_x*met_x+met_y*met_y);
		// MET Significance
		METSig = metsig;
		
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
		HT = 0;
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
		MHT = TMath::Sqrt(MHT_x*MHT_x+MHT_y*MHT_y);

		// metLD
		metLD = 0.00397*MET+ 0.00265*MHT;
		
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
		
		visMtot = p4tot.M();
		
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
	
		float deta_l1T = lep_eta->at(ilep1)-loose_tau_eta->at(0);
		float dphi_l1T = TVector2::Phi_mpi_pi(lep_phi->at(ilep1)-loose_tau_phi->at(0));
		dRlep1tau = TMath::Sqrt(deta_l1T*deta_l1T+dphi_l1T*dphi_l1T);
		
		float deta_l2T = lep_eta->at(ilep2)-loose_tau_eta->at(0);
		float dphi_l2T = TVector2::Phi_mpi_pi(lep_phi->at(ilep2)-loose_tau_phi->at(0));
		dRlep2tau = TMath::Sqrt(deta_l2T*deta_l2T+dphi_l2T*dphi_l2T);

		lep1dxy = lep_vtx_dxy->at(ilep1);
		lep1dz = lep_vtx_dz->at(ilep1);
		lep2dxy = lep_vtx_dxy->at(ilep2);
		lep2dz = lep_vtx_dz->at(ilep2);
		lepmaxdxy = std::max(lep_vtx_dxy->at(ilep1),lep_vtx_dxy->at(ilep2));

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
		visM_taulep1 = (tau_p4+lep1_p4).M();
		visM_taulep2 = (tau_p4+lep2_p4).M();
		visM_taulepmaxdxy = (tau_p4+lep_maxdxy_p4).M();
		visM_taulepmindR = (tau_p4+lep_mindR_p4).M();

		// Di-lepton mass
		dilepmass = (lep1_p4+lep2_p4).M();

		dRlep1jet1 = lep1_p4.DeltaR(jets[0]);
		dRlep1jet2 = lep1_p4.DeltaR(jets[1]);
		dRlep2jet1 = lep2_p4.DeltaR(jets[0]);
		dRlep2jet2 = lep2_p4.DeltaR(jets[1]);

		// Fill target tree		
		tree_mva -> Fill();
		
	} // end of event loop
	
	tree_mva -> Write();
}
