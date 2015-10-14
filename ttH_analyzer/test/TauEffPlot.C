#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TPaveText.h"

using namespace std;


void TauEffPlot(const TString signal =
				"/uscms/home/ztao/work/CU_ttH_WD/Outputs/CU_ttH_EDA_output_sig.root"
				)
{
	TFile* f_sig = new TFile(signal);

	f_sig->cd("ttHsyncExercise");
	
	// Get historamgs
	TH1D* h_num_genHadTau = (TH1D*)f_sig->Get("ttHsyncExercise/h_num_genHadTau");
	TH1D* h_genHadTau_pt = (TH1D*)f_sig->Get("ttHsyncExercise/h_genHadTau_pt");
	TH1D* h_genHadTau_eta = (TH1D*)f_sig->Get("ttHsyncExercise/h_genHadTau_eta");
	TH1D* h_genHadTau_phi = (TH1D*)f_sig->Get("ttHsyncExercise/h_genHadTau_phi");

	TH1D* h_nTau_noniso = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_num_selectedTau_noniso");
	TH1D* h_eff_Tau_noniso_pt = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_selectedTau_noniso_genpt");
	TH1D* h_eff_Tau_noniso_eta = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_selectedTau_noniso_geneta");
	TH1D* h_eff_Tau_noniso_phi = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_selectedTau_noniso_genphi");

	TH1D* h_nTau_loose = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_num_selectedTau_loose");
	TH1D* h_eff_Tau_loose_pt = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_selectedTau_loose_genpt");
	TH1D* h_eff_Tau_loose_eta = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_selectedTau_loose_geneta");
	TH1D* h_eff_Tau_loose_phi = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_selectedTau_loose_genphi");

	TH1D* h_nTau_medium = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_num_selectedTau_medium");
	TH1D* h_eff_Tau_medium_pt = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_selectedTau_medium_genpt");
	TH1D* h_eff_Tau_medium_eta = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_selectedTau_medium_geneta");
	TH1D* h_eff_Tau_medium_phi = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_selectedTau_medium_genphi");

	TH1D* h_nTau_tight = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_num_selectedTau_tight");
	TH1D* h_eff_Tau_tight_pt = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_selectedTau_tight_genpt");
	TH1D* h_eff_Tau_tight_eta = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_selectedTau_tight_geneta");
	TH1D* h_eff_Tau_tight_phi = (TH1D*)f_sig->
		Get("ttHsyncExercise/h_selectedTau_tight_genphi");

	h_num_genHadTau -> Sumw2();
	h_genHadTau_pt  -> Sumw2();
	h_genHadTau_eta -> Sumw2();
	h_genHadTau_phi -> Sumw2();
	h_nTau_noniso        -> Sumw2();
	h_eff_Tau_noniso_pt  -> Sumw2();
	h_eff_Tau_noniso_eta -> Sumw2();
	h_eff_Tau_noniso_phi -> Sumw2();
	h_nTau_loose        -> Sumw2();
	h_eff_Tau_loose_pt  -> Sumw2();
	h_eff_Tau_loose_eta -> Sumw2();
	h_eff_Tau_loose_phi -> Sumw2();
	h_nTau_medium        -> Sumw2();
	h_eff_Tau_medium_pt  -> Sumw2();
	h_eff_Tau_medium_eta -> Sumw2();
	h_eff_Tau_medium_phi -> Sumw2();
	h_nTau_tight        -> Sumw2();
	h_eff_Tau_tight_pt  -> Sumw2();
	h_eff_Tau_tight_eta -> Sumw2();
	h_eff_Tau_tight_phi -> Sumw2();
	
	// Rebin
	h_genHadTau_pt -> Rebin(5);
	h_eff_Tau_loose_pt  -> Rebin(5);
	h_eff_Tau_medium_pt -> Rebin(5);
	h_eff_Tau_tight_pt  -> Rebin(5);
	h_eff_Tau_noniso_pt -> Rebin(5);

	h_genHadTau_eta -> Rebin(5);
	h_eff_Tau_loose_eta  -> Rebin(5);
	h_eff_Tau_medium_eta -> Rebin(5);
	h_eff_Tau_tight_eta  -> Rebin(5);
	h_eff_Tau_noniso_eta -> Rebin(5);

	h_genHadTau_phi -> Rebin(4);
	h_eff_Tau_loose_phi  -> Rebin(4);
	h_eff_Tau_medium_phi -> Rebin(4);
	h_eff_Tau_tight_phi  -> Rebin(4);
	h_eff_Tau_noniso_phi -> Rebin(4);

	// Efficiency
	h_eff_Tau_noniso_pt  -> Divide(h_genHadTau_pt);
	h_eff_Tau_noniso_eta -> Divide(h_genHadTau_eta);
	h_eff_Tau_noniso_phi -> Divide(h_genHadTau_phi);
	
	h_eff_Tau_loose_pt  -> Divide(h_genHadTau_pt);
	h_eff_Tau_loose_eta -> Divide(h_genHadTau_eta);
	h_eff_Tau_loose_phi -> Divide(h_genHadTau_phi);

	h_eff_Tau_medium_pt  -> Divide(h_genHadTau_pt);
	h_eff_Tau_medium_eta -> Divide(h_genHadTau_eta);
	h_eff_Tau_medium_phi -> Divide(h_genHadTau_phi);

	h_eff_Tau_tight_pt  -> Divide(h_genHadTau_pt);
	h_eff_Tau_tight_eta -> Divide(h_genHadTau_eta);
	h_eff_Tau_tight_phi -> Divide(h_genHadTau_phi);

	TCanvas c;
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(1);
	
	h_eff_Tau_loose_pt -> SetTitle("ttH, H->#tau#tau (M125_13TeV_Spring15)");
	h_eff_Tau_loose_pt -> GetXaxis() -> SetTitle("gen #tau_{h} p_{T} [GeV]");
	h_eff_Tau_loose_pt -> GetYaxis() -> SetTitle("Selection Efficiency");
	h_eff_Tau_loose_pt -> SetMinimum(0);
	h_eff_Tau_loose_pt -> SetLineColor(3);
	h_eff_Tau_loose_pt -> Draw("E");

	h_eff_Tau_medium_pt -> SetLineColor(2);
	h_eff_Tau_tight_pt  -> SetLineColor(1);
	//h_eff_Tau_noniso_pt  -> SetLineColor(4);

	h_eff_Tau_medium_pt ->Draw("E same");
	h_eff_Tau_tight_pt  ->Draw("E same");
	//h_eff_Tau_noniso_pt  ->Draw("E same");

	TLegend* leg_pt = new TLegend(0.6,0.68,0.88,0.85);
	leg_pt->AddEntry((TObject*)0, "p_{T}>20GeV, |#eta|<2.5","");
	//leg_pt->AddEntry(h_eff_Tau_noniso_pt,"nonIso","l");
	leg_pt->AddEntry(h_eff_Tau_loose_pt,"loose","l");
	leg_pt->AddEntry(h_eff_Tau_medium_pt,"medium","l");
	leg_pt->AddEntry(h_eff_Tau_tight_pt,"tight","l");
	leg_pt->Draw("same");

	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/tauEffpt.pdf");

	h_eff_Tau_loose_eta -> SetTitle("ttH, H->#tau#tau (M125_13TeV_Spring15)");
	h_eff_Tau_loose_eta -> GetXaxis() -> SetTitle("gen #tau_{h} #eta");
	h_eff_Tau_loose_eta -> GetYaxis() -> SetTitle("Selection Efficiency");
	h_eff_Tau_loose_eta -> SetMinimum(0);
	h_eff_Tau_loose_eta -> SetLineColor(3);
	h_eff_Tau_loose_eta -> Draw("E");

	h_eff_Tau_medium_eta -> SetLineColor(2);
	h_eff_Tau_tight_eta  -> SetLineColor(1);
	//h_eff_Tau_noniso_eta  -> SetLineColor(4);

	h_eff_Tau_medium_eta ->Draw("E same");
	h_eff_Tau_tight_eta  ->Draw("E same");
	//h_eff_Tau_noniso_eta  ->Draw("E same");

	TLegend* leg_eta = new TLegend(0.6,0.68,0.88,0.85);
	leg_eta->AddEntry((TObject*)0, "p_{T}>20GeV, |#eta|<2.5","");
	//leg_eta->AddEntry(h_eff_Tau_noniso_eta,"nonIso","l");
	leg_eta->AddEntry(h_eff_Tau_loose_eta,"loose","l");
	leg_eta->AddEntry(h_eff_Tau_medium_eta,"medium","l");
	leg_eta->AddEntry(h_eff_Tau_tight_eta,"tight","l");
	leg_eta->Draw("same");

	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/tauEffeta.pdf");

	h_eff_Tau_loose_phi -> SetTitle("ttH, H->#tau#tau (M125_13TeV_Spring15)");
	h_eff_Tau_loose_phi -> GetXaxis() -> SetTitle("gen #tau_{h} #phi");
	h_eff_Tau_loose_phi -> GetYaxis() -> SetTitle("Selection Efficiency");
	h_eff_Tau_loose_phi -> SetMinimum(0);
	h_eff_Tau_loose_phi -> SetLineColor(3);
	h_eff_Tau_loose_phi -> Draw("E");

	h_eff_Tau_medium_phi -> SetLineColor(2);
	h_eff_Tau_tight_phi  -> SetLineColor(1);
	//h_eff_Tau_loose_phi  -> SetLineColor(4);
	
	h_eff_Tau_medium_phi ->Draw("E same");
	h_eff_Tau_tight_phi  ->Draw("E same");
	h_eff_Tau_noniso_phi  ->Draw("E same");
	
	TLegend* leg_phi = new TLegend(0.6,0.68,0.88,0.85);
	leg_phi->AddEntry((TObject*)0, "p_{T}>20GeV, |#eta|<2.5","");
	//leg_phi->AddEntry(h_eff_Tau_noniso_phi,"nonIso","l");
	leg_phi->AddEntry(h_eff_Tau_loose_phi,"loose","l");
	leg_phi->AddEntry(h_eff_Tau_medium_phi,"medium","l");
	leg_phi->AddEntry(h_eff_Tau_tight_phi,"tight","l");
	leg_phi->Draw("same");

	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/tauEffphi.pdf");


	int nevt = h_num_genHadTau -> GetEntries();
	h_num_genHadTau -> Scale(1.0/nevt);
	h_nTau_noniso -> Scale(1.0/nevt);
	h_nTau_loose  -> Scale(1.0/nevt);
	h_nTau_medium -> Scale(1.0/nevt);
	h_nTau_tight  -> Scale(1.0/nevt);
	
	TCanvas c2;
	gStyle->SetOptStat(10);

	h_nTau_tight->SetTitle("ttH, H->#tau#tau (M125_13TeV_Spring15)");
	h_nTau_tight->GetXaxis()->SetTitle("n_{#tau_{h}}");
	h_nTau_tight->GetYaxis()->SetTitle("Fraction of events");
	h_nTau_tight->SetLineColor(1);
	h_nTau_tight->Draw("E");

	h_nTau_noniso -> SetLineColor(5);
	h_nTau_loose  -> SetLineColor(3);
	h_nTau_medium -> SetLineColor(2);
	h_num_genHadTau -> SetLineColor(4);

	h_nTau_noniso  -> Draw("E same");
	h_nTau_loose  -> Draw("E same");
	h_nTau_medium -> Draw("E same");
	h_num_genHadTau -> Draw("E same");

	TLegend* legn = new TLegend(0.6,0.68,0.88,0.85);
	legn->AddEntry((TObject*)0, "p_{T}>20 GeV  |#eta|<2.5", "");
	legn->AddEntry(h_nTau_noniso, "noniso","l");
	legn->AddEntry(h_nTau_loose,"loose","l");
	legn->AddEntry(h_nTau_medium, "medium","l");
	legn->AddEntry(h_nTau_tight, "tight","l");
	legn->AddEntry(h_num_genHadTau,"gen","l");
	legn->Draw("same");

	c2.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/tauEffnum.pdf");
	
	
	delete f_sig;
}
