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


void Draw_mtopmtautau(const TString input =
					  "/uscms/home/ztao/work/CU_ttH_WD/Outputs/gen_histograms.root" )
{
	// Open root file and get histograms
	TFile* f = new TFile(input);
	
	TH2F* h_mtautau_mtop1 = (TH2F*)f->Get("mTTmtop1");
	TH2F* h_mtautau_mtop2 = (TH2F*)f->Get("mTTmtop2");
	TH2F* h_ptautau_mtop1 = (TH2F*)f->Get("pTTmtop1");
	TH2F* h_ptautau_mtop2 = (TH2F*)f->Get("pTTmtop2");
	TH2F* h_mtop1TT_mtop2TT = (TH2F*)f->Get("mtop1TT_mtop2TT");
	TH2F* h_mplus_mminus = (TH2F*)f->Get("mplus_mminus");
	TH2F* h_mminus_mplus = (TH2F*)f->Get("mminus_mplus");
	TH2F* h_mttTT_mplus = (TH2F*)f->Get("mttTT_mplus");
	TH2F* h_mttTT_mminus = (TH2F*)f->Get("mttTT_mminus");

	TH1F* h_top_higgs_dRmin = (TH1F*)f->Get("dRmin");
	TH1F* h_top_higgs_dRmax = (TH1F*)f->Get("dRmax");

	TCanvas c;
	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(1);
	
	h_mtautau_mtop1->GetXaxis()->SetTitle("m_{#tau#tau} [GeV]");
	h_mtautau_mtop1->GetYaxis()->SetTitle("m_{top1} [GeV]");
	h_mtautau_mtop1->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mtautau_mtop1.pdf");
	
	h_mtautau_mtop2->GetXaxis()->SetTitle("m_{#tau#tau} [GeV]");
	h_mtautau_mtop2->GetYaxis()->SetTitle("m_{top2} [GeV]");
	h_mtautau_mtop2->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mtautau_mtop2.pdf");

	h_ptautau_mtop1->GetXaxis()->SetTitle("|p_{#tau#tau}| [GeV]");
	h_ptautau_mtop1->GetYaxis()->SetTitle("m_{top1} [GeV]");
	h_ptautau_mtop1->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/ptautau_mtop1.pdf");
	
	h_ptautau_mtop2->GetXaxis()->SetTitle("|p_{#tau#tau}| [GeV]");
	h_ptautau_mtop2->GetYaxis()->SetTitle("m_{top2} [GeV]");
	h_ptautau_mtop2->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/ptautau_mtop2.pdf");

	h_mtop1TT_mtop2TT->GetXaxis()->SetTitle("m_{top1+#tau+#tau} [GeV]");
	h_mtop1TT_mtop2TT->GetYaxis()->SetTitle("m_{top2+#tau+#tau} [GeV]");
	h_mtop1TT_mtop2TT->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mtop1TT_mtop2TT.pdf");

	h_mplus_mminus->GetXaxis()->SetTitle("max(m_{top1+#tau+#tau},m_{top2+#tau+#tau}) [GeV]");
	h_mplus_mminus->GetYaxis()->SetTitle("min(m_{top1+#tau+#tau},m_{top2+#tau+#tau}) [GeV]");
	h_mplus_mminus->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mplus_mminus.pdf");

	h_mminus_mplus->GetYaxis()->SetTitle("max(m_{top1+#tau+#tau},m_{top2+#tau+#tau}) [GeV]");
	h_mminus_mplus->GetXaxis()->SetTitle("min(m_{top1+#tau+#tau},m_{top2+#tau+#tau}) [GeV]");
	h_mminus_mplus->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mminus_mplus.pdf");

	h_mttTT_mplus->GetXaxis()->SetTitle("m_{top+top+#tau+#tau}");
	h_mttTT_mplus->GetYaxis()->SetTitle("max(m_{top1+#tau+#tau},m_{top2+#tau+#tau})");
	h_mttTT_mplus->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mttTT_mplus.pdf");

	h_mttTT_mminus->GetXaxis()->SetTitle("m_{top+top+#tau+#tau}");
	h_mttTT_mminus->GetYaxis()->SetTitle("min(m_{top1+#tau+#tau},m_{top2+#tau+#tau})");
	h_mttTT_mminus->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mttTT_mminus.pdf");
	
	gStyle->SetOptStat(1);
	h_top_higgs_dRmax->GetXaxis()->SetTitle("#DeltaR(top,Higgs)");
	h_top_higgs_dRmax->SetTitle("ttH, H->#tau#tau [GEN]");
	h_top_higgs_dRmax->SetLineColor(1);
	h_top_higgs_dRmax->Draw();
	h_top_higgs_dRmax->SetLineColor(2);
	h_top_higgs_dRmin->Draw("same");
	
	TLegend* leg = new TLegend(0.6,0.58,0.88,0.75);
	leg->AddEntry(h_top_higgs_dRmin,"min #DeltaR","l");
	leg->AddEntry(h_top_higgs_dRmax,"max #DeltaR","l");
	leg->Draw("same");
	
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dR_top_Higgs.pdf");
	
	delete f;
}

void MakeTauROCPlot(const TString input =
					"/uscms/home/ztao/work/CU_ttH_WD/Outputs/taueff.root")
{
	// Open root file and get histograms
	TFile* f = new TFile(input);

	TH1F* h_taueff_noniso = (TH1F*)f->Get("h_taueff_noniso");
	TH1F* h_taueff_loose = (TH1F*)f->Get("h_taueff_loose");
	TH1F* h_taueff_medium = (TH1F*)f->Get("h_taueff_medium");
	TH1F* h_taueff_tight = (TH1F*)f->Get("h_taueff_tight");

	TH1F* h_bkgtaupt_noniso = (TH1F*)f->Get("h_bkgtaupt_noniso");
	TH1F* h_bkgtaupt_loose = (TH1F*)f->Get("h_bkgtaupt_loose");
	TH1F* h_bkgtaupt_medium = (TH1F*)f->Get("h_bkgtaupt_medium");
	TH1F* h_bkgtaupt_tight = (TH1F*)f->Get("h_bkgtaupt_tight");

	const int nBins_eff = h_taueff_loose->GetNbinsX();
	const int nBins_bkg = h_bkgtaupt_loose->GetNbinsX();

	float eff_noniso[nBins_eff];
	float eff_loose[nBins_eff];
	float eff_medium[nBins_eff];
	float eff_tight[nBins_eff];
	float eff_err_noniso[nBins_eff];
	float eff_err_loose[nBins_eff];
	float eff_err_medium[nBins_eff];
	float eff_err_tight[nBins_eff];

	float bkg_noniso[nBins_bkg];
	float bkg_loose[nBins_bkg];
	float bkg_medium[nBins_bkg];
	float bkg_tight[nBins_bkg];
	float bkg_err_noniso[nBins_bkg];
	float bkg_err_loose[nBins_bkg];
	float bkg_err_medium[nBins_bkg];
	float bkg_err_tight[nBins_bkg];

	for (int i = 0; i < nBins_eff; ++i ) {
		eff_noniso[i] = h_taueff_noniso->GetBinContent(i+1);
		eff_err_noniso[i] = h_taueff_noniso->GetBinError(i+1);
		eff_loose[i] = h_taueff_loose->GetBinContent(i+1);
		eff_err_loose[i] = h_taueff_loose->GetBinError(i+1);
		eff_medium[i] = h_taueff_medium->GetBinContent(i+1);
		eff_err_medium[i] = h_taueff_medium->GetBinError(i+1);
		eff_tight[i] = h_taueff_tight->GetBinContent(i+1);
		eff_err_tight[i] = h_taueff_tight->GetBinError(i+1);
	}

	for (int i = 0; i < nBins_bkg; ++i ) {
		bkg_noniso[i] = h_bkgtaupt_noniso->GetBinContent(i+1);
		bkg_err_noniso[i] = h_bkgtaupt_noniso->GetBinError(i+1);
		bkg_loose[i] = h_bkgtaupt_loose->GetBinContent(i+1);
		bkg_err_loose[i] = h_bkgtaupt_loose->GetBinError(i+1);
		bkg_medium[i] = h_bkgtaupt_medium->GetBinContent(i+1);
		bkg_err_medium[i] = h_bkgtaupt_medium->GetBinError(i+1);
		bkg_tight[i] = h_bkgtaupt_tight->GetBinContent(i+1);
		bkg_err_tight[i] = h_bkgtaupt_tight->GetBinError(i+1);
	}
	
	TCanvas c;

	//gPad->SetLogy();
    gPad->SetGrid();
	
	TGraphErrors* gr_noniso =
		new TGraphErrors(nBins_eff, eff_noniso, bkg_noniso, eff_err_noniso, bkg_err_noniso);
	TGraphErrors* gr_loose =
		new TGraphErrors(nBins_eff, eff_loose, bkg_loose, eff_err_loose, bkg_err_loose);
	TGraphErrors* gr_medium =
		new TGraphErrors(nBins_eff, eff_medium, bkg_medium, eff_err_medium, bkg_err_medium);
	TGraphErrors* gr_tight =
		new TGraphErrors(nBins_eff, eff_tight, bkg_tight, eff_err_tight, bkg_err_tight);

	gr_noniso->SetLineColor(kBlack);
	gr_noniso->SetMarkerColor(kBlack);
	gr_loose->SetLineColor(kGreen);
    gr_loose->SetMarkerColor(kGreen);
	gr_medium->SetLineColor(kBlue);
    gr_medium->SetMarkerColor(kBlue);
	gr_tight->SetLineColor(kRed);
    gr_tight->SetMarkerColor(kRed);
	
	gr_loose->SetTitle("Background vs signal efficiency");
	gr_loose->GetXaxis()->SetTitle("Signal efficiency #tau_{h}");
	gr_loose->GetYaxis()->SetTitle("Background (Number of fake taus)");
	gr_loose->SetMarkerStyle(7);
	gr_medium->SetMarkerStyle(7);
	gr_tight->SetMarkerStyle(7);
	gr_noniso->SetMarkerStyle(7);

	gr_loose->Draw("APZ0");
	gr_medium->Draw("same PZ0");
	gr_tight->Draw("same PZ0");
	gr_noniso->Draw("same PZ0");
	
	TLegend* leg=new TLegend(0.23,0.7,0.50,0.9);
	leg->AddEntry(gr_noniso,"nonIso","p");
	leg->AddEntry(gr_loose,"loose","p");
	leg->AddEntry(gr_medium,"medium","p");
	leg->AddEntry(gr_tight,"tight","p");
	leg->Draw("same");

	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/tauROC.pdf");

	delete f;
}

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
	h_eff_Tau_loose_pt -> SetLineColor(1);
	h_eff_Tau_loose_pt -> Draw("E");

	h_eff_Tau_medium_pt -> SetLineColor(2);
	h_eff_Tau_tight_pt  -> SetLineColor(3);
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
	h_eff_Tau_loose_eta -> SetLineColor(1);
	h_eff_Tau_loose_eta -> Draw("E");

	h_eff_Tau_medium_eta -> SetLineColor(2);
	h_eff_Tau_tight_eta  -> SetLineColor(3);
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
	h_eff_Tau_loose_phi -> SetLineColor(1);
	h_eff_Tau_loose_phi -> Draw("E");

	h_eff_Tau_medium_phi -> SetLineColor(2);
	h_eff_Tau_tight_phi  -> SetLineColor(3);
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
	h_nTau_tight->SetLineColor(3);
	h_nTau_tight->Draw("E");

	h_nTau_noniso -> SetLineColor(5);
	h_nTau_loose  -> SetLineColor(1);
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



void HistDrawer() {
  
  /// MC Truth plots
  //Draw_mtopmtautau();

  /// Tau selection efficiency 
  TauEffPlot();

  /// Event Selection plots
  //MakeTauROCPlot();                                                                                                                                                                               
}
