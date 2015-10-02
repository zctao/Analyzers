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

using namespace std;

void Draw_mtopmtautau(const TString input = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/histograms.root" ) {
	// Open root file and get histograms
	TFile* f = new TFile(input);
	
	TH2F* h_mtautau_mtop1 = (TH2F*)f->Get("mTTmtop1");
	TH2F* h_mtautau_mtop2 = (TH2F*)f->Get("mTTmtop2");
	TH2F* h_ptautau_mtop1 = (TH2F*)f->Get("pTTmtop1");
	TH2F* h_ptautau_mtop2 = (TH2F*)f->Get("pTTmtop2");
	TH2F* h_mtop1TT_mtop2TT = (TH2F*)f->Get("mtop1TT_mtop2TT");
	TH2F* h_mplus_mminus = (TH2F*)f->Get("mplus_mminus");
	TH2F* h_mttTT_mplus = (TH2F*)f->Get("mttTT_mplus");

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

	h_mttTT_mplus->GetXaxis()->SetTitle("m_{top+top+#tau+#tau}");
	h_mttTT_mplus->GetYaxis()->SetTitle("max(m_{top1+#tau+#tau},m_{top2+#tau+#tau})");
	h_mttTT_mplus->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/mttTT_mplus.pdf");
	
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

void TauEff(const TString input = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/taueff.root")
{
	// Open root file and get histograms
	TFile* f = new TFile(input);

	TH1F* h_taus_selected_loose = (TH1F*)f->Get("h_taus_selected_loose");
	TH1F* h_taus_selected_medium = (TH1F*)f->Get("h_taus_selected_medium");
	TH1F* h_taus_selected_tight = (TH1F*)f->Get("h_taus_selected_tight");
	TH1F* h_taus_gen = (TH1F*)f->Get("h_taus_gen");

	TH1F* h_taus_selected_loose_eta = (TH1F*)f->Get("h_taus_selected_loose_eta");
	TH1F* h_taus_selected_medium_eta = (TH1F*)f->Get("h_taus_selected_medium_eta");
	TH1F* h_taus_selected_tight_eta = (TH1F*)f->Get("h_taus_selected_tight_eta");
	TH1F* h_taus_gen_eta = (TH1F*)f->Get("h_taus_gen_eta");
	
	TCanvas c;
	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(10);
	
	h_taus_selected_tight->SetTitle("ttH, H->#tau#tau (M125_13TeV_Spring15)");
	h_taus_selected_tight->GetXaxis()->SetTitle("n_{#tau_{h}}");
	h_taus_selected_tight->GetYaxis()->SetTitle("Fraction of events / 1.0");
	h_taus_selected_tight->SetLineColor(1);
	h_taus_selected_tight->Draw("E");
	
	h_taus_selected_medium->SetLineColor(2);
	h_taus_selected_loose->SetLineColor(3);
	h_taus_gen->SetLineColor(4);

	h_taus_selected_medium->Draw("E same");
	h_taus_selected_loose->Draw("E same");
	h_taus_gen->Draw("E same");
	
	TLegend* leg = new TLegend(0.6,0.68,0.88,0.85);
	leg->AddEntry((TObject*)0, "p_{T}>30 GeV", "");
	leg->AddEntry(h_taus_selected_tight,"PAT tight","l");
	leg->AddEntry(h_taus_selected_medium,"PAT medium","l");
	leg->AddEntry(h_taus_selected_loose,"PAT loose","l");
	leg->AddEntry(h_taus_gen,"GEN","l");
	leg->Draw("same");

	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/tauEff.pdf");

	TCanvas c2;
	h_taus_selected_tight_eta->SetTitle("ttH, H->#tau#tau (M125_13TeV_Spring15)");
	h_taus_selected_tight_eta->GetXaxis()->SetTitle("n_{#tau_{h}}");
	h_taus_selected_tight_eta->GetYaxis()->SetTitle("Fraction of events / 1.0");
	h_taus_selected_tight_eta->SetLineColor(1);
	h_taus_selected_tight_eta->Draw("E");
	
	h_taus_selected_medium_eta->SetLineColor(2);
	h_taus_selected_loose_eta->SetLineColor(3);
	h_taus_gen_eta->SetLineColor(4);

	h_taus_selected_medium_eta->Draw("E same");
	h_taus_selected_loose_eta->Draw("E same");
	h_taus_gen_eta->Draw("E same");
	
	TLegend* leg2 = new TLegend(0.6,0.68,0.88,0.85);
	leg2->AddEntry((TObject*)0, "p_{T}>25 GeV  |#eta|<2.3", "");
	leg2->AddEntry(h_taus_selected_tight_eta,"PAT tight","l");
	leg2->AddEntry(h_taus_selected_medium_eta,"PAT medium","l");
	leg2->AddEntry(h_taus_selected_loose_eta,"PAT loose","l");
	leg2->AddEntry(h_taus_gen_eta,"GEN","l");
	leg2->Draw("same");

	c2.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/tauEff_eta.pdf");
	
	delete f;
}

void HistDrawer() {
	Draw_mtopmtautau();
	//TauEff();
}
