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


void GenPlotter(const TString input =
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

	TH2F* h_ctH_ctbarH_lab = (TH2F*)f->Get("ctH_ctbarH_lab");
	TH2F* h_ctH_ctbarH_com = (TH2F*)f->Get("ctH_ctbarH_com");
	TH2F* h_dRtH_dRtbarH_lab = (TH2F*)f->Get("dRtH_dRtbarH_lab");
	TH2F* h_dRtH_dRtbarH_com = (TH2F*)f->Get("dRtH_dRtbarH_com");

	TH1F* h_top_higgs_dRmin = (TH1F*)f->Get("dRmin");
	TH1F* h_top_higgs_dRmax = (TH1F*)f->Get("dRmax");

	TCanvas c;
	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(10);
	
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

	h_ctH_ctbarH_lab->GetXaxis()->SetTitle("cos#theta(tbar,H)");
	h_ctH_ctbarH_lab->GetYaxis()->SetTitle("cos#theta(t,H)");
	h_ctH_ctbarH_lab->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/ctH_ctbarH_lab.pdf");

	h_ctH_ctbarH_com->GetXaxis()->SetTitle("cos#theta(tbar,H)");
	h_ctH_ctbarH_com->GetYaxis()->SetTitle("cos#theta(t,H)");
	h_ctH_ctbarH_com->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/ctH_ctbarH_com.pdf");

	h_dRtH_dRtbarH_lab->GetXaxis()->SetTitle("dR(tbar,H)");
	h_dRtH_dRtbarH_lab->GetYaxis()->SetTitle("dR(t,H)");
	h_dRtH_dRtbarH_lab->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dRtH_dRtbarH_lab.pdf");

	h_dRtH_dRtbarH_com->GetXaxis()->SetTitle("dR(tbar,H)");
	h_dRtH_dRtbarH_com->GetYaxis()->SetTitle("dR(t,H)");
	h_dRtH_dRtbarH_com->Draw("colz");
	c.SaveAs("/uscms/home/ztao/work/CU_ttH_WD/Outputs/dRtH_dRtbarH_com.pdf");
	
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

void HistDrawer() {
  
  /// MC Truth plots
  GenPlotter();
                                                                                        
}
