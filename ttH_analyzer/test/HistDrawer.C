#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;

void HistDrawer(const TString file = "/uscms/home/ztao/work/CU_ttH_WD/Outputs/histograms.root") {
	// Open root file and get histograms
	TFile* f = new TFile(file);
	
	TH2F* h_mtautau_mtop1 = (TH2F*)f->Get("mTTmtop1");
	TH2F* h_mtautau_mtop2 = (TH2F*)f->Get("mTTmtop2");
	TH2F* h_ptautau_mtop1 = (TH2F*)f->Get("pTTmtop1");
	TH2F* h_ptautau_mtop2 = (TH2F*)f->Get("pTTmtop2");

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

	delete f;
}
