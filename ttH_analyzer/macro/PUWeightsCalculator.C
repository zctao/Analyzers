#ifndef PUWeightsCalculator_c
#define PUWeightsCalculator_c

#include "TH1.h"
#include "TFile.h"

#include <iostream>

void PUWeightsCalculator(TString fname_data_pu = "../data/MyDataPileupHistogram.root",
					TString fname_out = "PU_weights_2016_ReReco_MCMoriond17_271036_284044.root")
{

	using namespace std;
	
	// MC pileup scenario (Moriond 17)
	// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py
	int MC_Pileup_bin[] =
		{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74};
	double MC_Pileup_value[] =
		{1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,1.73032e-05 ,1.435e-05 ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05};

	assert(sizeof(MC_Pileup_bin)/sizeof(MC_Pileup_bin[0])==
		   sizeof(MC_Pileup_value)/sizeof(MC_Pileup_value[0]));

	TH1D* h_mc_pileup = new TH1D("h_MC_M17","",60,0,60);

	for (int i = 0; i < 60; ++i) {
		int ibin = MC_Pileup_bin[i];
		double puValue = MC_Pileup_value[i];
		h_mc_pileup -> SetBinContent(ibin, puValue);
	}

	// Data pileup scenario
	TFile* f_data = new TFile(fname_data_pu);
	TH1D* h_data_pileup = (TH1D*)f_data->Get("pileup");
	h_data_pileup -> SetName("h_data_pileup");
	h_data_pileup->SetDirectory(0);

	// close file
	f_data->Close();
	delete f_data;

	// normalize
	double data_scale = 1./ h_data_pileup -> Integral();
	h_data_pileup->Scale(data_scale);

	TFile* f_out = new TFile(fname_out, "RECREATE");

	int nbins = h_mc_pileup -> GetNbinsX();
	assert(h_data_pileup->GetNbinsX()==nbins);

	TH1D* h_ratio_data_mc = new TH1D("h_ratio_data_MC","",nbins,0,nbins);

	for (int i = 0; i < nbins; ++i) {
		double puData = h_data_pileup->GetBinContent(i);
		double puMC = h_mc_pileup->GetBinContent(i);
		double weight = puMC>0. ? puData/puMC : 0.;

		h_ratio_data_mc -> SetBinContent(i, weight);
	}
		
	h_mc_pileup -> Write();
	h_data_pileup -> Write();
	h_ratio_data_mc -> Write();

}

#endif
