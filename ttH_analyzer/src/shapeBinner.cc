#ifndef shapeBinner_cc
#define shapeBinner_cc

//#include "Analyzers/ttH_analyzer/interface/shapeBinner.h"
#include "../interface/shapeBinner.h"

shapeBinner::shapeBinner(float relErr_bkg1, float relErr_bkg2, TString filename)
{
	_relErrThreshold_bkg1 = relErr_bkg1;
	_relErrThreshold_bkg2 = relErr_bkg2;

	TFile* _inputfile = new TFile(filename, "read");
	_fine_datacards = getHistograms(_inputfile);
	
}

shapeBinner::~shapeBinner()
{
	delete _inputfile;
}

std::vector<TH1*> shapeBinner::getHistograms(TFile* f)
{
	std::vector<TH1*> results;
	
	TIter next(f->GetListOfKeys());
	TKey* key;

	while ((key = (TKey*)next())) {
		TClass* c1 = gROOT->GetClass(key->GetClassName());
		if (!c1->InheritsFrom("TH1")) continue;

		TString hname = key->GetName();
		if (!hname.BeginsWith("TMVA_fine_inclusive_")) continue;

		TH1* h = (TH1*) key->ReadObj();
		hname.ReplaceAll("TMVA_fine_inclusive_","x_");
		h->SetName(hname);
		results.push_back(h);
	}

	return results;
}

void shapeBinner::renameHistograms()
{

}

std::vector<double> shapeBinner::computeBinEdges(TH1* h_sig, TH1* h_bkg1, TH1* h_bkg2)
{
	// get current bin edges
	int nbins = h_sig->GetNbinsX();
	assert(nbins==h_bkg1->GetNbinsX() and nbins==h_bkg2->GetNbinsX());		
	std::vector<double> binEdges_orig;
	binEdges_orig.resize(nbins);
	h_sig->GetXaxis()->GetLowEdge(&binEdges_orig[0]);

	std::vector<double> binEdges_rebinned;

	double binError_sig = 0., binContent_sig = 0.;
	double binError_bkg1 = 0., binContent_bkg1 = 0.;
	double binError_bkg2 = 0., binContent_bkg2 = 0.;
	
	for (int ibin = nbins; ibin > 0; ibin--) { // start from right most bin 

		double lowedge = binEdges_orig.at(ibin-1);
	
		binContent_sig += h_sig->GetBinContent(ibin);		
		binContent_bkg1 += h_bkg1->GetBinContent(ibin);		
		binContent_bkg2 += h_bkg2->GetBinContent(ibin);
	
		binError_sig = addBinErrors(binError_sig, h_sig->GetBinError(ibin));
		binError_bkg1 = addBinErrors(binError_bkg1, h_bkg1->GetBinError(ibin));
		binError_bkg2 = addBinErrors(binError_bkg2, h_bkg2->GetBinError(ibin));

		bool goodBkg1 = binError_bkg1 < binContent_bkg1 * _relErrThreshold_bkg1;
		bool goodBkg2 = binError_bkg2 < binContent_bkg2 * _relErrThreshold_bkg2;
		
		if (goodBkg1 and goodBkg2) {
			
			binEdges_rebinned.push_back(lowedge);  // reversed order
			
			binError_sig = 0.;
			binError_bkg1 = 0.;
			binError_bkg2 = 0.;
			binContent_sig = 0.;
			binContent_bkg1 = 0.;
			binContent_bkg2 = 0.;
		}
	}

	if (binEdges_orig.at(0)!=binEdges_rebinned.back())
		binEdges_rebinned.push_back(binEdges_orig.at(0));
	
	std::reverse(binEdges_rebinned.begin(),binEdges_rebinned.end());

	// add the upper edge
	double upedge = h_sig->GetXaxis()->GetBinUpEdge(nbins);
	binEdges_rebinned.push_back(upedge);
	
	return binEdges_rebinned;
}

double shapeBinner::addBinErrors(double binError1, double binError2)
{
	return TMath::Sqrt(binError1*binError1 + binError2*binError2);
}

void shapeBinner::rebinHistograms()
{
	std::cout << "fine_datacards size : " << _fine_datacards.size() << std::endl;
	
	// Get signal, reducible and irreducible histograms
	
	THStack hs_signal("sig","");
	THStack hs_bkg_irreducible("bkg_irr","");	
	TH1 *h_bkg_reducible;

	for (auto h : _fine_datacards) {
		TString hname = h->GetName();

		if (hname.Contains("_CMS_ttHl_")) continue; // systematics
		
		if (hname.Contains("ttH")) {
			hs_signal.Add(h);
		}
		else if (hname.Contains("TTW") or hname.Contains("TTZ")) {
			hs_bkg_irreducible.Add(h);
		}
		else if (hname.Contains("fakes")) {
			h_bkg_reducible = (TH1*)h->Clone();
		}
	}

	TH1* h_signal = (TH1*)(hs_signal.GetStack()->Last())->Clone();
	TH1* h_bkg_irreducible = (TH1*)(hs_bkg_irreducible.GetStack()->Last())->Clone();

	std::cout << "start computing bin edges" << std::endl;
    _binEdges = computeBinEdges(h_signal, h_bkg_irreducible, h_bkg_reducible);

	std::cout << "number of bins : " << _binEdges.size()-1 << std::endl;
	
	delete h_signal;
	delete h_bkg_reducible;
	delete h_bkg_irreducible;
	
	// rebin all shapes with the new bin edges and output to root file
	TFile *output = new TFile("rebinned_datacards.root", "recreate");
	
	int nbins = _binEdges.size()-1;
	for (auto h : _fine_datacards) {
		TString hname = h->GetName();
		TH1* h_rebin = h->Rebin(nbins, hname, &_binEdges[0]);
		h_rebin->Write();
	}

	return;
}

std::vector<double> shapeBinner::showBinEdges()
{
	return _binEdges;
}

#endif
