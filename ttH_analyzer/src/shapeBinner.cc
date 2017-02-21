#ifndef shapeBinner_cc
#define shapeBinner_cc

//#include "Analyzers/ttH_analyzer/interface/shapeBinner.h"
#include "../interface/shapeBinner.h"

shapeBinner::shapeBinner(float relErr_bkg1, float relErr_bkg2, TString filename, bool extraBinningUniform)
{
	_relErrThreshold_bkg1 = relErr_bkg1;
	_relErrThreshold_bkg2 = relErr_bkg2;

	TFile* _inputfile = new TFile(filename, "read");
	_fine_datacards = getHistograms(_inputfile);
	_makeUniformBins = extraBinningUniform;

	_purities.clear();
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
	std::vector<double> purites_rebinned;

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

		bool goodBkg1 =
			(binError_bkg1 < binContent_bkg1 * _relErrThreshold_bkg1) and
			(binContent_bkg1 > 0);
		bool goodBkg2 =
			(binError_bkg2 < binContent_bkg2 * _relErrThreshold_bkg2) and
			(binContent_bkg2 > 0);
		bool goodSig = binContent_sig > 0;
		
		if (goodBkg1 and goodBkg2 and goodSig) {
			
			binEdges_rebinned.push_back(lowedge);  // reversed order

			// purity
			_purities.push_back(binContent_sig / 
								(binContent_bkg1+binContent_bkg2));
			// significance
			computeSignificance(binContent_sig, binContent_bkg1+binContent_bkg2);
			
			// reset
			binError_sig = 0.;
			binError_bkg1 = 0.;
			binError_bkg2 = 0.;
			binContent_sig = 0.;
			binContent_bkg1 = 0.;
			binContent_bkg2 = 0.;
		}
	}

	if (binEdges_orig.at(0)!=binEdges_rebinned.back()) {
		binEdges_rebinned.push_back(binEdges_orig.at(0));
		_purities.push_back(binContent_sig / 
							(binContent_bkg1+binContent_bkg2));
	}
	
	std::reverse(binEdges_rebinned.begin(),binEdges_rebinned.end());
	std::reverse(_purities.begin(),_purities.end());
	std::reverse(_pvalue.begin(),_pvalue.end());
	std::reverse(_SoverSqrtB.begin(),_SoverSqrtB.end());
	std::reverse(_SoverSqrtSplusB.begin(),_SoverSqrtSplusB.end());
	
	// add the upper edge
	double upedge = h_sig->GetXaxis()->GetBinUpEdge(nbins);
	binEdges_rebinned.push_back(upedge);
	
	return binEdges_rebinned;
}

double shapeBinner::addBinErrors(double binError1, double binError2)
{
	return TMath::Sqrt(binError1*binError1 + binError2*binError2);
}

std::vector<double> shapeBinner::makeUniformBins(int nbins, double xlow, double xup)
{
	std::vector<double> binEdges;

	for (int i = 0; i <= nbins; ++i) {
		binEdges.push_back(xlow + i*(xup-xlow)/nbins);
	}

	return binEdges;
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
		else if (hname.Contains("TTW") or hname.Contains("TTZ") or
				 hname.Contains("Rares") or hname.Contains("EWK")) {
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

	int nbins = _binEdges.size()-1;
	std::cout << "number of bins : " << nbins << std::endl;
	
	delete h_signal;
	delete h_bkg_reducible;
	delete h_bkg_irreducible;
	
	// rebin all shapes with the new bin edges and output to root file
	TString outname = "rebinned_datacards";

	ostringstream param1;
	ostringstream param2;
	param1 << _relErrThreshold_bkg1;
	param2 << _relErrThreshold_bkg2;

	outname += ("_"+TString::Itoa(nbins,10)+"bins_"+param1.str()+"_"+param2.str()+".root");

	TFile *output = new TFile(outname, "recreate");
	
	for (auto h : _fine_datacards) {
		TString hname = h->GetName();
		
		//TH1* h_rebin = h->Rebin(nbins, hname, &_binEdges[0]);
		TH1* h_rebin = getRebinnedHistogram1d(h, _binEdges);
		
		makeBinContentsPositive(h_rebin,0);
		h_rebin->Write();
	}
	std::cout << "Output file : " << outname << std::endl;
	
	output->Close();
	
	if (_makeUniformBins) {
		TString outname_u = "rebinned_datacards_uniform_"+TString::Itoa(nbins,10)+"bins.root";
		TFile *output_u = new TFile(outname_u,"recreate");
		
		std::vector<double> binEdgesUniform =
			makeUniformBins(nbins,_binEdges.at(0),_binEdges.back());

		for (auto h : _fine_datacards) {
			TString hname = h->GetName();
			TH1* h_rebin_u = h->Rebin(nbins, hname, &binEdgesUniform[0]);
			makeBinContentsPositive(h_rebin_u,0);
			h_rebin_u->Write();
		}
		std::cout << "Output file : " << outname_u << std::endl;
		
		output_u->Close();
	}
	
	return;
}

void shapeBinner::computeSignificance(float s, float b)
{
	_pvalue.push_back(TMath::Poisson(s+b,b));
	_SoverSqrtB.push_back(s/TMath::Sqrt(b));
	_SoverSqrtSplusB.push_back(s/TMath::Sqrt(s+b));
}

double shapeBinner::square(double x)
{
	return x*x;
}

std::vector<double> shapeBinner::getBinEdges()
{
	return _binEdges;
}

std::vector<double> shapeBinner::getPurities()
{
	return _purities;
}

std::vector<float> shapeBinner::getSignificance()
{
	return _pvalue;
}

TH1* shapeBinner::getRebinnedHistogram1d(const TH1* histoOriginal,
										 std::vector<double> binEdges_rebinned)
{
	int asize = binEdges_rebinned.size();
	const TArrayD BinEdges(asize, &binEdges_rebinned[0]);

	return getRebinnedHistogram1d(histoOriginal, BinEdges);
}

TH1* shapeBinner::getRebinnedHistogram1d(const TH1* histoOriginal,
										 const TArrayD& binEdges_rebinned)
{
	std::string histoRebinnedName = std::string(histoOriginal->GetName());//.append("_rebinned");
	TH1* histoRebinned = new TH1D(histoRebinnedName.data(),
								  histoOriginal->GetTitle(), 
								  binEdges_rebinned.GetSize() - 1,
								  binEdges_rebinned.GetArray());
	histoRebinned->Sumw2();
	const TAxis* axis_original = histoOriginal->GetXaxis();
	int numBins_original = axis_original->GetNbins();
	for ( int idxBin = 1; idxBin <= numBins_original; ++idxBin ) {
		double binContent_original = histoOriginal->GetBinContent(idxBin);
		double binError_original = histoOriginal->GetBinError(idxBin);
		double binCenter_original = axis_original->GetBinCenter(idxBin);
		int binIndex_rebinned = histoRebinned->FindBin(binCenter_original);
		double binContent_rebinned = histoRebinned->GetBinContent(binIndex_rebinned);
		binContent_rebinned += binContent_original;
		histoRebinned->SetBinContent(binIndex_rebinned, binContent_rebinned);   
		double binError_rebinned = histoRebinned->GetBinError(binIndex_rebinned);
		binError_rebinned = TMath::Sqrt(binError_rebinned*binError_rebinned + binError_original*binError_original);
		histoRebinned->SetBinError(binIndex_rebinned, binError_rebinned);
	}
	return histoRebinned;
}

// functions for removing negative bins from C. Veelken
double shapeBinner::compIntegral(TH1* histogram, bool includeUnderflowBin, bool includeOverflowBin)
{
	double sumBinContent = 0.;
	int numBins = histogram->GetNbinsX();
	int firstBin = ( includeUnderflowBin ) ? 0 : 1;
	int lastBin = ( includeOverflowBin  ) ? (numBins + 1) : numBins;
	
	for ( int iBin = firstBin; iBin <= lastBin; ++iBin ) {
		sumBinContent += histogram->GetBinContent(iBin);
	}
	
	return sumBinContent;
}

void shapeBinner::makeBinContentsPositive(TH1* histogram, int verbosity)

{
	if ( verbosity ) {
		std::cout << "<makeBinContentsPositive>:" << std::endl;
		std::cout << " integral(" << histogram->GetName() << ") = " << histogram->Integral() << std::endl;
	}
	
	double integral_original = compIntegral(histogram, true, true);
	
	if ( integral_original < 0. ) integral_original = 0.;
	
	if ( verbosity ) {
		std::cout << " integral_original = " << integral_original << std::endl;
	}
	
	int numBins = histogram->GetNbinsX();
	
	for ( int iBin = 0; iBin <= (numBins + 1); ++iBin ) {
		double binContent_original = histogram->GetBinContent(iBin);
		double binError2_original = square(histogram->GetBinError(iBin));
		
		if ( binContent_original < 0. ) {
			double binContent_modified = 0.;
			double binError2_modified = binError2_original + square(binContent_original - binContent_modified);
			
			assert(binError2_modified >= 0.);
			
			if ( verbosity ) {
				std::cout << "bin #" << iBin << " (x =  " << histogram->GetBinCenter(iBin) << "): binContent = " << binContent_original << " +/- " << TMath::Sqrt(binError2_original) << " --> setting it to binContent = " << binContent_modified << " +/- " << TMath::Sqrt(binError2_modified) << std::endl;
			}

			histogram->SetBinContent(iBin, binContent_modified);
			histogram->SetBinError(iBin, TMath::Sqrt(binError2_modified));
		}
	}

	double integral_modified = compIntegral(histogram, true, true);
	
	if ( integral_modified < 0. ) integral_modified = 0.;
	
	if ( verbosity ) {
		std::cout << " integral_modified = " << integral_modified << std::endl;
	}
	
	if ( integral_modified > 0. ) {
		double sf = integral_original/integral_modified;
		
		if ( verbosity ) {
			std::cout << "--> scaling histogram by factor = " << sf << std::endl;
		}
		
		histogram->Scale(sf);
	} else {
		for ( int iBin = 0; iBin <= (numBins + 1); ++iBin ) {
			histogram->SetBinContent(iBin, 0.);
		}
	}

	if ( verbosity ) {
		std::cout << " integral(" << histogram->GetName() << ") = " << histogram->Integral() << std::endl;
	}
}

#endif
