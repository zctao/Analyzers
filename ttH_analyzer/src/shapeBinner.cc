#ifndef shapeBinner_cc
#define shapeBinner_cc

//#include "Analyzers/ttH_analyzer/interface/shapeBinner.h"
#include "../interface/shapeBinner.h"

shapeBinner::shapeBinner(float relErr_bkg1, float relErr_bkg2, TFile* inputfile, bool uniform, int nbins)
{
	_relErrThreshold_bkg1 = relErr_bkg1;
	_relErrThreshold_bkg2 = relErr_bkg2;
	
	_fine_datacards = getHistograms(inputfile);
	_makeUniformBins = uniform;
	_nbins = nbins;
}

shapeBinner::shapeBinner(float relErr_bkg1, float relErr_bkg2, TString filename, bool uniform, int nbins)
{
	_relErrThreshold_bkg1 = relErr_bkg1;
	_relErrThreshold_bkg2 = relErr_bkg2;

	TFile* _inputfile = new TFile(filename, "read");
	_fine_datacards = getHistograms(_inputfile);
	
	_makeUniformBins = uniform;
	_nbins = nbins;
	
}

shapeBinner::~shapeBinner()
{
	if (_inputfile)
		delete _inputfile;
}

template<typename T>
std::ostream &operator<< (ostream& out, const std::vector<T>& v)
{
	out << "{ ";
	size_t last = v.size() - 1;
	for (size_t i = 0; i < v.size(); ++i) {
		out << v[i];
		if (i != last)
			out << ", ";
	}
	out << " }";
	return out;
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
			/*
			// purity
			_purities.push_back(binContent_sig / 
								(binContent_bkg1+binContent_bkg2));
			// significance
			computeSignificance(binContent_sig, binContent_bkg1+binContent_bkg2);
			*/
			
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
	}
	
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

std::vector<double> shapeBinner::makeUniformBins(int nbins, double xlow, double xup)
{
	std::vector<double> binEdges;

	for (int i = 0; i <= nbins; ++i) {
		binEdges.push_back(xlow + i*(xup-xlow)/nbins);
	}

	return binEdges;
}

void shapeBinner::rebinHistograms(int verbosity)
{
	if (verbosity)
		std::cout << "fine_datacards size : " << _fine_datacards.size()
				  << std::endl;
	
	// Get signal, reducible and irreducible histograms
	
	THStack hs_signal("sig","");
	THStack hs_bkg_irreducible("bkg_irr","");
	THStack hs_bkg_reducible("bkg_red","");

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
			hs_bkg_reducible.Add(h);
		}
	}

	TH1* h_signal = (TH1*)(hs_signal.GetStack()->Last())->Clone();
	TH1* h_bkg_irreducible = (TH1*)(hs_bkg_irreducible.GetStack()->Last())->Clone();
	TH1 *h_bkg_reducible = (TH1*)(hs_bkg_reducible.GetStack()->Last())->Clone();

	if (verbosity)
		std::cout << "start computing bin edges" << std::endl;

	if (not _makeUniformBins) {
		_binEdges = computeBinEdges(h_signal, h_bkg_irreducible, h_bkg_reducible);
		_nbins = _binEdges.size()-1;
	}
	else {
		double xlow = h_signal->GetXaxis()->GetXmin();
		double xhigh = h_signal->GetXaxis()->GetXmax();
		//std::cout << xlow << " " << xhigh << std::endl;
		assert(_nbins>0 and _makeUniformBins);
		_binEdges = makeUniformBins(_nbins, xlow, xhigh);
	}

	if (verbosity)
		std::cout << "number of bins : " << _nbins << std::endl;

	TH1* h_sig_rebin = h_signal->Rebin(_nbins, "", &_binEdges[0]);
	TH1* h_bkg1_rebin = h_bkg_irreducible->Rebin(_nbins, "", &_binEdges[0]);
	TH1* h_bkg2_rebin = h_bkg_reducible->Rebin(_nbins, "", &_binEdges[0]);
	analyzeBins(h_sig_rebin, h_bkg1_rebin, h_bkg2_rebin);
	
	delete h_signal;
	delete h_bkg_reducible;
	delete h_bkg_irreducible;
	delete h_sig_rebin;
	delete h_bkg1_rebin;
	delete h_bkg2_rebin;
	
	// rebin all shapes with the new bin edges and output to root file
	TString outname = "rebinned_datacards_"+TString::Itoa(_nbins,10)+"bins_";

	if (not _makeUniformBins) {
		ostringstream param1;
		ostringstream param2;
		param1 << _relErrThreshold_bkg1;
		param2 << _relErrThreshold_bkg2;
		
		outname += param1.str()+"_"+param2.str()+".root";
	}
	else {
		outname += "uniform.root";
	}

	TFile *output = new TFile(outname, "recreate");
	
	for (auto h : _fine_datacards) {
		TString hname = h->GetName();
		
		TH1* h_rebin = h->Rebin(_nbins, hname, &_binEdges[0]);
		//TH1* h_rebin = getRebinnedHistogram1d(h, _binEdges);
		
		makeBinContentsPositive(h_rebin,0);
		h_rebin->Write();
	}
	
	if (verbosity)
		std::cout << "Output file : " << outname << std::endl;
	
	output->Close();
	
	return;
}

void shapeBinner::analyzeBins(TH1* h_sig, TH1* h_bkg1, TH1* h_bkg2)
{
	_purities.clear();
	_relErrors_sig.clear();
	_relErrors_bkg1.clear();
	_relErrors_bkg2.clear();
	_pvalue.clear();
	
	int nbins = h_sig->GetNbinsX();
	assert(nbins == h_bkg1->GetNbinsX() and nbins == h_bkg2->GetNbinsX());

	for (int ibin = 1; ibin <= nbins; ++ibin) {
		double binContent_sig  = h_sig  -> GetBinContent(ibin);
		double binContent_bkg1 = h_bkg1 -> GetBinContent(ibin);
		double binContent_bkg2 = h_bkg2 -> GetBinContent(ibin);
		double binError_sig  = h_sig  -> GetBinError(ibin);
		double binError_bkg1 = h_bkg1 -> GetBinError(ibin);
		double binError_bkg2 = h_bkg2 -> GetBinError(ibin);

		// signal purity
		double binContent_sum = binContent_sig+binContent_bkg1+binContent_bkg2;
		_purities.push_back(binContent_sig / binContent_sum);

		// relative errors
		double relerr_sig =
			(binContent_sig > 0) ? (binError_sig/binContent_sig) : 0;
		double relerr_bkg1 =
			(binContent_bkg1 > 0) ? (binError_bkg1/binContent_bkg1) : 0;
		double relerr_bkg2 =
			(binContent_bkg2 > 0) ? (binError_bkg2/binContent_bkg2) : 0;
		
		_relErrors_sig.push_back(relerr_sig);
		_relErrors_bkg1.push_back(relerr_bkg1);
		_relErrors_bkg2.push_back(relerr_bkg2);

		// significance
		computeSignificance(binContent_sig, binContent_bkg1 + binContent_bkg2);
	}
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

std::vector<double> shapeBinner::getRelativeErrors(int ch)
{
	switch(ch) {
	case 0:  // signal
		return _relErrors_sig;
		break;
	case 1:  // irreducible background
		return _relErrors_bkg1;
		break;
	case 2:  // reducible background
		return _relErrors_bkg2;
		break;
	default:
		std::cerr << "getRelativeErrors: Not valid option!" << std::endl;
		std::vector<double> empty;
		return empty;
	}
}

std::vector<float> shapeBinner::getSignificance()
{
	return _pvalue;
}

void shapeBinner::printResultsAll()
{
	using std::cout;
	using std::endl;

	if (not _makeUniformBins) {
		ostringstream param1;
		ostringstream param2;
		param1 << _relErrThreshold_bkg1;
		param2 << _relErrThreshold_bkg2;
		cout << "## " << param1.str() << " " << param2.str() << " => " << _nbins
			 << " bins" << endl;
	}
	else
		cout << "uniform " << _nbins << " bins" << endl;
	
	cout << "Bin edges ";
	cout << _binEdges << endl;

	cout << "Signal purities ";
	cout << _purities << endl;
	
	cout << "Relative errors " << endl;
	cout << "Signal ";
	cout << _relErrors_sig << endl;
	cout << "Irreducible background ";
	cout << _relErrors_bkg1 << endl;
	cout << "Reducible background ";
	cout << _relErrors_bkg2 << endl;
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
