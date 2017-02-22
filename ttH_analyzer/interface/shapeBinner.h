#if defined(__ROOTCLING__) || defined(__ACLIC__)

#ifndef shapeBinner_h
#define shapeBinner_h

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TString.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TMath.h"
#include "TArrayD.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>

class shapeBinner
{
 public:

	// constructor and destructor
	shapeBinner(float, float, TString, bool uniform=false, int nbins=-1);
	shapeBinner(float, float, TFile*, bool uniform=false, int nbins=-1);
	~shapeBinner();

	// member function	
	void rebinHistograms(int verbosity=1);
	std::vector<double> getBinEdges();
	std::vector<double> getPurities();
	std::vector<float> getSignificance();
	std::vector<double> getRelativeErrors(int);
	void printResultsAll();
	int getNbins(){return _nbins;};
	
 protected:

	// binning
	std::vector<double> computeBinEdges(TH1*, TH1*, TH1*);
	std::vector<double> makeUniformBins(int, double, double);

	// bin analyzer
	void analyzeBins(TH1*, TH1*, TH1*);
	void computeSignificance(float,float);

	// utilities
	std::vector<TH1*> getHistograms(TFile*);
	double addBinErrors(double, double);
	double square(double);

	// rebin function from C. Veelken
	TH1* getRebinnedHistogram1d(const TH1*, const TArrayD&);
	// overload
	TH1* getRebinnedHistogram1d(const TH1*, std::vector<double>);
	
	// functions for removing negative bins from C. Veelken
	double compIntegral(TH1*, bool, bool);
	void makeBinContentsPositive(TH1*, int);
	
 private:

	std::vector<TH1*> _fine_datacards;
	TFile* _inputfile;
	
	float _relErrThreshold_bkg1;
	float _relErrThreshold_bkg2;

	std::vector<double> _binEdges;
	std::vector<double> _purities;
	std::vector<double> _relErrors_sig;
	std::vector<double> _relErrors_bkg1;
	std::vector<double> _relErrors_bkg2;

	bool _makeUniformBins;
	int _nbins;

	std::vector<float> _pvalue;
	std::vector<float> _punzi;
	std::vector<float> _approxpunzi;
	std::vector<float> _SoverSqrtB;
	std::vector<float> _SoverSqrtSplusB;
};

#endif
#endif
