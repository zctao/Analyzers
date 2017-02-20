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

#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>

class shapeBinner
{
 public:

	// constructor and destructor
	shapeBinner(float, float, TString, bool extraBinningUniform=false);
	~shapeBinner();

	// member function
	std::vector<TH1*> getHistograms(TFile*);	
	void renameHistograms();
	std::vector<double> computeBinEdges(TH1*, TH1*, TH1*);
	std::vector<double> makeUniformBins(int, double, double);
	void rebinHistograms();
	std::vector<double> getBinEdges();
	std::vector<double> getPurities();
	std::vector<float> getSignificance();
	
	// functions for removing negative bins from C. Veelken
	double compIntegral(TH1*, bool, bool);
	void makeBinContentsPositive(TH1*, int);
	
 protected:

	double addBinErrors(double, double);
	double square(double);
	void computeSignificance(float,float);
	
 private:

	std::vector<TH1*> _fine_datacards;
	TFile* _inputfile;
	
	float _relErrThreshold_bkg1;
	float _relErrThreshold_bkg2;

	std::vector<double> _binEdges;
	std::vector<double> _purities;

	bool _makeUniformBins;

	std::vector<float> _pvalue;
	std::vector<float> _punzi;
	std::vector<float> _approxpunzi;
	std::vector<float> _SoverSqrtB;
	std::vector<float> _SoverSqrtSplusB;
};

#endif
