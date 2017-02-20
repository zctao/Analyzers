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
#include <algorithm>

class shapeBinner
{
 public:

	// constructor and destructor
	shapeBinner(float, float, TString);
	~shapeBinner();

	// member function
	std::vector<TH1*> getHistograms(TFile*);	
	void renameHistograms();
	std::vector<double> computeBinEdges(TH1*, TH1*, TH1*);
	void rebinHistograms();
	std::vector<double> showBinEdges();
	
 protected:

	double addBinErrors(double, double);
	
 private:

	float _relErrThreshold_bkg1;
	float _relErrThreshold_bkg2;

	std::vector<double> _binEdges;
	std::vector<TH1*> _fine_datacards;
	
	TFile* _inputfile;
	
};

#endif
