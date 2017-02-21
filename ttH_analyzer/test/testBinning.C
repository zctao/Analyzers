#include "TROOT.h"
#include "TString.h"
#include "TFile.h"

#include "../interface/shapeBinner.h"

#include <iostream>

void testBinning()
{
	using std::cout;
	using std::endl;
	
	//gROOT->ProcessLine(".L ../src/shapeBinner.cc+");

	float parameters[]={0.15, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26,0.27, 0.28, 0.35, 0.40, 0.50};

	TString inputfile = "limits.root";
	TFile* f = new TFile(inputfile, "read");
	
	for (auto p : parameters) {

		shapeBinner* sb = new shapeBinner(p,p,f);
		sb->rebinHistograms(0);
		cout << "-------------------------------------------------------" << endl;
		sb->printResultsAll();
		cout << "Expected limit : " << endl;

		int nbins = sb->getNbins();
		shapeBinner* sbu = new shapeBinner(p,p,f,true,nbins);
		sbu->rebinHistograms(0);
		cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
		sbu->printResultsAll();
		cout << "Expected limit : " << endl;
		
		//delete sb;
		//delete sbu;
	}
}
