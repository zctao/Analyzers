#ifndef _Cross_Sections_h
#define _Cross_Sections_h
#include <map>
#include <string>

namespace xsection{  // unit: pb
	std::map<std::string, float> xsection =
		{{"ttH",0.2586},          // ttH, H->nonbb
		 {"TTW",0.2043},          // TTWJetsToLNu
		 {"TTZ",0.2529},          // TTZToLLNuNu
		 {"TTJets_ll",87.3},      // TTJets_DiLept
		 {"TTJets_lt",182},       // TTJets_SingleLeptFromT
		 {"TTJets_ltbar",182}     // TTJets_SingleLeptFromTbar
		};
}

#endif
