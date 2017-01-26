#ifndef _Cross_Sections_h
#define _Cross_Sections_h
#include <map>
#include <string>

namespace xsection{  // unit: pb
	std::map<std::string, float> xsection =
		{{"ttH",0.215},          // ttH, H->nonbb
		 {"TTW",0.204},          // TTWJetsToLNu
		 {"TTZ",0.253},          // TTZToLLNuNu
		 {"TTGJets",3.70},       // TTGJets	 
		 {"TGJets",2.97},        // TGJets
		 {"WG",586},             // WGToLNuG
		 {"ZG",131},             // ZGTo2LG
		 {"WW",10.5},            // WWTo2L2Nu
		 {"WpWp",0.0371},        // WpWpJJ
		 {"WZ",4.10},            // WZTo3LNu
		 {"ZZ",1.26},            // ZZTo4L
		 {"WWW",0.209},          // WWW_4F
		 {"WWZ",0.165},          // WWZ
		 {"WZZ",0.0557},         // WZZ
		 {"ZZZ",0.014},          // ZZZ
		 {"tZq",0.0758},         // tZq_ll_4f
		 {"TTTT",0.0091},        // TTTT
		 {"TTJets_ll",87.3},     // TTJets_DiLept
		 {"TTJets_lt",182},      // TTJets_SingleLeptFromT
		 {"TTJets_ltbar",182},   // TTJets_SingleLeptFromTbar
		 {"DYJets_M10to50",18600},// DYJetsToLL_M-10to50
		 {"DYJets_M50",5770},    // DYJetsToLL_M-50
		 {"WJets",61500},        // WJetsToLNu
		 {"ST_sLep",3.75},       // ST_s-channel_4f_lep
		 {"ST_tT",70.7},         // ST_t-channel_top_4f
		 {"ST_tTbar",70.7},      // ST_t-channel_antitop_4f
		 {"ST_tWT",35.6},        // ST_tW_top_5f
		 {"ST_tWTbar",35.6},     // ST_tW_antitop_5f
		 {"WWds",} // WW_DoubleScattering
		};
}

#endif
