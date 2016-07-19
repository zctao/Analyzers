#ifndef kinMVA_ttV_2lss_h
#define kinMVA_ttV_2lss_h

#include "Analyzers/ttH_analyzer/interface/kinMVA_2lss.h"

class kinMVA_ttV_2lss : public kinMVA_2lss
{
 public:

	kinMVA_ttV_2lss();
	~kinMVA_ttV_2lss();

	void Set_up_Reader(TMVA::Reader*);

	float Get_lep1_conePt() const {return lep1_conePt;}
	float Get_lep2_conePt() const {return lep2_conePt;}

	void Calculate_mvaVars(const std::vector<miniLepton>&,
						   const std::vector<pat::Tau>&,
						   const std::vector<pat::Jet>&,
						   const pat::MET&);
	
 private:
	float lep1_conePt;
	float lep2_conePt;
};

#endif
