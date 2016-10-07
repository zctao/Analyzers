#ifndef kinMVA_ttbar_2lss_h
#define kinMVA_ttbar_2lss_h

#include "Analyzers/ttH_analyzer/interface/kinMVA_2lss.h"

class kinMVA_ttbar_2lss : public kinMVA_2lss
{
 public:

	kinMVA_ttbar_2lss();
	~kinMVA_ttbar_2lss();

	void Set_up_Reader(TMVA::Reader*);

	void Calculate_mvaVars(const std::vector<miniLepton>&,
						   const std::vector<pat::Tau>&,
						   const std::vector<pat::Jet>&,
						   const pat::MET&);

	float Get_avg_dr_jet() const {return avg_dr_jet;}
	float Get_MET() const {return met;}

	float Set_avg_dr_jet(const vector<pat::Jet>&);
	float Set_MET(double pfmet) { return min(pfmet,400.);}
	
 private:
	float avg_dr_jet;
	float met;
};

#endif
