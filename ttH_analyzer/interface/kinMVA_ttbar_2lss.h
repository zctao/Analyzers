#ifndef kinMVA_ttbar_2lss_h
#define kinMVA_ttbar_2lss_h

#include "Analyzers/ttH_analyzer/interface/kinMVA_2lss.h"

class kinMVA_ttbar_2lss : public kinMVA_2lss
{
 public:

	kinMVA_ttbar_2lss();
	~kinMVA_ttbar_2lss();

	void Set_up_Reader(TMVA::Reader*);

	void Calculate_mvaVars(const CU_ttH_EDA_event_vars&, int);

	float Get_avg_dr_jet() {return avg_dr_jet;}
	float Get_MET() {return MET;}

	void Set_avg_dr_jet(const vector<pat::Jet>&);
	void Set_MET(double pfmet) { MET = min(pfmet,400.);}
	
 private:
	float avg_dr_jet;
	float MET;
};

#endif
