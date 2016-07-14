#ifndef kinMVA_ttV_2lss_h
#define kinMVA_ttV_2lss_h

#include "Analyzers/ttH_analyzer/interface/kinMVA_2lss.h"

class kinMVA_ttV_2lss : public kinMVA_2lss
{
 public:

	kinMVA_ttV_2lss();
	~kinMVA_ttV_2lss();

	void Set_up_Reader(TMVA::Reader*);

	void Calculate_mvaVars(const CU_ttH_EDA_event_vars&, int);
};

#endif
