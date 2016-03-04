#include "DataFormats/Common/interface/Wrapper.h"
#include "Analyzers/ttH_analyzer/interface/CU_ttH_EDA_Ntuple.h"

namespace { 
  struct dictionary 
  {
    CU_ttH_EDA_Ntuple dummyNtuple1;
    edm::Wrapper<CU_ttH_EDA_Ntuple> dummyNtuple2; 
  };
}
