#ifndef Types_enum_h
#define Types_enum_h

/*
 * enum for analysis type.
 */

enum analysis_types {
	Analyze_lepton_jet,
	Analyze_dilepton,
	Analyze_1l2tau,
	Analyze_2lss1tau,
	Analyze_3l
};

/*
 * enum for event selection type.
 */

// enum for event selection region
enum Selection_types {
	Signal_2lss1tau,
	Signal_1l2tau,
	Signal_3l,
	Control_2los1tau,
	Control_1lfakeable,
	Control_WZ
};


#endif
