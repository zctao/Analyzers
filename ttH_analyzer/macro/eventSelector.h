#ifndef passTauCharge_h
#define passTauCharge_h

#include "TString.h"

#include <iostream>

bool passTauCharge(int lepCharge, int tauCharge, TString region="signal")
{
	if (region=="signal") {
		return lepCharge + tauCharge == 0;
	}
	else if (region=="control") {
		return lepCharge == tauCharge;
	}
	else {
		std::cout << "WARNING: Selection region not defined." << std::endl;
	}

	return false;
}

#endif
