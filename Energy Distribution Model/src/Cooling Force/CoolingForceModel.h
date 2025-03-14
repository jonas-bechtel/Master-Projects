#pragma once
#include <TVector3.h>

namespace cf
{
	TVector3 CoolingForce(TVector3 relativeVelocity, double kT_trans, double electronDensity, int ionCharge);

	double CoulombLogarithm(double relativeVelocity, double kT_trans, double electronDensity, int ionCharge);
	double B_min(double relativeVelocity, double kT_trans, int ionCharge);
	double B_max(double relativeVelocity, double kT_trans, double electronDensity);
	double DebyeScreeningLength(double kT_trans, double electronDensity);
	double PlasmaFrequency(double electronDensity);
}
