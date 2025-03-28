#pragma once
#include <TVector3.h>

struct NumericalIntegrationParameter
{
	double relativeVelocity = 0; // longitudinal
	double ionCharge = 1;
	double electronDensity = 3.1e11;
	double kT_trans = 0.002;
	double kT_long = 187e-6;

	float relativeVelocityRange[2] = { -100000.0f, 100000.0f };
	int numberPoints = 300;

	std::string String() const;
	void FromString(std::string& input);
	void ShowWindow(bool& show);
};

namespace CoolingForceModel
{
	TVector3 CoolingForce(TVector3 relativeVelocity, double kT_trans, double electronDensity, int ionCharge, bool onlyVRelLongInLC);

	double NumericalIntegrand(double* vels, double* params);
	double FlattenedMaxwellDistribution(double vTrans, double vLong, double sigmaTrans, double sigmaLong);
	double CoulombLogarithm(double relativeVelocity, double kT_trans, double electronDensity, int ionCharge);
	double B_min(double relativeVelocity, double kT_trans, int ionCharge);
	double B_max(double relativeVelocity, double kT_trans, double electronDensity);
	double DebyeScreeningLength(double kT_trans, double electronDensity);
	double PlasmaFrequency(double electronDensity);
}
