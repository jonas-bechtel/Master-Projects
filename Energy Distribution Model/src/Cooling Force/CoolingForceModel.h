#pragma once
#include <TVector3.h>

struct NumericalIntegrationParameter
{
	// mean relative velocity, no electron Temperature effects
	TVector3 relativeVelocity = { 0, 0, 1 };
	//double relativeVelocity = 0; // longitudinal
	int ionCharge = 1;
	double electronDensity = 3.1e11;
	double kT_trans = 0.002;
	double kT_long = 187e-6;

	float relativeVelocityRange[2] = { -100000.0f, 100000.0f };
	int numberPoints = 50;
	static double precision;

	std::string String() const;
	void FromString(std::string& input);
	void ShowWindow(bool& show);
};

namespace CoolingForceModel
{
	TVector3 CoolingForce(const NumericalIntegrationParameter& params);
	double ForceZ(const NumericalIntegrationParameter& params);

	double NumericalIntegrandPolar(double* vels, double* params);
	double NumericalIntegrandCartesian(double* vels, double* params);
	double FlattenedMaxwellDistributionPolar(double vTrans, double vLong, double deltaTrans, double deltaLong);
	double FlattenedMaxwellDistributionCartesian(double vX, double vY, double vZ, double sigmaX, double sigmaY, double sigmaZ);
	double CoulombLogarithm(double relativeVelocity, double kT_trans, double electronDensity, int ionCharge);
	double B_min(double relativeVelocity, double kT_trans, int ionCharge);
	double B_max(double relativeVelocity, double kT_trans, double electronDensity);
	double DebyeScreeningLength(double kT_trans, double electronDensity);
	double PlasmaFrequency(double electronDensity);
	double GetConstantsFactor();
}
