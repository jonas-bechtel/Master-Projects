#include "pch.h"
#include "CoolingForceModel.h"
#include "Constants.h"

// returns the 3d force in eV
TVector3 cf::CoolingForce(TVector3 relativeVelocity, double kT_trans, double electronDensity, int ionCharge)
{
	if (electronDensity == 0)
		return TVector3(0, 0, 0);

	double factor = pow(TMath::Qe(), 3) / (4 * TMath::Pi() * pow(PhysicalConstants::epsilon_0, 2) * PhysicalConstants::electronMass);
	double L_C = CoulombLogarithm(relativeVelocity.Mag(), kT_trans, electronDensity, ionCharge);
	return -factor * ionCharge * ionCharge * electronDensity * L_C / pow(relativeVelocity.Mag(), 3) * relativeVelocity;
}

double cf::CoulombLogarithm(double relativeVelocity, double kT_trans, double electronDensity, int ionCharge)
{
	return log(B_max(relativeVelocity, kT_trans, electronDensity) / B_min(relativeVelocity, kT_trans, ionCharge));
}

double cf::B_min(double relativeVelocity, double kT_trans, int ionCharge)
{
	double deltaE_trans = sqrt(2 * kT_trans * TMath::Qe() / PhysicalConstants::electronMass);
	return ionCharge * TMath::Qe() * TMath::Qe() / (4 * TMath::Pi() * PhysicalConstants::epsilon_0 * PhysicalConstants::electronMass * pow(std::max(relativeVelocity, deltaE_trans), 2));
}

double cf::B_max(double relativeVelocity, double kT_trans, double electronDensity)
{
	return std::max(relativeVelocity / PlasmaFrequency(electronDensity), DebyeScreeningLength(kT_trans, electronDensity));
}

double cf::DebyeScreeningLength(double kT_trans, double electronDensity)
{
	return sqrt(PhysicalConstants::epsilon_0 * kT_trans * TMath::Qe() / (electronDensity * TMath::Qe() * TMath::Qe()));
}

double cf::PlasmaFrequency(double electronDensity)
{
	return sqrt(electronDensity * TMath::Qe() * TMath::Qe() / (PhysicalConstants::epsilon_0 * PhysicalConstants::electronMass));
}
