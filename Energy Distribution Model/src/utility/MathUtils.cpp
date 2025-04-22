#include "pch.h"
#include "MathUtils.h"
#include "Constants.h"

namespace Math
{
	double LongTempToVelocitySpread(double kT_long)
	{
		return TMath::Sqrt(kT_long * TMath::Qe() / PhysicalConstants::electronMass);
	}

	double TransTempToVelocitySpread(double kT_trans)
	{
		return TMath::Sqrt(2 * kT_trans * TMath::Qe() / PhysicalConstants::electronMass);
	}

	double VelocitySpreadToLongTemp(double longVel)
	{
		return longVel * longVel * PhysicalConstants::electronMass / TMath::Qe();
	}

	double VelocitySpreadToTransTemp(double transVel)
	{
		return transVel * transVel * PhysicalConstants::electronMass / TMath::Qe() / 2;
	}
}

