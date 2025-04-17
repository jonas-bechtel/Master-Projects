#pragma once
namespace Math
{
	double LongTempToVelocitySpread(double kT_long);
	double TransTempToVelocitySpread(double kT_trans);
	double VelocitySpreadToLongTemp(double longVel);
	double VelocitySpreadToTransTemp(double transVel);
}

