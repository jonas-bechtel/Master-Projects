#pragma once
#include <TVector3.h>
#include "Math.h"

namespace CoolingForce
{
namespace Model
{
	enum class Type
	{
		NonMagOriginal, Parkhomchuk, NonMagNumeric3D, DerbenovSkrinsky, Count
	};

	struct Parameter
	{
		Type model;

		TVector3 relativeVelocity = { 0,0,0 };
		int ionCharge = 1;
		double electronDensity = 3.1e11;
		double kT_trans = 0.002;
		double kT_long = 187e-6;
		double velSpreadTrans = Math::TransTempToVelocitySpread(kT_trans);
		double velSpreadLong = Math::LongTempToVelocitySpread(kT_long);

		double coolerTime = 1e-6;
		double magneticField = 0.01; // T
		double effectiveVelocity = 1;
		double relTolerance = 1e-8;

		bool useGSL = true;

		// Derbenov Skrinsky Model
		double smoothingFactor = 2;
		bool magneticOnly = false;

		static float relativeVelocityRange[2];
		static int numberPoints;

		static bool showLinesLong;
		static bool showLinesTrans;

		std::string String() const;
		void FromString(std::string& input);
		void ShowWindow(bool& show);
		void ShowValues();
		void ShowVelocityLines();
	};

	TVector3 CoolingForce(const Parameter& params);
	double ForceZ(const Parameter& params);

	double NumericalIntegrandPolar(double* vels, double* params);
	double NumericalIntegrandCartesian(double* vels, double* params);
	double FlattenedMaxwellDistributionPolar(double vTrans, double vLong, double deltaTrans, double deltaLong);
	double FlattenedMaxwellDistributionCartesian(double vX, double vY, double vZ, double sigmaX, double sigmaY, double sigmaZ);
	double CoulombLogarithm(double relativeVelocity, double kT_trans, double electronDensity, int ionCharge);
	double B_min(double relativeVelocity, double kT_trans, int ionCharge);
	double B_max(double relativeVelocity, double kT_trans, double electronDensity);
	double TransverseDebyeScreeningLength(double kT_trans, double electronDensity);
	double LongitudinalDebyeScreeningLength(double kT_long, double electronDensity);
	double CyclotronFrequency(double magneticField);
	double PlasmaFrequency(double electronDensity);
	double CyclotronRadius(double magneticField, double transverseVelocity);

	namespace JSPEC
	{
		double ForceZ_Parkhomchuk(const Parameter& parameter);
		double ForceZ_NonMagNumeric3D(const Parameter& parameter);
		double ForceZ_DerbenovSkrinsky(const Parameter& parameter);
	}

	/*namespace Adiabatic
	{
		double L_C_nm(double u_ad, double kT_trans, int ionCharge, double magneticField);
		double L_C_ad(double u_ad, double kT_long, int ionCharge, double magneticField, double electronDensity);
		double B_min_nm(double u_ad, double kT_trans, int ionCharge);
		double B_max_nm(double u_ad, double magneticField);
		double B_min_ad(double u_ad, double magneticField, int ionCharge, double kT_long);
		double B_max_ad(double u_ad, double electronDensity, double kT_long);
	}*/
}
}
