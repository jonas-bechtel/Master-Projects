#pragma once
#include <TVector3.h>
#include "MathUtils.h"

namespace CoolingForce
{
namespace Model
{
	enum class Type
	{
		NonMagOriginal, JSPEC_Parkhomchuk, JSPEC_NonMagNumeric3D, JSPEC_DerbenovSkrinsky,
		Betacool_Parkhomchuk, Betacool_NonMag, Betacool_NonMagNumeric3D, Betacool_DerbenovSkrinsky, Betacool_Toeppfler, Count
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
		bool magnetic = true;
		bool adiabatic = true;
		bool fast = true;

		inline static float relativeVelocityRange[2] = { -30000.0f, 30000.0f };
		inline static int numberPoints = 100;

		inline static bool showLinesLong = false;
		inline static bool showLinesTrans = false;

		std::string String() const;
		void FromString(std::string& input);
		void ShowWindow(bool& show);
		void ShowValues();
		void ShowVelocityLines();
	};

	// cooling force, direction = 0: longitudinal, 1: transverse
	double Force(const Parameter& params, int direction = 0);

	TVector3 CoolingForce(const Parameter& params);
	std::array<double, 2> Force_Original(const Parameter& params);
	//double ForceTrans_Original(const Parameter& params);

	double NumericalIntegrandPolar(double* vels, double* params);
	//double NumericalIntegrandPolarTrans(double* vels, double* params);
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
		std::array<double,2> Force_Parkhomchuk(const Parameter& parameter);
		std::array<double, 2> Force_NonMagNumeric3D(const Parameter& parameter);
		std::array<double, 2> Force_DerbenovSkrinsky(const Parameter& parameter);
	}

	namespace Betacool
	{
		std::array<double, 2> ForceZ_Parkhomchuk(const Parameter& parameter);
		std::array<double, 2> ForceZ_NonMag(const Parameter& parameter);
		std::array<double, 2> ForceZ_NonMagNumeric3D(const Parameter& parameter);
		std::array<double, 2> ForceZ_DerbenovSkrinsky(const Parameter& parameter);
		std::array<double, 2> ForceZ_Toepffler(const Parameter& parameter);
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
