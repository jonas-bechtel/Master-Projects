#include "pch.h"
#include "CoolingForceModel.h"
#include "Constants.h"

double NumericalIntegrationParameter::precision = 1e-7;

namespace CoolingForceModel
{
	static double constantsFactor =  pow(TMath::Qe(), 3) / (4 * TMath::Pi() * pow(PhysicalConstants::epsilon_0, 2) * PhysicalConstants::electronMass);
}

// returns the 3d force in eV at one position for one relative velocity
TVector3 CoolingForceModel::CoolingForce(const NumericalIntegrationParameter& params)
{
	if (params.electronDensity == 0)
		return TVector3(0, 0, 0);

	double factor = constantsFactor;// pow(TMath::Qe(), 3) / (4 * TMath::Pi() * pow(PhysicalConstants::epsilon_0, 2) * PhysicalConstants::electronMass);
	double L_C = CoulombLogarithm(params.relativeVelocity.Mag(), params.kT_trans, params.electronDensity, params.ionCharge);
	return -factor * pow(params.ionCharge, 2) * params.electronDensity * L_C / pow(params.relativeVelocity.Mag(), 3) * params.relativeVelocity;
}

// integrates the longitudinal cooling force for a given density and temperature and mean relative velocity
double CoolingForceModel::ForceZ(const NumericalIntegrationParameter& params)
{
	if (params.electronDensity == 0)
		return 0.0;

	double deltaTrans = sqrt(2 * params.kT_trans * TMath::Qe() / PhysicalConstants::electronMass);
	double deltaLong = sqrt(params.kT_long * TMath::Qe() / PhysicalConstants::electronMass);

	int numberParams = ceil(sizeof(NumericalIntegrationParameter) / sizeof(double));
	//std::cout << "inside: " << params.String() << std::endl;
	thread_local TF2 func("func", CoolingForceModel::NumericalIntegrandPolar, 0.0, 5.0 * deltaTrans, -5.0 * deltaLong, 5.0 * deltaLong, numberParams, 2);

	func.SetParameters((double*)&params);

	double result = func.Integral(0.0, 5.0 * deltaTrans, -5.0 * deltaLong, 5.0 * deltaLong, params.precision);

	return result;
}

double CoolingForceModel::NumericalIntegrandPolar(double* vels, double* params)
{
	double vTrans = vels[0];
	double vLong = vels[1];

	NumericalIntegrationParameter& parameter = *(NumericalIntegrationParameter*)params;
	double relativeVelocity = parameter.relativeVelocity.z();
	int ionCharge = parameter.ionCharge;
	double electronDensity = parameter.electronDensity;
	double kT_trans = parameter.kT_trans;
	double kT_long = parameter.kT_long;

	/*std::cout << relativeVelocity << std::endl;
	std::cout << ionCharge << std::endl;
	std::cout << electronDensity << std::endl;
	std::cout << kT_trans << std::endl;
	std::cout << kT_long << std::endl;*/

	double factor = ionCharge * ionCharge * electronDensity * constantsFactor;
	double v_r = sqrt(pow(relativeVelocity - vLong, 2) + pow(vTrans, 2));
	double L_C = CoulombLogarithm(v_r, kT_trans, electronDensity, ionCharge);
	double deltaTrans = sqrt(2 * kT_trans * TMath::Qe() / PhysicalConstants::electronMass);
	double deltaLong = sqrt(kT_long * TMath::Qe() / PhysicalConstants::electronMass);
	double f_v = FlattenedMaxwellDistributionPolar(vTrans, vLong, deltaTrans, deltaLong);

	//std::cout << factor << std::endl;
	//std::cout << L_C << std::endl;
	//std::cout << deltaTrans << std::endl;
	//std::cout << deltaLong << std::endl;
	//std::cout << f_v << std::endl;
	double result = -factor * 2 * TMath::Pi() * L_C * f_v * (relativeVelocity - vLong) * vTrans /
		pow(sqrt(pow(relativeVelocity - vLong, 2) + vTrans * vTrans), 3);

	//std::cout << result << std::endl;

	return result;
}

double CoolingForceModel::NumericalIntegrandCartesian(double* vels, double* params)
{
	double vX = vels[0];
	double vY = vels[1];
	double vZ = vels[2];

	NumericalIntegrationParameter& parameter = *(NumericalIntegrationParameter*)params;
	TVector3 meanRelativeVelocity = parameter.relativeVelocity;
	int ionCharge = parameter.ionCharge;
	double electronDensity = parameter.electronDensity;
	double kT_trans = parameter.kT_trans;
	double kT_long = parameter.kT_long;

	double factor = ionCharge * ionCharge * electronDensity * constantsFactor;
	TVector3 relativeVelocity = meanRelativeVelocity - TVector3(vX, vY, vZ);
	//relativeVelocity.Print();
	//meanRelativeVelocity.Unit().Print();
	//std::cout << relativeVelocity.Dot(meanRelativeVelocity.Unit()) << std::endl << std::endl;
	double v_r = relativeVelocity.Mag();
	//double v_r = sqrt(pow(meanRelativeVelocity.z() - vZ, 2) + vY * vY + vX * vX);
	double L_C = CoulombLogarithm(v_r, kT_trans, electronDensity, ionCharge);
	double sigmaX = sqrt(kT_trans * TMath::Qe() / PhysicalConstants::electronMass);
	double sigmaY = sigmaX;
	double sigmaZ = sqrt(kT_long * TMath::Qe() / PhysicalConstants::electronMass);
	double f_v = FlattenedMaxwellDistributionCartesian(vX, vY, vZ, sigmaX, sigmaY, sigmaZ);

	//std::cout << factor << std::endl;
	/*std::cout << v_r << std::endl;
	std::cout << L_C << std::endl;
	std::cout << sigmaX << std::endl;
	std::cout << sigmaZ << std::endl;
	std::cout << f_v << std::endl;*/
	if (relativeVelocity.z() == 0) 
		return 0;

	double result = -factor * L_C * f_v * relativeVelocity.z() / pow(v_r, 3);
	if (std::isnan(result))
	{
		std::cout << "result is nan" << std::endl;
	}
	//std::cout << result << std::endl;

	return result;
}

double CoolingForceModel::FlattenedMaxwellDistributionPolar(double vTrans, double vLong, double deltaTrans, double deltaLong)
{
	return 1 / (pow(deltaTrans, 2) * sqrt(2) * deltaLong * pow(TMath::Pi(), 1.5)) *
		exp(-((pow(vTrans, 2) / pow(deltaTrans, 2))
			+ (pow(vLong, 2) / (2 * pow(deltaLong, 2)))));
}

double CoolingForceModel::FlattenedMaxwellDistributionCartesian(double vX, double vY, double vZ, double sigmaX, double sigmaY, double sigmaZ)
{
	return 1 / (sigmaX * sigmaY * sigmaZ * pow(2 * TMath::Pi(), 1.5)) *
		exp(-0.5 * (pow(vX, 2) / pow(sigmaX, 2)
				 +  pow(vY, 2) / pow(sigmaY, 2)
				 +  pow(vZ, 2) / pow(sigmaZ, 2)));
}

double CoolingForceModel::CoulombLogarithm(double relativeVelocity, double kT_trans, double electronDensity, int ionCharge)
{
	return log(B_max(abs(relativeVelocity), kT_trans, electronDensity) / B_min(abs(relativeVelocity), kT_trans, ionCharge));
}

double CoolingForceModel::B_min(double relativeVelocity, double kT_trans, int ionCharge)
{
	double deltaE_trans = sqrt(2 * kT_trans * TMath::Qe() / PhysicalConstants::electronMass);
	return ionCharge * TMath::Qe() * TMath::Qe() / (4 * TMath::Pi() * PhysicalConstants::epsilon_0 * PhysicalConstants::electronMass * pow(std::max(relativeVelocity, deltaE_trans), 2));
}

double CoolingForceModel::B_max(double relativeVelocity, double kT_trans, double electronDensity)
{
	return std::max(relativeVelocity / PlasmaFrequency(electronDensity), DebyeScreeningLength(kT_trans, electronDensity));
}

double CoolingForceModel::DebyeScreeningLength(double kT_trans, double electronDensity)
{
	return sqrt(PhysicalConstants::epsilon_0 * kT_trans * TMath::Qe() / (electronDensity * TMath::Qe() * TMath::Qe()));
}

double CoolingForceModel::PlasmaFrequency(double electronDensity)
{
	return sqrt(electronDensity * TMath::Qe() * TMath::Qe() / (PhysicalConstants::epsilon_0 * PhysicalConstants::electronMass));
}

double CoolingForceModel::GetConstantsFactor()
{
	return 0.0;
}

std::string NumericalIntegrationParameter::String() const
{
	std::ostringstream density;
	density << std::scientific << std::setprecision(3) << electronDensity;
	std::ostringstream Tlong;
	Tlong << std::scientific << std::setprecision(3) << kT_long;
	std::ostringstream Ttrans;
	Ttrans << std::scientific << std::setprecision(3) << kT_trans;

	std::string result = "Q=" + std::to_string(ionCharge) + "_ne=" + density.str()
		+ "_Tlong=" + Tlong.str() + "_Ttrans=" + Ttrans.str();

	return result;
}

void NumericalIntegrationParameter::FromString(std::string& input)
{
	std::istringstream stream(input);
	std::string token;

	while (std::getline(stream, token, '_')) 
	{
		if (token.find("Q=") == 0) 
		{
			ionCharge = std::stoi(token.substr(2));
		}
		else if (token.find("ne=") == 0) 
		{
			electronDensity = std::stod(token.substr(3));
		}
		else if (token.find("Tlong=") == 0) 
		{
			kT_long = std::stod(token.substr(6)); 
		}
		else if (token.find("Ttrans=") == 0) 
		{
			kT_trans = std::stod(token.substr(7));
		}
	}
}

void NumericalIntegrationParameter::ShowWindow(bool& show)
{
	if (!show)
	{
		return;
	}
	if (ImGui::Begin("Numerical integration parameter", &show, ImGuiWindowFlags_NoDocking))
	{
		ImGui::PushItemWidth(150.0f);
		ImGui::InputFloat2("velocity range", relativeVelocityRange, "%.1e");
		ImGui::InputInt("number points", &numberPoints);
		ImGui::InputDouble("relative precision", &precision, 0, 0, "%.1e");
		ImGui::Separator();
		ImGui::InputDouble("kT long", &kT_long);
		ImGui::InputDouble("kT trans", &kT_trans);
		ImGui::InputInt("ion charge", &ionCharge);
		ImGui::InputDouble("electron density", &electronDensity, 0, 0, "%.2e");
		ImGui::PopItemWidth();
	}
	ImGui::End();
}
