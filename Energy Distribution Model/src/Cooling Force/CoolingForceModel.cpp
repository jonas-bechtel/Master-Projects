#include "pch.h"
#include "CoolingForceModel.h"
#include "Constants.h"

// returns the 3d force in eV
TVector3 CoolingForceModel::CoolingForce(TVector3 relativeVelocity, double kT_trans, double electronDensity, int ionCharge)
{
	if (electronDensity == 0)
		return TVector3(0, 0, 0);

	double factor = pow(TMath::Qe(), 3) / (4 * TMath::Pi() * pow(PhysicalConstants::epsilon_0, 2) * PhysicalConstants::electronMass);
	double L_C = CoulombLogarithm(relativeVelocity.Mag(), kT_trans, electronDensity, ionCharge);
	return -factor * ionCharge * ionCharge * electronDensity * L_C / pow(relativeVelocity.Mag(), 3) * relativeVelocity;
}

double CoolingForceModel::NumericalIntegrand(double* vels, double* params)
{
	double vTrans = vels[0];
	double vLong = vels[1];

	double relativeVelocity = params[0];
	double ionCharge = params[1];
	double electronDensity = params[2];
	double kT_trans = params[3];
	double kT_long = params[4];

	//std::cout << relativeVelocity << std::endl;
	//std::cout << ionCharge << std::endl;
	//std::cout << electronDensity << std::endl;
	//std::cout << kT_trans << std::endl;
	//std::cout << kT_long << std::endl;

	double factor = ionCharge * ionCharge * pow(TMath::Qe(), 3) * electronDensity 
		/ (4 * TMath::Pi() * pow(PhysicalConstants::epsilon_0, 2)  * PhysicalConstants::electronMass);
	double L_C = CoulombLogarithm(vLong, kT_trans, electronDensity, ionCharge);
	double deltaTrans = sqrt(2 * kT_trans * TMath::Qe() / PhysicalConstants::electronMass);
	double deltaLong = sqrt(kT_long * TMath::Qe() / PhysicalConstants::electronMass);
	double f_v = FlattenedMaxwellDistribution(vTrans, vLong, deltaTrans, deltaLong);

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

double CoolingForceModel::FlattenedMaxwellDistribution(double vTrans, double vLong, double deltaTrans, double deltaLong)
{
	return 1 / (pow(deltaTrans, 2) * sqrt(2) * deltaLong * pow(TMath::Pi(), 1.5)) *
		exp(-((pow(vTrans, 2) / pow(deltaTrans, 2))
			+ (pow(vLong, 2) / (2 * pow(deltaLong, 2)))));
}

double CoolingForceModel::CoulombLogarithm(double relativeVelocity, double kT_trans, double electronDensity, int ionCharge)
{
	return log(B_max(relativeVelocity, kT_trans, electronDensity) / B_min(relativeVelocity, kT_trans, ionCharge));
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
		ImGui::InputDouble("kT long", &kT_long);
		ImGui::InputDouble("kT trans", &kT_trans);
		ImGui::InputDouble("ion charge", &ionCharge);
		ImGui::InputDouble("electron density", &electronDensity, 0, 0, "%.2e");
		ImGui::PopItemWidth();
	}
	ImGui::End();
}
