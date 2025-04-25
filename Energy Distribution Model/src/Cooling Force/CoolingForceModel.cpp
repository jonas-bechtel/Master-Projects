#include "pch.h"
#include "CoolingForceModel.h"
#include "Constants.h"
#include "MathUtils.h"

// JSPEC force functions
#include "force.h"
#include "beam.h"

// Betacool force functions
#include "xForce.h"


namespace CoolingForce
{
namespace Model
{
	static double constantsFactor = pow(TMath::Qe(), 3) / (4 * TMath::Pi() * pow(PhysicalConstants::epsilon_0, 2) * PhysicalConstants::electronMass);
	static const char* modelTypes[] = { "NonMagOriginal", "JSPEC_Parkhomchuk", "JSPEC_NonMagNumeric3D", "JSPEC_DerbenovSkrinsky",
		"Betacool_Parkhomchuk", "Betacool_NonMag", "Betacool_NonMagNumeric3D", "Betacool_DerbenovSkrinsky", "Betacool_Toeppfler" };
	static int currentModelIndex = 0;
	float Parameter::relativeVelocityRange[2] = { -100000.0f, 100000.0f };
	int Parameter::numberPoints = 100;
	bool Parameter::showLinesLong = false;
	bool Parameter::showLinesTrans = false;

	std::string Parameter::String() const
	{
		std::ostringstream density;
		density << std::scientific << std::setprecision(3) << electronDensity;
		std::ostringstream Tlong;
		Tlong << std::scientific << std::setprecision(3) << kT_long;
		std::ostringstream Ttrans;
		Ttrans << std::scientific << std::setprecision(3) << kT_trans;
		std::ostringstream Bfield;
		Bfield << std::setprecision(3) << magneticField;
		std::ostringstream effVel;
		effVel << std::scientific << std::setprecision(1) << effectiveVelocity;

		std::string result = std::string(modelTypes[(int)model]) + "_Q=" + std::to_string(ionCharge) + "_ne=" + density.str()
			+ "_Tlong=" + Tlong.str() + "_Ttrans=" + Ttrans.str();
		if (model == Type::JSPEC_Parkhomchuk || model == Type::Betacool_Parkhomchuk)
		{
			result += "_B=" + Bfield.str() + "_Veff=" + effVel.str();
		}
		if (model == Type::JSPEC_DerbenovSkrinsky || model == Type::Betacool_DerbenovSkrinsky || model == Type::Betacool_Toeppfler)
		{
			result += "_B=" + Bfield.str();
		}
		return result;
	}

	void Parameter::FromString(std::string& input)
	{
		std::istringstream stream(input);
		std::string token;

		// get first elemetn which is the model type
		std::getline(stream, token, '_');
		for (int i = 0 ; i < (int)Type::Count; i++)
		{
			if (std::string(modelTypes[i]) == token)
			{
				model = (Type)i;
			}
		}
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
			else if (token.find("B=") == 0)
			{
				magneticField = std::stod(token.substr(2));
			}
			else if (token.find("Veff=") == 0)
			{
				effectiveVelocity = std::stod(token.substr(5));
			}
		}
	}

	void Parameter::ShowWindow(bool& show)
	{
		if (!show)
		{
			return;
		}
		if (ImGui::Begin("Numerical integration parameter", &show, ImGuiWindowFlags_NoDocking))
		{
			ImGui::PushItemWidth(140.0f);
			ImGui::InputFloat2("velocity range", relativeVelocityRange, "%.1e");
			ImGui::InputInt("number points", &numberPoints);
			ImGui::Separator();
			if (ImGui::Combo("models", &currentModelIndex, modelTypes, IM_ARRAYSIZE(modelTypes)))
			{
				model = static_cast<Model::Type>(currentModelIndex);
			}
			ImGui::InputInt("ion charge", &ionCharge);
			ImGui::InputDouble("electron density", &electronDensity, 0, 0, "%.2e");
			ImGui::BeginGroup();
			if (ImGui::InputDouble("kT long", &kT_long))
			{
				velSpreadLong = Math::LongTempToVelocitySpread(kT_long);

			}
			if (ImGui::InputDouble("kT trans", &kT_trans))
			{
				velSpreadTrans = Math::TransTempToVelocitySpread(kT_trans);
			}
			ImGui::EndGroup();
			ImGui::SameLine();
			ImGui::BeginGroup();
			if (ImGui::InputDouble("vel long", &velSpreadLong, 0, 0, "%.0f"))
			{
				kT_long = Math::VelocitySpreadToLongTemp(velSpreadLong);

			}
			if (ImGui::InputDouble("vel trans", &velSpreadTrans, 0, 0, "%.0f"))
			{
				kT_trans = Math::VelocitySpreadToTransTemp(velSpreadTrans);
			}
			ImGui::EndGroup();
			ImGui::SameLine();
			ImGui::BeginGroup();
			ImGui::Checkbox("show##long", &showLinesLong);
			ImGui::Checkbox("show##trans", &showLinesTrans);
			ImGui::EndGroup();

			if (model == Model::Type::JSPEC_Parkhomchuk || model == Model::Type::Betacool_Parkhomchuk)
			{
				ImGui::InputDouble("effective vel trans", &effectiveVelocity);
				ImGui::InputDouble("magnetic field", &magneticField);
			}
			if (model == Model::Type::NonMagOriginal || model == Model::Type::JSPEC_NonMagNumeric3D)
			{
				ImGui::InputDouble("relative tolerance", &relTolerance, 0, 0, "%.1e");
			}
			if (model == Model::Type::JSPEC_DerbenovSkrinsky || model == Model::Type::Betacool_DerbenovSkrinsky)
			{
				ImGui::InputDouble("magnetic field", &magneticField);
				ImGui::InputDouble("smoothing factor", &smoothingFactor);
				ImGui::Checkbox("magnetic", &magnetic);
				ImGui::Checkbox("adiabatic", &adiabatic);
				ImGui::Checkbox("fast", &fast);
			}
			if (model == Model::Type::Betacool_Toeppfler)
			{
				ImGui::InputDouble("magnetic field", &magneticField);
			}

			ImGui::PopItemWidth();
		}
		ImGui::End();
	}

	void Parameter::ShowValues()
	{
		ImGui::BeginGroup();
		ImGui::Text("model:");
		ImGui::Text("electron density:");
		ImGui::Text("ion charge:");
		ImGui::Text("longitudinal kT:");
		ImGui::Text("transverse kT:");

		if (model == Model::Type::JSPEC_Parkhomchuk || model == Model::Type::Betacool_Parkhomchuk)
		{
			ImGui::Text("effective vel trans:");
			ImGui::Text("magnetic field:");
		}
		if (model == Model::Type::NonMagOriginal || model == Model::Type::JSPEC_NonMagNumeric3D)
		{
			ImGui::Text("relative tolerance");
		}
		if (model == Model::Type::JSPEC_DerbenovSkrinsky || model == Model::Type::Betacool_DerbenovSkrinsky)
		{
			ImGui::Text("magnetic field:");
			ImGui::Text("smoothing factor:");
			ImGui::Text("magnetic:");
			ImGui::Text("adiabatic:");
			ImGui::Text("fast:");
		}
		if (model == Model::Type::Betacool_Toeppfler)
		{
			ImGui::Text("magnetic field:");
		}
		
		ImGui::EndGroup();

		ImGui::SameLine();
		ImGui::BeginGroup();
		ImGui::Text("%s", modelTypes[(int)model]);
		ImGui::Text("%.2e [1/m^3]", electronDensity);
		ImGui::Text("%d", ionCharge);
		ImGui::Text("%.2e [eV]", kT_long);
		ImGui::Text("%.2e [eV]", kT_trans);

		if (model == Model::Type::JSPEC_Parkhomchuk || model == Model::Type::Betacool_Parkhomchuk)
		{
			ImGui::Text("%.0f [m/s]", effectiveVelocity);
			ImGui::Text("%.3f [T]", magneticField);
		}
		if (model == Model::Type::NonMagOriginal || model == Model::Type::JSPEC_NonMagNumeric3D)
		{
			ImGui::Text("%.1e", relTolerance);
		}
		if (model == Model::Type::JSPEC_DerbenovSkrinsky || model == Model::Type::Betacool_DerbenovSkrinsky)
		{
			ImGui::Text("%.3f [T]", magneticField);
			ImGui::Text("%.1f", smoothingFactor);
			ImGui::Text("%s", magnetic ? "True" : "False");
			ImGui::Text("%s", adiabatic ? "True" : "False");
			ImGui::Text("%s", fast ? "True" : "False");
		}
		if (model == Model::Type::Betacool_Toeppfler)
		{
			ImGui::Text("%.3f [T]", magneticField);
		}
		ImGui::EndGroup();
	}

	void Parameter::ShowVelocityLines()
	{
		double negVelLongTemp = -velSpreadLong;
		double negVelTransTemp = -velSpreadTrans;
		bool changed1 = false;
		bool changed2 = false;
		ImVec4 colorLong = ImVec4(0.6f, 0.05f, 0.18f, 1.0f);
		ImVec4 colorTrans = ImVec4(0.8f, 0.31f, 0.22f, 1.0f);

		if (showLinesLong)
		{
			changed1 |= ImPlot::DragLineX(0, &velSpreadLong, colorLong, 1, ImPlotDragToolFlags_NoFit);
			ImPlot::TagX(velSpreadLong, colorLong, "kT long");
			changed2 |= ImPlot::DragLineX(1, &negVelLongTemp, colorLong, 1, ImPlotDragToolFlags_NoFit);
			ImPlot::TagX(negVelLongTemp, colorLong, "kT long");

		}
		if (showLinesTrans)
		{
			changed1 |= ImPlot::DragLineX(2, &velSpreadTrans, colorTrans, 1, ImPlotDragToolFlags_NoFit);
			ImPlot::TagX(velSpreadTrans, colorTrans, "kT trans");
			changed2 |= ImPlot::DragLineX(3, &negVelTransTemp, colorTrans, 1, ImPlotDragToolFlags_NoFit);
			ImPlot::TagX(negVelTransTemp, colorTrans, "kT trans");
		}
		
		if (changed2)
		{
			velSpreadLong = -negVelLongTemp;
			velSpreadTrans = -negVelTransTemp;
		}
		if (changed1 || changed2)
		{
			kT_long = Math::VelocitySpreadToLongTemp(velSpreadLong);
			kT_trans = Math::VelocitySpreadToTransTemp(velSpreadTrans);
		}
	}

	double ForceZ(const Parameter& params)
	{
		if (params.electronDensity == 0)
			return 0.0;

		double result = 0;
		switch (params.model)
		{
		case Model::Type::NonMagOriginal:
			result = Model::ForceZ_Original(params);
			break;
		case Model::Type::JSPEC_Parkhomchuk:
			result = Model::JSPEC::ForceZ_Parkhomchuk(params);
			break;
		case Model::Type::JSPEC_NonMagNumeric3D:
			result = Model::JSPEC::ForceZ_NonMagNumeric3D(params);
			break;
		case Model::Type::JSPEC_DerbenovSkrinsky:
			result = Model::JSPEC::ForceZ_DerbenovSkrinsky(params);
			break;
		case Model::Type::Betacool_NonMag:
			result = Model::Betacool::ForceZ_NonMag(params);
			break;
		case Model::Type::Betacool_DerbenovSkrinsky:
			result = Model::Betacool::ForceZ_DerbenovSkrinsky(params);
			break;
		case Model::Type::Betacool_Toeppfler:
			result = Model::Betacool::ForceZ_Toepffler(params);
			break;
		case Model::Type::Betacool_Parkhomchuk:
			result = Model::Betacool::ForceZ_Parkhomchuk(params);
			break;
		case Model::Type::Betacool_NonMagNumeric3D:
			result = Model::Betacool::ForceZ_NonMagNumeric3D(params);
			break;
		}
		
		return result;
	}


	// returns the 3d force in eV at one position for one relative velocity
	TVector3 CoolingForce(const Parameter& params)
	{
		if (params.electronDensity == 0)
			return TVector3(0, 0, 0);

		double factor = constantsFactor;// pow(TMath::Qe(), 3) / (4 * TMath::Pi() * pow(PhysicalConstants::epsilon_0, 2) * PhysicalConstants::electronMass);
		double L_C = CoulombLogarithm(params.relativeVelocity.Mag(), params.kT_trans, params.electronDensity, params.ionCharge);
		return -factor * pow(params.ionCharge, 2) * params.electronDensity * L_C / pow(params.relativeVelocity.Mag(), 3) * params.relativeVelocity;
	}

	// integrates the longitudinal cooling force for a given density and temperature and mean relative velocity
	double ForceZ_Original(const Parameter& params)
	{
		if (params.electronDensity == 0)
			return 0.0;

		double deltaTrans = Math::TransTempToVelocitySpread(params.kT_trans);
		double deltaLong = Math::LongTempToVelocitySpread(params.kT_long);

		int numberParams = ceil(sizeof(Model::Parameter) / sizeof(double));
		//std::cout << "inside: " << params.String() << std::endl;
		thread_local TF2 func("func", Model::NumericalIntegrandPolar, 0.0, 5.0 * deltaTrans, -5.0 * deltaLong, 5.0 * deltaLong, numberParams, 2);

		func.SetParameters((double*)&params);

		double result = func.Integral(0.0, 5.0 * deltaTrans, -5.0 * deltaLong, 5.0 * deltaLong, params.relTolerance);

		return result;
	}

	double NumericalIntegrandPolar(double* vels, double* params)
	{
		double vTrans = vels[0];
		double vLong = vels[1];

		Parameter& parameter = *(Parameter*)params;
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
		double L_c = CoulombLogarithm(v_r, kT_trans, electronDensity, ionCharge);
		double deltaTrans = Math::TransTempToVelocitySpread(kT_trans);
		double deltaLong = Math::LongTempToVelocitySpread(kT_long);
		double f_v = FlattenedMaxwellDistributionPolar(vTrans, vLong, deltaTrans, deltaLong);

		//std::cout << factor << std::endl;
		//std::cout << L_C << std::endl;
		//std::cout << deltaTrans << std::endl;
		//std::cout << deltaLong << std::endl;
		//std::cout << f_v << std::endl;
		double result = -factor * 2 * TMath::Pi() * L_c * f_v * (relativeVelocity - vLong) * vTrans /
			pow(sqrt(pow(relativeVelocity - vLong, 2) + vTrans * vTrans), 3);

		if (ionCharge < 0)
			result += -factor * 4 * TMath::Pi() * f_v * (relativeVelocity - vLong) * vTrans /
			pow(sqrt(pow(relativeVelocity - vLong, 2) + vTrans * vTrans), 3);

		//std::cout << result << std::endl;

		return result;
	}

	double NumericalIntegrandCartesian(double* vels, double* params)
	{
		double vX = vels[0];
		double vY = vels[1];
		double vZ = vels[2];

		Parameter& parameter = *(Parameter*)params;
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

	double FlattenedMaxwellDistributionPolar(double vTrans, double vLong, double deltaTrans, double deltaLong)
	{
		return 1 / (pow(deltaTrans, 2) * sqrt(2) * deltaLong * pow(TMath::Pi(), 1.5)) *
			exp(-((pow(vTrans, 2) / pow(deltaTrans, 2))
				+ (pow(vLong, 2) / (2 * pow(deltaLong, 2)))));
	}

	double FlattenedMaxwellDistributionCartesian(double vX, double vY, double vZ, double sigmaX, double sigmaY, double sigmaZ)
	{
		return 1 / (sigmaX * sigmaY * sigmaZ * pow(2 * TMath::Pi(), 1.5)) *
			exp(-0.5 * (pow(vX, 2) / pow(sigmaX, 2)
				+ pow(vY, 2) / pow(sigmaY, 2)
				+ pow(vZ, 2) / pow(sigmaZ, 2)));
	}

	double CoulombLogarithm(double relativeVelocity, double kT_trans, double electronDensity, int ionCharge)
	{
		double bMax = B_max(abs(relativeVelocity), kT_trans, electronDensity);
		double bMin = B_min(abs(relativeVelocity), kT_trans, ionCharge);
		if (ionCharge < 0)
			bMin *= 2;
		return log(bMax / bMin);
	}

	double B_min(double relativeVelocity, double kT_trans, int ionCharge)
	{
		double deltaE_trans = sqrt(2 * kT_trans * TMath::Qe() / PhysicalConstants::electronMass);
		return abs(ionCharge) * TMath::Qe() * TMath::Qe() / (4 * TMath::Pi() * PhysicalConstants::epsilon_0 * PhysicalConstants::electronMass * pow(std::max(relativeVelocity, deltaE_trans), 2));
	}

	double B_max(double relativeVelocity, double kT_trans, double electronDensity)
	{
		return std::max(relativeVelocity / PlasmaFrequency(electronDensity), TransverseDebyeScreeningLength(kT_trans, electronDensity));
	}

	/*namespace Adiabatic
	{
		double L_C_nm(double u_ad, double kT_trans, int ionCharge, double magneticField)
		{
			return 0.0;
		}
		double L_C_ad(double u_ad, double kT_long, int ionCharge, double magneticField, double electronDensity);

		double B_min_nm(double u_ad, double kT_trans, int ionCharge)
		{
			double deltaE_trans = sqrt(2 * kT_trans * TMath::Qe() / PhysicalConstants::electronMass);
			return abs(ionCharge) * TMath::Qe() * TMath::Qe() / (4 * TMath::Pi() * PhysicalConstants::epsilon_0 * PhysicalConstants::electronMass * pow(std::max(u_ad, deltaE_trans), 2));
		}

		double B_max_nm(double u_ad, double magneticField)
		{
			return u_ad / CyclotronFrequency(magneticField);
		}

		double B_min_ad(double u_ad, double magneticField, int ionCharge, double kT_long)
		{
			double first = u_ad / CyclotronFrequency(magneticField);
			double deltaE_long = sqrt(2 * kT_long * TMath::Qe() / PhysicalConstants::electronMass);
			double second = abs(ionCharge) * TMath::Qe() * TMath::Qe() / (4 * TMath::Pi() * PhysicalConstants::epsilon_0 * PhysicalConstants::electronMass * pow(std::max(u_ad, deltaE_long), 2));

			return std::max(first, second);
		}

		double B_max_ad(double u_ad, double electronDensity, double kT_long)
		{
			return std::max(u_ad / PlasmaFrequency(electronDensity), LongitudinalDebyeScreeningLength(kT_long, electronDensity));
		}
	}*/

	double TransverseDebyeScreeningLength(double kT_trans, double electronDensity)
	{
		return sqrt(PhysicalConstants::epsilon_0 * kT_trans * TMath::Qe() / (electronDensity * TMath::Qe() * TMath::Qe()));
	}

	double LongitudinalDebyeScreeningLength(double kT_long, double electronDensity)
	{
		return sqrt(PhysicalConstants::epsilon_0 * kT_long * TMath::Qe() / (electronDensity * TMath::Qe() * TMath::Qe()));
	}

	double PlasmaFrequency(double electronDensity)
	{
		return sqrt(electronDensity * TMath::Qe() * TMath::Qe() / (PhysicalConstants::epsilon_0 * PhysicalConstants::electronMass));
	}

	double CyclotronFrequency(double magneticField)
	{
		return TMath::Qe() * magneticField / PhysicalConstants::electronMass;
	}

	double CyclotronRadius(double magneticField, double transverseVelocity)
	{
		return transverseVelocity / CyclotronFrequency(magneticField);
	}

	namespace JSPEC
	{
		static thread_local UniformCylinder beam(1, 1);
		static thread_local ForcePark parkhomchukModel;
		static thread_local ForceNonMagNumeric3D nonMagNum3DModel(1000);
		static thread_local ForceDSM derbenovSkrinskyModel;

		double ForceZ_Parkhomchuk(const Parameter& parameter)
		{
			double longSigma = Math::LongTempToVelocitySpread(parameter.kT_long);
			double transSigma = Math::TransTempToVelocitySpread(parameter.kT_trans);

			beam.set_v_rms(transSigma, longSigma);
			parkhomchukModel.set_v_eff(parameter.effectiveVelocity);
			parkhomchukModel.set_mag_field(parameter.magneticField);
			parkhomchukModel.set_time_cooler(parameter.coolerTime);
			std::vector<double> v_tr = { parameter.relativeVelocity.Perp() };
			std::vector<double> v_l = { parameter.relativeVelocity.z() };
			std::vector<double> n_e = { parameter.electronDensity };
			std::vector<double> f_tr;
			std::vector<double> f_l;
			parkhomchukModel.friction_force(parameter.ionCharge, 1, v_tr, v_l, n_e, beam, f_tr, f_l);

			return f_l.at(0) / TMath::Qe();
		}

		double ForceZ_NonMagNumeric3D(const Parameter& parameter)
		{
			double longSigma = Math::LongTempToVelocitySpread(parameter.kT_long);
			double transSigma = Math::TransTempToVelocitySpread(parameter.kT_trans);

			beam.set_v_rms(transSigma, longSigma);
			//nonMagNum3DModel.set_mag_field(parameter.magneticField);
			nonMagNum3DModel.set_time_cooler(parameter.coolerTime);
			nonMagNum3DModel.set_esprel(parameter.relTolerance);
			nonMagNum3DModel.set_gsl(true);
			std::vector<double> v_tr = { parameter.relativeVelocity.Perp() };
			std::vector<double> v_l = { parameter.relativeVelocity.z() };
			std::vector<double> n_e = { parameter.electronDensity };
			std::vector<double> f_tr;
			std::vector<double> f_l;
			nonMagNum3DModel.friction_force(parameter.ionCharge, 1, v_tr, v_l, n_e, beam, f_tr, f_l);
			
			return f_l.at(0) / TMath::Qe();
		}

		double ForceZ_DerbenovSkrinsky(const Parameter& parameter)
		{
			double longSigma = Math::LongTempToVelocitySpread(parameter.kT_long);
			double transSigma = Math::TransTempToVelocitySpread(parameter.kT_trans);

			beam.set_v_rms(transSigma, longSigma);
			derbenovSkrinskyModel.set_mag_field(parameter.magneticField);
			derbenovSkrinskyModel.set_time_cooler(parameter.coolerTime);
			derbenovSkrinskyModel.set_smooth_factor(parameter.smoothingFactor);
			derbenovSkrinskyModel.set_mag_only(!parameter.adiabatic && !parameter.fast);
			derbenovSkrinskyModel.set_steps(500);
			derbenovSkrinskyModel.set_grid(70, 70, 30);
			std::vector<double> v_tr = { parameter.relativeVelocity.Perp() };
			std::vector<double> v_l = { parameter.relativeVelocity.z() };
			std::vector<double> n_e = { parameter.electronDensity };
			std::vector<double> f_tr;
			std::vector<double> f_l;
			derbenovSkrinskyModel.friction_force(parameter.ionCharge, 1, v_tr, v_l, n_e, beam, f_tr, f_l);

			return f_l.at(0) / TMath::Qe();
		}
	}

	namespace Betacool
	{
		static thread_local xFrParam betacoolParams;
		static thread_local xForce betacoolForce;

		double ForceZ_Parkhomchuk(const Parameter& parameter)
		{
			betacoolParams.mfield = parameter.magneticField;
			betacoolParams.Z = parameter.ionCharge;
			betacoolParams.n_e = parameter.electronDensity;
			betacoolParams.tau = parameter.coolerTime;
			betacoolParams.V_long_e = Math::LongTempToVelocitySpread(parameter.kT_long);
			betacoolParams.V_tr_e = Math::TransTempToVelocitySpread(parameter.kT_trans);
			betacoolParams.V_eff_e = parameter.effectiveVelocity;
			betacoolForce.v[0] = parameter.relativeVelocity.x();
			betacoolForce.v[1] = parameter.relativeVelocity.y();
			betacoolForce.v[2] = parameter.relativeVelocity.z();
			betacoolForce.Vtr = parameter.relativeVelocity.Perp();
			
			betacoolForce.Parhom(betacoolParams);

			return betacoolForce.f[2].v;
		}

		double ForceZ_NonMag(const Parameter& parameter)
		{
			betacoolParams.Z = parameter.ionCharge;
			betacoolParams.n_e = parameter.electronDensity;
			betacoolParams.tau = parameter.coolerTime;
			betacoolParams.V_long_e = Math::LongTempToVelocitySpread(parameter.kT_long);
			betacoolParams.V_tr_e = Math::TransTempToVelocitySpread(parameter.kT_trans);
			betacoolForce.v[0] = parameter.relativeVelocity.x();
			betacoolForce.v[1] = parameter.relativeVelocity.y();
			betacoolForce.v[2] = parameter.relativeVelocity.z();
			betacoolForce.Vtr = parameter.relativeVelocity.Perp();
			// steps for manual integration
			betacoolForce.dt = 40;
			betacoolForce.dl = 40;
			betacoolForce.nfi = 20;

			betacoolForce.NonMag(betacoolParams);

			return betacoolForce.f[2].v;
		}

		double ForceZ_NonMagNumeric3D(const Parameter& parameter)
		{
			betacoolParams.Z = parameter.ionCharge;
			betacoolParams.n_e = parameter.electronDensity;
			betacoolParams.tau = parameter.coolerTime;
			betacoolParams.V_long_e = Math::LongTempToVelocitySpread(parameter.kT_long);
			betacoolParams.V_tr_e = Math::TransTempToVelocitySpread(parameter.kT_trans);
			betacoolParams.V_tr_x = betacoolParams.V_tr_e / sqrt(2);
			betacoolParams.V_tr_y = betacoolParams.V_tr_e / sqrt(2);
			betacoolForce.v[0] = parameter.relativeVelocity.x();
			betacoolForce.v[1] = parameter.relativeVelocity.y();
			betacoolForce.v[2] = parameter.relativeVelocity.z();
			betacoolForce.Vtr = parameter.relativeVelocity.Perp();
			betacoolForce.D3dl = 30;
			betacoolForce.D3dx = 30;
			betacoolForce.D3dy = 30;

			betacoolForce.D4(betacoolParams);

			return betacoolForce.f[2].v;
		}

		double ForceZ_DerbenovSkrinsky(const Parameter& parameter)
		{
			betacoolParams.mfield = parameter.magneticField;
			betacoolParams.Z = parameter.ionCharge;
			betacoolParams.n_e = parameter.electronDensity;
			betacoolParams.Smoos = parameter.smoothingFactor;
			betacoolParams.tau = parameter.coolerTime;
			betacoolParams.V_long_e = Math::LongTempToVelocitySpread(parameter.kT_long);
			betacoolParams.V_tr_e = Math::TransTempToVelocitySpread(parameter.kT_trans);
			betacoolForce.v[0] = parameter.relativeVelocity.x();
			betacoolForce.v[1] = parameter.relativeVelocity.y();
			betacoolForce.v[2] = parameter.relativeVelocity.z();
			betacoolForce.Vtr = parameter.relativeVelocity.Perp();
			betacoolForce.Magnetized = parameter.magnetic;
			betacoolForce.Fast = parameter.fast;
			betacoolForce.Adiabatic = parameter.adiabatic;
			// steps for manual integration
			betacoolForce.dt = 70;
			betacoolForce.dl = 70;
			betacoolForce.nfi = 30;
			betacoolForce.N_M = 70;

			//betacoolForce.Pestrikov = 1;
			//betacoolForce.Constant = 1;
			//betacoolForce.nfiP = 70;

			betacoolForce.DerSkr(betacoolParams);

			return betacoolForce.f[2].v;
		}

		double ForceZ_Toepffler(const Parameter& parameter)
		{
			betacoolParams.mfield = parameter.magneticField;
			betacoolParams.Z = parameter.ionCharge;
			betacoolParams.n_e = parameter.electronDensity;
			betacoolParams.tau = parameter.coolerTime;
			betacoolParams.V_long_e = Math::LongTempToVelocitySpread(parameter.kT_long);
			betacoolParams.V_tr_e = Math::TransTempToVelocitySpread(parameter.kT_trans);
			betacoolForce.v[0] = parameter.relativeVelocity.x();
			betacoolForce.v[1] = parameter.relativeVelocity.y();
			betacoolForce.v[2] = parameter.relativeVelocity.z();
			betacoolForce.Vtr = parameter.relativeVelocity.Perp();
			betacoolForce.N_M = 30;
			// for toepper
			betacoolForce.TFast = 1;
			betacoolForce.Tight = 1;
			betacoolForce.Stretched = 1;
			betacoolForce.Tdl = 30;
			betacoolForce.Tdt = 30;
			betacoolForce.Tnfi = 30;

			betacoolForce.Toepffer(betacoolParams);

			return betacoolForce.f[2].v;
		}
	}
}
}



