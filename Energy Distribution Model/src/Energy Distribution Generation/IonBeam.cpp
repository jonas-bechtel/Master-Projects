#include "pch.h"

#include "IonBeam.h"
#include "ElectronBeam.h"
#include "MCMC.h"
#include "EnergyDistribution.h"
#include "Constants.h"

namespace IonBeam
{
	static IonBeamParameters parameter;

	// 3D Hist with main data
	static HistData3D histData;

	// optional parameters
	static bool doubleGaussian;
	static double amplitude2;
	static float sigma2[2];
	static bool limitZRange = false;
	static float limitedZRange[2] = { -0.4f, 0.4f };

	static float SliceZ = 0.0f;

	IonBeamParameters GetParameters()
	{
		return parameter;
	}

	void SetParameters(const IonBeamParameters& params)
	{
		parameter = params;
	}

	void Init()
	{
		histData = HistData3D(new TH3D("ion beam", "ion beam", 100, -0.04, 0.04, 100, -0.04, 0.04, 200, -0.7, 0.7));
		UpdateHistData();
	}

	void CreateFromReference(TH3D* reference)
	{
		if (!reference)
		{
			std::cout << "ion beam reference is null" << std::endl;
			return;
		}
		histData = HistData3D((TH3D*)reference->Clone("ion density"));
		histData.GetHist()->SetTitle("ion density");
		
		UpdateHistData();
	}

	void UpdateHistData()
	{
		TH3D* beam = histData.GetHist();

		if (!beam) 
			return;

		beam->Reset();
		int nXBins = beam->GetXaxis()->GetNbins();
		int nYBins = beam->GetYaxis()->GetNbins();
		int nZBins = beam->GetZaxis()->GetNbins();

		for (int i = 1; i <= nXBins; i++) 
		{
			for (int j = 1; j <= nYBins; j++)
			{
				for (int k = 1; k <= nZBins; k++) 
				{
					// Calculate the coordinates for this bin
					double x = beam->GetXaxis()->GetBinCenter(i);
					double y = beam->GetYaxis()->GetBinCenter(j);
					double z = beam->GetZaxis()->GetBinCenter(k);

					double value = GetValue(x, y, z);

					beam->SetBinContent(i, j, k, value);
				}
			}
		}
		histData.UpdateSlice(SliceZ);
	}

	void ShowWindow()
	{
		if (ImGui::Begin("Ion Beam"))
		{
			if (ImGui::BeginChild("settings", ImVec2(300, -1), ImGuiChildFlags_ResizeX))
			{
				ShowParameterControls();
				ImGui::SeparatorText("cooling force specific options");
				ShowCoolingForceParameterControls();
				ImGui::Separator();

				if (ImGui::SliderFloat("slice z", &SliceZ, -0.7f, 0.7f))
				{
					histData.UpdateSlice(SliceZ);
				}

			}
			ImGui::EndChild();
			ImGui::SameLine();
			ShowPlots();

		}
		ImGui::End();
	}

	void ShowPlots()
	{
		if (ImPlot::BeginSubplots("##ion beam subplots", 2, 3, ImVec2(-1, -1), ImPlotSubplotFlags_ShareItems))
		{
			if (ImPlot::BeginPlot("Projection X"))
			{
				histData.PlotProjectionX(0);
				ImPlot::EndPlot();
			}

			if (ImPlot::BeginPlot("Projection Y"))
			{
				histData.PlotProjectionY(0);
				ImPlot::EndPlot();
			}

			if (ImPlot::BeginPlot("Projection Z"))
			{
				histData.PlotProjectionZ(0);
				ImPlot::EndPlot();
			}
			histData.PlotSlice(0);

			ImPlot::EndSubplots();
		}
	}

	void ShowParameterControls()
	{
		bool somethingChanged = false;

		ImGui::BeginGroup();

		ImGui::PushItemWidth(170.0f);

		somethingChanged |= ImGui::InputFloat2("shift in x and y [m]", parameter.shift, "%.4f");
		somethingChanged |= ImGui::InputFloat2("horizontal, vertical angles [rad]", parameter.angles, "%.4f");

		ImGui::Separator();
		ImGui::BeginDisabled(!doubleGaussian);
		somethingChanged |= ImGui::InputDouble("amplitude", parameter.amplitude, 0.0f, 0.0f, "%.4f");
		ImGui::EndDisabled();
		somethingChanged |= ImGui::InputFloat2("sigmas x and y [m]", parameter.sigma, "%.5f");
		
		if (ImGui::Button("Set Emittance Values"))
		{
			somethingChanged = true;
			parameter.sigma.set({ 0.01044, 0.00455 });
		}

		ImGui::Separator();
		somethingChanged |= ImGui::Checkbox("use second gaussian", &doubleGaussian);
		ImGui::BeginDisabled(!doubleGaussian);
		somethingChanged |= ImGui::InputDouble("amplitude 2", &amplitude2, 0.0f, 0.0f, "%.4f");
		somethingChanged |= ImGui::InputFloat2("sigmas 2 x and y [m]", sigma2, "%.5f");
		
		if (ImGui::Button("Set Lucias Values"))
		{
			somethingChanged = true;
			parameter.amplitude.set(10.1);
			amplitude2 = 8.1;

			parameter.sigma.set({ 9.5e-3f, 5.71e-3f });
			sigma2[0] = 1.39e-3f;
			sigma2[1] = 2.15e-3f;

		}

		ImGui::EndDisabled();
		
		if (somethingChanged)
		{
			UpdateHistData();
		}

		ImGui::PopItemWidth();
		ImGui::EndGroup();
	}

	void ShowCoolingForceParameterControls()
	{
		ImGui::BeginGroup();
		ImGui::PushItemWidth(170.0f);
		ImGui::BeginDisabled(!limitZRange);
		ImGui::InputFloat2("##limited range", limitedZRange);
		ImGui::EndDisabled();
		ImGui::SameLine();
		ImGui::Checkbox("limit z range", &limitZRange);
		ImGui::PopItemWidth();
		ImGui::EndGroup();
	}

	TVector3 GetDirection()
	{
		float angleX = parameter.angles.get().x;
		float angleY = parameter.angles.get().y;

		return TVector3(angleX, angleY, 1).Unit();
	}

	TVector3 GetVelocity()
	{
		return GetDirection() * GetVelocityMagnitude();
	}

	double GetValue(double x, double y, double z)
	{
		if (parameter.sigma.get().x == 0 || parameter.sigma.get().y == 0)
		{
			return 0;
		}

		// apply shift of ion beam
		x -= parameter.shift.get().x;
		y -= parameter.shift.get().y;

		// apply the angles with small angle approximation
		x -= parameter.angles.get().x * z;
		y -= parameter.angles.get().y * z;

		double value = 0.0;

		double normalisation1 = 1 / (parameter.sigma.get().x * parameter.sigma.get().y * 2 * TMath::Pi());
		double gauss1 = normalisation1 * exp(-0.5 * ((x * x) / pow(parameter.sigma.get().x, 2) + (y * y) / pow(parameter.sigma.get().y, 2)));

		if (doubleGaussian)
		{
			double amplitudeSum = parameter.amplitude + amplitude2;

			if (sigma2[0] == 0 || sigma2[1] == 0 || amplitudeSum == 0)
			{
				return 0;
			}

			double normalisation2 = 1.0 / (2 * TMath::Pi() * sigma2[0] * sigma2[1]);
			double gauss2 = normalisation2 * exp(-0.5 * ((x * x) / pow(sigma2[0], 2) + (y * y) / pow(sigma2[1], 2)));

			value = (parameter.amplitude * gauss1 + amplitude2 * gauss2) / amplitudeSum;
		}
		else
		{
			value = gauss1;
		}

		return value;
	}

	double GetVelocityMagnitude()
	{
		return TMath::Sqrt(2 * ElectronBeam::GetCoolingEnergy() * TMath::Qe() / PhysicalConstants::electronMass);
	}

	float GetSigmaX()
	{
		return parameter.sigma.get().x;
	}

	float GetSigmaY()
	{
		return parameter.sigma.get().y;
	}

	float* GetLimitedRange()
	{
		return limitedZRange;
	}

	bool IsRangeLimited()
	{
		return limitZRange;
	}

	TH3D* Get()
	{
		return histData.GetHist();
	}

	std::string GetTags()
	{
		std::string tags = "";
		if (doubleGaussian) tags += "ion-2gaus, ";
		if (limitZRange) tags += Form("z samples %.3f - %.3f, ", limitedZRange[0], limitedZRange[1]);
		
		return tags;
	}
}


