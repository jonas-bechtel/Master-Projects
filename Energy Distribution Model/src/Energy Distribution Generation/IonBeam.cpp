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
	static TH3D* beam;
	static TH3D* beamSmall;

	// optional parameters
	static bool doubleGaussian;
	static double amplitude2;
	static float sigma2[2];
	static bool limitZRange = false;
	static float limitedZRange[2] = { -0.4f, 0.4f };

	// plotting data
	static std::vector<double> xAxis;
	static std::vector<double> yAxis;
	static std::vector<double> zAxis;

	static std::vector<double> projectionValuesX;
	static std::vector<double> projectionValuesY;
	static std::vector<double> projectionValuesZ;

	static HeatMapData slice;
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
		beam = new TH3D("ion beam", "ion beam", 100, -0.04, 0.04, 100, -0.04, 0.04, 200, -0.7, 0.7);
		UpdatePlotData();
	}

	void CreateFromReference(TH3D* reference)
	{
		if (!reference)
		{
			std::cout << "ion beam reference is null" << std::endl;
		}

		delete beam;
		beam = (TH3D*)reference->Clone("ion density");
		beam->SetTitle("ion density");
		
		beam->Reset();
		UpdateMainData();
	}

	void UpdateMainData()
	{
		int nXBins = beam->GetXaxis()->GetNbins();
		int nYBins = beam->GetYaxis()->GetNbins();
		int nZBins = beam->GetZaxis()->GetNbins();

		for (int i = 1; i <= nXBins; i++) 
		{
			for (int j = 1; j <= nYBins; j++)
			{
				for (int k = 1; k <= nZBins; k++) 
				{
					if(parameter.sigma.get().x == 0 || parameter.sigma.get().y == 0)
					{
						beam->SetBinContent(i, j, k, 0);
						continue;
					}
					
					// Calculate the coordinates for this bin
					double x = beam->GetXaxis()->GetBinCenter(i);
					double y = beam->GetYaxis()->GetBinCenter(j);
					double z = beam->GetZaxis()->GetBinCenter(k);

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
							beam->SetBinContent(i, j, k, 0);
							continue;
						}

						double normalisation2 = 1.0 / (2 * TMath::Pi() * sigma2[0] * sigma2[1]);
						double gauss2 = normalisation2 * exp(-0.5 * ((x * x) / pow(sigma2[0], 2) + (y * y) / pow(sigma2[1], 2)));

						value = (parameter.amplitude * gauss1 + amplitude2 * gauss2) / amplitudeSum;
					}
					else
					{
						value = gauss1;
					}

					beam->SetBinContent(i, j, k, value);
				}
			}
		}
		
		//TH3D* beamSmalltest = (TH3D*)beam->Rebin3D(2, 2, 4, "ion beam small");
		////remove all bin values that are below 100
		//for (int i = 0; i <= beamSmalltest->GetNbinsX() + 1; i++)
		//{
		//	for (int j = 0; j <= beamSmalltest->GetNbinsY() + 1; j++)
		//	{
		//		for (int k = 0; k <= beamSmalltest->GetNbinsZ() + 1; k++)
		//		{
		//			if (beamSmalltest->GetBinContent(i, j, k) <= 100.0)
		//			{
		//				//std::cout << "Bin (" << i << ", " << j << ", " << k << "): " << beamSmalltest->GetBinContent(i, j, k) << "\n";
		//				beamSmalltest->SetBinContent(i, j, k, 0.0);
		//			}
  //                  //std::cout << "Bin (" << i << ", " << j << ", " << k << "): " << beamSmalltest->GetBinContent(i, j, k) << "\n";
		//		}
		//	}
		//}

		//beamSmalltest->SaveAs("ion_beam_small.C");
		//delete beamSmalltest;
	}

	void UpdatePlotData()
	{
		xAxis.clear();
		yAxis.clear();
		zAxis.clear();

		projectionValuesX.clear();
		projectionValuesY.clear();
		projectionValuesZ.clear();

		beam->Reset();
		UpdateMainData();

		xAxis.reserve(beam->GetNbinsX());
		yAxis.reserve(beam->GetNbinsY());
		zAxis.reserve(beam->GetNbinsZ());

		projectionValuesX.reserve(beam->GetNbinsX());
		projectionValuesY.reserve(beam->GetNbinsY());
		projectionValuesZ.reserve(beam->GetNbinsZ());

		for (int i = 1; i <= beam->GetNbinsX(); i++)
		{
			xAxis.push_back(beam->GetXaxis()->GetBinCenter(i));
		}
		for (int i = 1; i <= beam->GetNbinsY(); i++)
		{
			yAxis.push_back(beam->GetYaxis()->GetBinCenter(i));
		}
		for (int i = 1; i <= beam->GetNbinsZ(); i++)
		{
			zAxis.push_back(beam->GetZaxis()->GetBinCenter(i));
		}

		TH1D* projectionX = beam->ProjectionX();
		TH1D* projectionY = beam->ProjectionY();
		TH1D* projectionZ = beam->ProjectionZ();

		for (int i = 1; i <= projectionX->GetNbinsX(); i++)
		{
			projectionValuesX.push_back(projectionX->GetBinContent(i));
		}
		for (int i = 1; i <= projectionY->GetNbinsX(); i++)
		{
			projectionValuesY.push_back(projectionY->GetBinContent(i));
		}
		for (int i = 1; i <= projectionZ->GetNbinsX(); i++)
		{
			projectionValuesZ.push_back(projectionZ->GetBinContent(i));
		}

		delete projectionX;
		delete projectionY;
		delete projectionZ;

		slice.FromTH3D(beam, SliceZ);
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
					slice.FromTH3D(beam, SliceZ);
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
				ImPlot::PlotLine("", xAxis.data(), projectionValuesX.data(), xAxis.size());
				ImPlot::EndPlot();
			}

			if (ImPlot::BeginPlot("Projection Y"))
			{
				ImPlot::PlotLine("", yAxis.data(), projectionValuesY.data(), yAxis.size());

				ImPlot::EndPlot();
			}

			if (ImPlot::BeginPlot("Projection Z"))
			{
				ImPlot::PlotLine("", zAxis.data(), projectionValuesZ.data(), zAxis.size());

				ImPlot::EndPlot();
			}

			slice.Plot("");

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
			UpdatePlotData();
		}

		ImGui::PopItemWidth();
		ImGui::EndGroup();
	}

	void ShowCoolingForceParameterControls()
	{
		ImGui::BeginGroup();
		ImGui::PushItemWidth(170.0f);
		ImGui::InputInt("number samples", parameter.numberSamples);
		ImGui::InputInt("charge", parameter.charge);
		ImGui::BeginDisabled(!limitZRange);
		ImGui::InputFloat2("##limited range", limitedZRange);
		ImGui::EndDisabled();
		ImGui::SameLine();
		ImGui::Checkbox("limit z range", &limitZRange);
		ImGui::PopItemWidth();
		ImGui::EndGroup();
	}

	std::vector<Point3D> GeneratePositions()
	{
		std::vector<Point3D> positions;
		positions.reserve(parameter.numberSamples);

		auto now = std::chrono::system_clock::now();       
		auto seconds = now.time_since_epoch().count();

		TRandom3 generator(seconds);

		float range[2] = { beam->GetZaxis()->GetXmin(), beam->GetZaxis()->GetXmax() };
		if (limitZRange)
		{
			range[0] = limitedZRange[0];
			range[1] = limitedZRange[1];
		}

		if (!doubleGaussian)
		{
			for (int i = 0; i < parameter.numberSamples; i++)
			{
				double z = generator.Uniform(range[0], range[1]);
				double x = generator.Gaus(parameter.shift.get().x + parameter.angles.get().x * z, parameter.sigma.get().x);
				double y = generator.Gaus(parameter.shift.get().y + parameter.angles.get().y * z, parameter.sigma.get().y);
				positions.emplace_back(x, y, z);
			}
		}

		return positions;
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

	double GetValue(const Point3D& point)
	{
		double x = point.x;
		double y = point.y;
		double z = point.z;

		// apply shift of ion beam
		x -= parameter.shift.get().x;
		y -= parameter.shift.get().y;

		// apply the angles with small angle approximation
		x -= parameter.angles.get().x * z;
		y -= parameter.angles.get().y * z;

		double value = 0;
		double normalisation = 1 / (parameter.sigma.get().x * parameter.sigma.get().y * 2 * TMath::Pi());
		value = normalisation * exp(-0.5 * ((x * x) / pow(parameter.sigma.get().x, 2) + (y * y) / pow(parameter.sigma.get().y, 2)));

		if (doubleGaussian)
		{
			value += amplitude2 * exp(-0.5 * ((x * x) / pow(sigma2[0], 2) + (y * y) / pow(sigma2[1], 2)));
		}

		return value;
	}

	double GetVelocityMagnitude()
	{
		return TMath::Sqrt(2 * ElectronBeam::GetCoolingEnergy() * TMath::Qe() / PhysicalConstants::electronMass);
	}

	int GetCharge()
	{
		return parameter.charge;
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
		return beam;
	}

	std::string GetTags()
	{
		std::string tags = "";
		if (doubleGaussian) tags += "ion-2gaus, ";
		if (limitZRange) tags += Form("z samples %.3f - %.3f, ", limitedZRange[0], limitedZRange[1]);
		
		return tags;
	}
}


