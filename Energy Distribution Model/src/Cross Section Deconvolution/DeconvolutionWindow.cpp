#include "pch.h"

#include "DeconvolutionWindow.h"
#include "BoltzmannDistribution.h"
#include "EnergyDistributionSet.h"
#include "EnergyDistributionWindow.h"

#include "FileUtils.h"
#include "Eigen/SVD"
#include <Constants.h>

namespace DeconvolutionWindow
{
	// main data
	static std::vector<CrossSection> crossSectionList;
	static std::vector<RateCoefficient> rateCoefficientList;
	static std::vector<PlasmaRateCoefficient> plasmaRateCoefficientList;
	static int currentCrossSectionIndex;
	static int currentRateCoefficientIndex;
	static int currentPlasmaRateCoefficientIndex;

	static CrossSectionBinningSettings binSettings;
	static FittingOptions fitSettings;

	// plot parameters
	static bool logX = true;
	static bool logY = true;
	static bool showMarkers = false;
	static bool showSubfunctions = false;
	 
	static bool showBoltzmannConvolutionWindow = false;
	static bool showBinningSettingsWindow = false;
	static bool showFitSettingsWindow = false;
	static bool showPlasmaRateWindow = false;
	 
	static char CSnameInput[64] = "cross section name";
	static char RCnameInput[64] = "rate coefficient name";

	// scale for 1/E cs
	static int scale = 1;

	static const ImVec4 inputTextColor = ImVec4(0.6f, 0.2f, 0.1f, 1.0f);

	void Init()
	{
		BoltzmannDistribution::Update(nullptr);

		// temp: load default rate coefficient
		RateCoefficient rc;
		rc.Load(FileUtils::GetMeasuredRateCoefficientFolder() / "amb-Ed_IeXXX_t0.5-14.0_pN0gt0_c00 - reversed.dat");

		AddRateCoefficientToList(rc);
	}

	void AddRateCoefficientToList(RateCoefficient& rc)
	{
		rateCoefficientList.emplace_back(std::move(rc));
		if(rateCoefficientList.size() == 1)
			currentRateCoefficientIndex = 0;
	}

	void RemoveRateCoefficient(int index)
	{
		rateCoefficientList.erase(rateCoefficientList.begin() + index);
		currentRateCoefficientIndex = std::min(currentRateCoefficientIndex, (int)rateCoefficientList.size() - 1);
	}

	void AddCrossSectionToList(CrossSection& cs)
	{
		crossSectionList.emplace_back(std::move(cs));
		currentCrossSectionIndex = crossSectionList.size() - 1;
	}

	void RemoveCrossSection(int index)
	{
		crossSectionList.erase(crossSectionList.begin() + index);
		currentCrossSectionIndex = std::min(currentCrossSectionIndex, (int)crossSectionList.size() - 1);
	}

	void AddPlasmaRateToList(PlasmaRateCoefficient& prc)
	{
		plasmaRateCoefficientList.emplace_back(std::move(prc));
		currentPlasmaRateCoefficientIndex = plasmaRateCoefficientList.size() - 1;
	}

	void RemovePlasmaRate(int index)
	{
		plasmaRateCoefficientList.erase(plasmaRateCoefficientList.begin() + index);
		currentPlasmaRateCoefficientIndex = std::min(currentPlasmaRateCoefficientIndex, (int)plasmaRateCoefficientList.size() - 1);

	}

	void ShowWindow()
	{
		if (ImGui::Begin("Deconvolution Window"))
		{
			ImGui::BeginGroup();
			EnergyDistributionWindow::ShowSetList();
			ShowRateCoefficientList();
			ShowCrossSectionList();
			ShowPlasmaRateList();
			ImGui::EndGroup();

			ImGui::SameLine();

			ImGui::BeginGroup();
			ShowSettings();
			ShowPlots();
			ImGui::EndGroup();

			CrossSection* currentCrosssection = crossSectionList.empty() ? nullptr : &crossSectionList.at(currentCrossSectionIndex);
			BoltzmannDistribution::ShowWindow(showBoltzmannConvolutionWindow, currentCrosssection);

			binSettings.ShowWindow(showBinningSettingsWindow);
			fitSettings.ShowWindow(showFitSettingsWindow);

			ShowPlasmaRateWindow();
		}
		ImGui::End();
	}

	void ShowSettings()
	{
		std::vector<EnergyDistributionSet>& setList = EnergyDistributionWindow::GetSetList();
		int currentSetIndex = EnergyDistributionWindow::GetCurrentSetIndex();

		ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX | ImGuiChildFlags_ResizeY;
		if (ImGui::BeginChild("DeconvolveSettings", ImVec2(100.0f, 0.0f), flags))
		{
			ImGui::Text("energy distribution set: "); ImGui::SameLine();
			ImGui::TextColored(inputTextColor, "%s",
				setList.empty() ? "" : setList.at(currentSetIndex).Label().c_str());

			ImGui::Text("target rate coefficient: "); ImGui::SameLine();
			ImGui::TextColored(inputTextColor, "%s",
				rateCoefficientList.empty() ? "" : rateCoefficientList.at(currentRateCoefficientIndex).GetLabel().c_str());

			ImGui::Checkbox("show binning settings", &showBinningSettingsWindow);
			ImGui::SameLine();
			ImGui::Checkbox("show Fit settings", &showFitSettingsWindow);

			ImGui::SetNextItemWidth(150.0f);
			ImGui::InputText("output CS name", CSnameInput, sizeof(CSnameInput));

			if (ImGui::Button("Deconvolve Cross Section"))
			{
				EnergyDistributionSet& currentSet = EnergyDistributionWindow::GetCurrentSet();
				RateCoefficient& currentRC = rateCoefficientList.at(currentRateCoefficientIndex);

				CrossSection cs;
				cs.SetLabel(CSnameInput);
				cs.Deconvolve(currentRC, currentSet, fitSettings, binSettings);
				cs.Save();
				AddCrossSectionToList(cs);
			}
		}
		ImGui::EndChild();
		
		ImGui::SameLine();
		if (ImGui::BeginChild("ConvolveSettings", ImVec2(100.0f, 0.0f), flags))
		{
			ImGui::Text("energy distribution set: "); ImGui::SameLine();
			ImGui::TextColored(inputTextColor, "%s",
				setList.empty() ? "" : setList.at(currentSetIndex).Label().c_str());

			ImGui::Text("cross section: "); ImGui::SameLine();
			ImGui::TextColored(inputTextColor, "%s",
				crossSectionList.empty() ? "" : crossSectionList.at(currentCrossSectionIndex).GetLabel().c_str());

			ImGui::SetNextItemWidth(150.0f);
			ImGui::InputText("output RC name", RCnameInput, sizeof(RCnameInput));

			if (ImGui::Button("Convolve Rate Coefficient"))
			{
				EnergyDistributionSet& currentSet = EnergyDistributionWindow::GetCurrentSet();
				CrossSection& currentCS = crossSectionList.at(currentCrossSectionIndex);

				RateCoefficient rc;
				rc.Convolve(currentCS, currentSet);
				rc.SetLabel(RCnameInput);
				rc.Save();
				AddRateCoefficientToList(rc);
			}
		}
		ImGui::EndChild();
	}

	void ShowPlots()
	{
		if (ImPlot::BeginPlot("rate coefficient"))
		{
			ImPlot::SetupAxis(ImAxis_X1, "detuning energy [eV]");
			ImPlot::SetupAxis(ImAxis_Y1, "rate coefficient");
			if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
			if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
			ImPlot::SetupLegend(ImPlotLocation_NorthEast);

			int i = 0;
			for (const RateCoefficient& rc : rateCoefficientList)
			{
				ImGui::PushID(i++);
				if (showSubfunctions)
				{
					rc.PlotSubfunctions();
				}
				rc.Plot(showMarkers);
				ImGui::PopID();
			}

			ImPlot::EndPlot();
		}

		ImGui::Checkbox("log X", &logX);
		ImGui::SameLine();
		ImGui::Checkbox("log Y", &logY);
		ImGui::SameLine();
		ImGui::Checkbox("show markers", &showMarkers);
		ImGui::SameLine();
		ImGui::Checkbox("show f_pl", &showBoltzmannConvolutionWindow);
		ImGui::SameLine();
		ImGui::Checkbox("show subfunctions", &showSubfunctions);
		ImGui::SameLine();
		ImGui::Checkbox("show plasma rates", &showPlasmaRateWindow);

		if (ImPlot::BeginPlot("cross section"))
		{
			ImPlot::SetupAxis(ImAxis_X1, "collision energy [eV]");
			ImPlot::SetupAxis(ImAxis_Y1, "sigma(E)");
			if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
			if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
			ImPlot::SetupLegend(ImPlotLocation_NorthEast);
			int i = 0;
			for (const CrossSection& cs : crossSectionList)
			{
				ImGui::PushID(i++);
				cs.Plot(showMarkers);
				ImGui::PopID();
			}

			if (showBoltzmannConvolutionWindow)
			{
				BoltzmannDistribution::Plot(showMarkers);
			}

			ImPlot::EndPlot();
		}
	}

	void ShowPlasmaRateWindow()
	{
		if (!showPlasmaRateWindow)
		{
			return;
		}
		if (ImGui::Begin("plasma rate window", &showPlasmaRateWindow, ImGuiWindowFlags_NoDocking))
		{
			if (ImPlot::BeginPlot("plasma rates"))
			{
				ImPlot::SetupAxis(ImAxis_X1, "Temperature [K]");
				ImPlot::SetupAxis(ImAxis_Y1, "plasma rate");
				if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
				if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
				ImPlot::SetupLegend(ImPlotLocation_NorthEast);
				int i = 0;
				for (const PlasmaRateCoefficient& prc : plasmaRateCoefficientList)
				{
					ImGui::PushID(i++);
					prc.Plot(showMarkers);
					ImGui::PopID();
				}

				ImPlot::EndPlot();
			}

			PlasmaRateCoefficient::ShowConvolutionParamterInputs();
		}
		ImGui::End();
	}

	void ShowRateCoefficientList()
	{
		ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX | ImGuiChildFlags_ResizeY;
		if (ImGui::BeginChild("merged beam rate coefficients", ImVec2(100, 100), flags))
		{
			if (ImGui::BeginListBox("##mbrclist", ImVec2(-1, 150)))
			{
				for (int i = 0; i < rateCoefficientList.size(); i++)
				{
					RateCoefficient& rc = rateCoefficientList.at(i);

					ImGui::PushID(i);
					bool selected = i == currentRateCoefficientIndex;
					if (ImGui::Selectable(rc.GetLabel().c_str(), selected, ImGuiSelectableFlags_AllowItemOverlap))
					{
						currentRateCoefficientIndex = i;
					}
					ImGui::SameLine();
					if (ImGui::SmallButton("x"))
					{
						RemoveRateCoefficient(i);
					}
					ImGui::PopID();

				}
				ImGui::EndListBox();
			}
			if (ImGui::Button("load measured rc"))
			{
				std::filesystem::path file = FileUtils::SelectFile(FileUtils::GetMeasuredRateCoefficientFolder(), {"*.dat"});
				if (!file.empty())
				{
					RateCoefficient rc;
					rc.Load(file);
					
					AddRateCoefficientToList(rc);
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("load fitted rc"))
			{
				std::vector<std::filesystem::path> files = FileUtils::SelectFiles(FileUtils::GetRateCoefficientFolder(), {"*.dat"});
				for (std::filesystem::path& file : files)
				{
					RateCoefficient rc;
					rc.Load(file);
					AddRateCoefficientToList(rc);
				}
			}
		}
		ImGui::EndChild();
	}

	void ShowCrossSectionList()
	{
		ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX | ImGuiChildFlags_ResizeY;
		if (ImGui::BeginChild("cross sections", ImVec2(100, 100), flags))
		{
			if (ImGui::BeginListBox("##cslist", ImVec2(-1, 150)))
			{
				for (int i = 0; i < crossSectionList.size(); i++)
				{
					CrossSection& cs = crossSectionList.at(i);

					ImGui::PushID(i);
					bool selected = i == currentCrossSectionIndex;
					if (ImGui::Selectable(cs.GetLabel().c_str(), selected, ImGuiSelectableFlags_AllowItemOverlap))
					{
						currentCrossSectionIndex = i;
					}
					ImGui::SameLine();
					if (ImGui::SmallButton("-> plasma rate"))
					{
						PlasmaRateCoefficient prc;
						prc.ConvolveFromErrorIterationArray(cs);
						prc.Save();
						AddPlasmaRateToList(prc);
					}
					ImGui::SameLine();
					if (ImGui::SmallButton("x"))
					{
						RemoveCrossSection(i);
					}
					ImGui::PopID();

				}
				ImGui::EndListBox();
			}
			if (ImGui::Button("load fitted cs"))
			{
				std::vector<std::filesystem::path> files = FileUtils::SelectFiles(FileUtils::GetCrossSectionFolder(), { "*.dat" });
				for (std::filesystem::path& file : files)
				{
					CrossSection cs;
					cs.Load(file);
					cs.SetLabel(file.filename().string());
					AddCrossSectionToList(cs);
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("show initial guess"))
			{
				CrossSection cs;
				cs.SetupBinning(binSettings, rateCoefficientList.at(currentRateCoefficientIndex));
				cs.SetInitialGuessValues(rateCoefficientList.at(currentRateCoefficientIndex));
				cs.SetLabel("initial guess");
				AddCrossSectionToList(cs);
			}
			if (ImGui::Button("create 1/E cs"))
			{
				CrossSection cs;
				cs.FillWithOneOverE(scale);
				AddCrossSectionToList(cs);
			}
			ImGui::SameLine();
			ImGui::SetNextItemWidth(50.0f);
			ImGui::InputInt("scale", &scale, 0, 0);

		}
		ImGui::EndChild();
	}

	void ShowPlasmaRateList()
	{
		ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX | ImGuiChildFlags_ResizeY;
		if (ImGui::BeginChild("plasma rate coefficients", ImVec2(100, 100), flags))
		{
			if (ImGui::BeginListBox("##plasmaRatelist", ImVec2(-1, 150)))
			{
				for (int i = 0; i < plasmaRateCoefficientList.size(); i++)
				{
					PlasmaRateCoefficient& prc = plasmaRateCoefficientList.at(i);

					ImGui::PushID(i);
					bool selected = i == currentPlasmaRateCoefficientIndex;
					if (ImGui::Selectable(prc.GetLabel().c_str(), selected, ImGuiSelectableFlags_AllowItemOverlap))
					{
						currentPlasmaRateCoefficientIndex = i;
					}

					ImGui::SameLine();
					if (ImGui::SmallButton("x"))
					{
						RemovePlasmaRate(i);
					}
					ImGui::PopID();

				}
				ImGui::EndListBox();
			}
			if (ImGui::Button("load fitted plasma rate"))
			{
				std::vector<std::filesystem::path> files = FileUtils::SelectFiles(FileUtils::GetPlasmaRateFolder(), { "*.dat" });
				for (std::filesystem::path& file : files)
				{
					PlasmaRateCoefficient prc;
					prc.Load(file);
					prc.SetLabel(file.filename().string());
					AddPlasmaRateToList(prc);
				}
			}
		}
		ImGui::EndChild();
	}
}


//
//double* DeconvolutionManager::DeconvolveWithGradientDescent(const RateCoefficient& rc, const EnergyDistributionSet& set, const CrossSection& cs)
//{
//	int n = set.distributions.size();
//	int p = cs.hist->GetNbinsX();
//
//	Eigen::MatrixXd PsiMatrix(n, p);
//	Eigen::VectorXd alphaVector(n);
//
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < p; j++)
//		{
//			PsiMatrix(i, j) = set.distributions[i].psi[j];
//		}
//		alphaVector[i] = rc.value[i];
//	}
//
//	Eigen::VectorXd x = Eigen::VectorXd::Map(cs.values.data(), cs.values.size());
//
//	for (int iter = 0; iter < iterations; iter++)
//	{
//		// Compute the residual: r = A*x - b
//		Eigen::VectorXd r = PsiMatrix * x - alphaVector;
//		
//		// Compute the gradient: grad = A^T * r
//		Eigen::VectorXd grad = PsiMatrix.transpose() * r;
//		//std::cout << "gradient: " << grad << std::endl;
//		//std::cout << "gradient modified: " << grad.cwiseMin(0.001 / learningRate).cwiseMax(-0.001 / learningRate) << std::endl;
//
//		// Update x using gradient descent step
//		x = x - learningRate * grad.cwiseMin(0.0001 / learningRate).cwiseMax(-0.0001 / learningRate);
//		//std::cout << "new x: " << x << std::endl;
//		x = x.cwiseAbs();
//	}
//
//	return x.data();
//}
