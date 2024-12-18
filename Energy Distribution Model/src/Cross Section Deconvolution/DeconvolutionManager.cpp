#include "pch.h"
#include "DeconvolutionManager.h"
#include "RateCoefficient.h"
#include "FileHandler.h"
#include "CrossSectionManager.h"

#include "Eigen/SVD"

DeconvolutionManager::DeconvolutionManager()
	: CrossSectionDeconvolutionModule("Deconvolution")
{
	deconvolutionManager = this;
}

void DeconvolutionManager::ShowRateCoefficientListWindow()
{
	if (ImGui::Begin("merged beam rate coefficients"))
	{
		if (ImGui::BeginListBox("##mbrclist", ImVec2(-1, 200)))
		{
			for (int i = 0; i < rateCoefficientList.size(); i++)
			{
				RateCoefficient& rc = rateCoefficientList.at(i);

				ImGui::PushID(i);
				bool selected = i == currentRateCoefficientIndex;
				if (ImGui::Selectable(rc.label.c_str(), selected, ImGuiSelectableFlags_AllowItemOverlap))
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
			std::filesystem::path file = FileHandler::GetInstance().SelectFile("data\\RateCoefficients\\", { "*.dat" });
			if (!file.empty())
			{
				RateCoefficient rc = FileHandler::GetInstance().LoadRateCoefficients(file);
				rc.measured = true;
				rc.file = file.filename();
				rc.label = file.filename().string();
				AddRateCoefficientToList(rc);

				PlotRateCoefficient();
			}
		}
		ImGui::SameLine();
		if (ImGui::Button("load fitted rc"))
		{
			std::filesystem::path file = FileHandler::GetInstance().SelectFile("Output\\Rate Coefficients\\", { "*.dat" });
			if (!file.empty())
			{
				RateCoefficient rc = FileHandler::GetInstance().LoadRateCoefficients(file);
				rc.measured = false;
				rc.file = file.filename();
				rc.label = file.filename().string();
				AddRateCoefficientToList(rc);

				PlotRateCoefficient();
			}
		}

		ImGui::End();
	}
}

void DeconvolutionManager::ShowCrossSectionListWindow()
{
	if (ImGui::Begin("cross sections"))
	{
		if (ImGui::BeginListBox("##cslist", ImVec2(-1, 200)))
		{
			for (int i = 0; i < crossSectionList.size(); i++)
			{
				CrossSection& cs = crossSectionList.at(i);

				ImGui::PushID(i);
				bool selected = i == currentCrossSectionIndex;
				if (ImGui::Selectable(cs.label.c_str(), selected, ImGuiSelectableFlags_AllowItemOverlap))
				{
					currentCrossSectionIndex = i;
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
			std::filesystem::path file = FileHandler::GetInstance().SelectFile("Output\\Cross Sections\\", { "*.dat" });
			if (!file.empty())
			{
				CrossSection cs = FileHandler::GetInstance().LoadCrossSection(file);
				cs.file = file.filename();
				cs.label = file.filename().string();
				AddCrossSectionToList(cs);
			}
		}
		ImGui::End();
	}
}

void DeconvolutionManager::ShowUI()
{
	ShowSettings();
	ShowPlots();
}

void DeconvolutionManager::ShowSettings()
{
	ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX | ImGuiChildFlags_ResizeY;
	if (ImGui::BeginChild("DeconvolveSettings", ImVec2(100.0f, 0.0f), flags))
	{
		ImGui::Text("energy distribution set: "); ImGui::SameLine();
		ImGui::TextColored(inputColor, "%s",
			energyDistributionSets.empty() ? "" : energyDistributionSets.at(currentSetIndex).Label().c_str());

		ImGui::Text("target rate coefficient: "); ImGui::SameLine();
		ImGui::TextColored(inputColor, "%s",
			rateCoefficientList.empty() ? "" : rateCoefficientList.at(currentRateCoefficientIndex).label.c_str());

		ImGui::Button("binning options");
		SetupBinningOptionsPopup();
		ImGui::OpenPopupOnItemClick("binning options", ImGuiPopupFlags_MouseButtonLeft);
		ImGui::SameLine();
		ImGui::Button("fitting options");
		SetupFitOptionsPopup();
		ImGui::OpenPopupOnItemClick("fit options", ImGuiPopupFlags_MouseButtonLeft);

		ImGui::SetNextItemWidth(150.0f);
		ImGui::InputText("output CS name", CSnameInput, sizeof(CSnameInput));

		if (ImGui::Button("Deconvolve Cross Section"))
		{
			EnergyDistributionSet& currentSet = energyDistributionSets.at(currentSetIndex);
			RateCoefficient& currentRC = rateCoefficientList.at(currentRateCoefficientIndex);

			if (!crossSectionList.empty() && crossSectionList.at(currentCrossSectionIndex).label == CSnameInput)
			{
				CrossSection& currentCS = crossSectionList.at(currentCrossSectionIndex);
				DeconvolveInPlace(currentRC, currentSet, currentCS);
			}
			else
			{
				CrossSection cs = Deconvolve(currentRC, currentSet);
				AddCrossSectionToList(cs);
			}
		}
		
		ImGui::EndChild();
	}
	ImGui::SameLine();
	if (ImGui::BeginChild("ConvolveSettings", ImVec2(100.0f, 0.0f), flags))
	{
		ImGui::Text("energy distribution set: "); ImGui::SameLine();
		ImGui::TextColored(inputColor, "%s", 
			energyDistributionSets.empty() ? "" : energyDistributionSets.at(currentSetIndex).Label().c_str());

		ImGui::Text("cross section: "); ImGui::SameLine();
		ImGui::TextColored(inputColor, "%s",
			crossSectionList.empty() ? "" : crossSectionList.at(currentCrossSectionIndex).label.c_str());

		ImGui::SetNextItemWidth(150.0f);
		ImGui::InputText("output RC name", RCnameInput, sizeof(RCnameInput));

		if (ImGui::Button("Convolve Rate Coefficient"))
		{
			EnergyDistributionSet& currentSet = energyDistributionSets.at(currentSetIndex);
			CrossSection& currentCS = crossSectionList.at(currentCrossSectionIndex);
			
			if (!rateCoefficientList.empty() && rateCoefficientList.at(currentRateCoefficientIndex).label == RCnameInput)
			{
				RateCoefficient& currentRC = rateCoefficientList.at(currentRateCoefficientIndex);
				ConvolveInPlace(currentCS, currentSet, currentRC);
			}
			else
			{
				RateCoefficient rc = Convolve(currentCS, currentSet);
				AddRateCoefficientToList(rc);
			}
			
		}
		ImGui::EndChild();
	}
}

void DeconvolutionManager::ShowPlots()
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
			if (showMarkers) ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
			ImPlot::PlotLine(rc.label.c_str(), rc.detuningEnergies.data(), rc.value.data(), rc.value.size());
			ImPlot::PlotErrorBars(rc.label.c_str(), rc.detuningEnergies.data(), rc.value.data(), rc.error.data(), rc.error.size());
			ImGui::PopID();
		}

		ImPlot::EndPlot();
	}

	ImGui::Checkbox("log X", &logX);
	ImGui::SameLine();
	ImGui::Checkbox("log Y", &logY);
	ImGui::SameLine();
	ImGui::Checkbox("show markers", &showMarkers);

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
			if (showMarkers) ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
			ImPlot::PlotLine(cs.label.c_str(), cs.energies.data(), cs.values.data(), cs.values.size());
			ImPlot::PlotErrorBars(cs.label.c_str(), cs.energies.data(), cs.values.data(), cs.errors.data(), cs.errors.size());

			ImPlot::PlotLine((cs.label + " init").c_str(), cs.energies.data(), cs.initialGuess.data(), cs.initialGuess.size());
			ImGui::PopID();
		}

		ImPlot::EndPlot();
	}
}

void DeconvolutionManager::SetupFitOptionsPopup()
{
	if (ImGui::BeginPopupContextItem("fit options"))
	{
		ImGui::Checkbox("ROOT fitting", &ROOT_fit);
		ImGui::SameLine();
		ImGui::Checkbox("limit to positive parameters", &limitROOTparameterRange);
		ImGui::Checkbox("SVD fitting", &SVD_fit);
		ImGui::Checkbox("GD fitting", &GD_fit);
		ImGui::InputInt("iterations", &iterations);
		ImGui::InputDouble("learning rate", &learningRate);
		
		ImGui::EndPopup();
	}

}

void DeconvolutionManager::SetupBinningOptionsPopup()
{
	if (ImGui::BeginPopupContextItem("binning options"))
	{
		ImGui::Combo("binning options", (int*)&currentSettings.scheme, binningOptions, IM_ARRAYSIZE(binningOptions));

		if (currentSettings.scheme == FactorBinning || currentSettings.scheme == PaperFactorMix)
		{
			ImGui::SameLine();
			ImGui::SetNextItemWidth(100.0f);
			ImGui::InputInt("number bins", &currentSettings.numberBins);
		}
		if (currentSettings.scheme == FactorBinning)
		{
			//ImGui::SameLine();
			//ImGui::Checkbox("limit bin size", &limitBinSize);
			//ImGui::SameLine();
			//ImGui::BeginDisabled(!limitBinSize);
			//ImGui::SetNextItemWidth(100.0f);
			//ImGui::InputDouble("min bin size", &minBinSize, 0.0, 0.0, "%.1e");
			//ImGui::EndDisabled();
		}
		if (currentSettings.scheme == PaperBinning)
		{
			ImGui::SameLine();
			ImGui::InputInt("max ration", &currentSettings.maxRatio);
			//ImGui::InputInt("factor", &binFactor);
		}
		ImGui::EndPopup();
	}
}

CrossSection DeconvolutionManager::Deconvolve(const RateCoefficient& rc, EnergyDistributionSet& set)
{
	if (rc.value.size() != set.distributions.size())
	{
		std::cout << "sizes of rate coefficients and energy distributions dont match: " << 
			rc.value.size() << " != " << set.distributions.size() << std::endl;
	}
	CrossSection cs = CrossSection();
	cs.energyDistriubtionSetFolder = set.folder / set.subFolder;
	cs.mergedBeamRateCoefficientFile = rc.file;
	cs.label = CSnameInput;
	cs.file = cs.label + ".dat";

	cs.SetupBinning(currentSettings, rc);

	set.CalculatePsisFromBinning(cs.hist);
	cs.SetupInitialGuess(rc, ROOT_fit);
	DeconvolveInPlace(rc, set, cs);

	return cs;
}

void DeconvolutionManager::DeconvolveInPlace(const RateCoefficient& rc, const EnergyDistributionSet& set, CrossSection& cs)
{
	if (ROOT_fit)
	{
		double* parameter = DeconvolveWithROOT(rc, set, cs);
		cs.SetValues(parameter, true);
	}
	else if (SVD_fit)
	{
		double* parameter = DeconvolveWithSVD(rc, set, cs);
		cs.SetValues(parameter);
	}
	else if (GD_fit)
	{
		double* parameter = DeconvolveWithGradientDescent(rc, set, cs);
		cs.SetValues(parameter);
	}

	FileHandler::GetInstance().SaveCrossSection(cs);
}

double* DeconvolutionManager::DeconvolveWithSVD(const RateCoefficient& rc, const EnergyDistributionSet& set, const CrossSection& cs)
{
	int n = set.distributions.size();
	int p = cs.hist->GetNbinsX();

	Eigen::MatrixXd PsiMatrix(n, p);
	Eigen::VectorXd alphaVector(n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < p; j++)
		{
			// fill matrix and vector with 0 if p > n
			PsiMatrix(i, j) = set.distributions[i].psi[j];
		}
		alphaVector[i] = rc.value[i];
	}

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(PsiMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);

	// Solve using the pseudoinverse (x = A^+ * b), where A^+ is the Moore-Penrose pseudoinverse
	Eigen::VectorXd result = svd.solve(alphaVector);

	return result.data();
}

double* DeconvolutionManager::DeconvolveWithROOT(const RateCoefficient& rc, const EnergyDistributionSet& set, const CrossSection& cs)
{
	TF1* fitFunction = new TF1("fit function", this, &DeconvolutionManager::ConvolveFit, 0, 100, cs.hist->GetNbinsX());

	//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSL");

	fitFunction->SetParameters(cs.values.data());
	if (limitROOTparameterRange)
	{
		for (int i = 0; i < cs.hist->GetNbinsX(); i++)
		{
			fitFunction->SetParLimits(i, 0, 1e30);
		}
	}

	rc.graph->Fit(fitFunction, "RN");

	return fitFunction->GetParameters();
}

double* DeconvolutionManager::DeconvolveWithGradientDescent(const RateCoefficient& rc, const EnergyDistributionSet& set, const CrossSection& cs)
{
	int n = set.distributions.size();
	int p = cs.hist->GetNbinsX();

	Eigen::MatrixXd PsiMatrix(n, p);
	Eigen::VectorXd alphaVector(n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < p; j++)
		{
			PsiMatrix(i, j) = set.distributions[i].psi[j];
		}
		alphaVector[i] = rc.value[i];
	}

	Eigen::VectorXd x = Eigen::VectorXd::Map(cs.values.data(), cs.values.size());

	for (int iter = 0; iter < iterations; iter++)
	{
		// Compute the residual: r = A*x - b
		Eigen::VectorXd r = PsiMatrix * x - alphaVector;
		
		// Compute the gradient: grad = A^T * r
		Eigen::VectorXd grad = PsiMatrix.transpose() * r;
		//std::cout << "gradient: " << grad << std::endl;
		//std::cout << "gradient modified: " << grad.cwiseMin(0.001 / learningRate).cwiseMax(-0.001 / learningRate) << std::endl;

		// Update x using gradient descent step
		x = x - learningRate * grad.cwiseMin(0.0001 / learningRate).cwiseMax(-0.0001 / learningRate);
		//std::cout << "new x: " << x << std::endl;
		x = x.cwiseAbs();
	}

	return x.data();
}

RateCoefficient DeconvolutionManager::Convolve(const CrossSection& cs, EnergyDistributionSet& set)
{
	RateCoefficient rc = RateCoefficient();
	rc.energyDistriubtionSetFolder = set.folder / set.subFolder;
	rc.crossSectionFile = cs.file;
	rc.measured = false;
	rc.label = RCnameInput;
	rc.file = rc.label + ".dat";

	for (const EnergyDistribution& eDist : set.distributions)
	{
		rc.detuningEnergies.push_back(eDist.eBeamParameter.detuningEnergy.get());
	}
	ConvolveInPlace(cs, set, rc);

	return rc;
}

void DeconvolutionManager::ConvolveInPlace(const CrossSection& cs, const EnergyDistributionSet& set, RateCoefficient& rc)
{
	rc.value.clear();
	rc.error.clear();
	for (const EnergyDistribution& eDist : set.distributions)
	{
		rc.value.push_back(0);
		rc.error.push_back(0);
		for (int i = 0; i < eDist.psi.size(); i++)
		{
			rc.value.back() += eDist.psi[i] * cs.values[i];
		}
	}
	rc.graph->Clear();
	for (int i = 0; i < rc.detuningEnergies.size(); i++)
	{
		rc.graph->AddPoint(rc.detuningEnergies.at(i), rc.value.at(i));
	}

	FileHandler::GetInstance().SaveRateCoefficients(rc);
}

double DeconvolutionManager::ConvolveFit(double* x, double* params)
{
	if (!x)
	{
		std::cout << "x is nullptr\n";
		return 0.0;
	}

	EnergyDistributionSet& set = energyDistributionSets.at(currentSetIndex);
	RateCoefficient& rc = rateCoefficientList.at(currentRateCoefficientIndex);

	double detuningEnergy = x[0];
	double sum = 0;
	
	// find correct distribution
	//std::cout << "Ed: " << detuningEnergy << "\n";
	int index = rc.GetIndexOfDetuningEnergy(detuningEnergy);
	if (index < 0)
	{
		std::cout << "did not find index of detuning energy " << detuningEnergy << std::endl;
		return 0.0;
	}

	EnergyDistribution& distribution = set.distributions.at(index);
	
	for (int i = 0; i < distribution.psi.size(); i++)
	{
		//std::cout << i << " ";
		//std::cout << distribution.GetPsis()[i] << " ";
		//std::cout << params[i] << "\n";
		sum += distribution.psi[i] * params[i] * params[i];
	}
	//std::cout << "sum " << sum << "\n";
	return sum;
}

void DeconvolutionManager::AddRateCoefficientToList(RateCoefficient& rc)
{
	rateCoefficientList.emplace_back(std::move(rc));
	currentRateCoefficientIndex = rateCoefficientList.size() - 1;
}

void DeconvolutionManager::RemoveRateCoefficient(int index)
{
	rateCoefficientList.erase(rateCoefficientList.begin() + index);
	currentRateCoefficientIndex = std::min(currentRateCoefficientIndex, (int)rateCoefficientList.size() - 1);
}

void DeconvolutionManager::AddCrossSectionToList(CrossSection& cs)
{
	crossSectionList.emplace_back(std::move(cs));
	currentCrossSectionIndex = crossSectionList.size() - 1;
}

void DeconvolutionManager::RemoveCrossSection(int index)
{
	crossSectionList.erase(crossSectionList.begin() + index);
	currentCrossSectionIndex = std::min(currentCrossSectionIndex, (int)crossSectionList.size() - 1);
}

void DeconvolutionManager::PlotRateCoefficient()
{
	m_mainCanvas->cd(1);

	int colors[5] = { kRed, kBlue, kGreen, kOrange, kMagenta };

	gPad->SetLogy();
	gPad->SetLogx();

	// Create a legend
	TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

	for (int i = 0; i < rateCoefficientList.size(); i++)
	{
		//if (!energyDistributionList[i]) return;

		rateCoefficientList[i].graph->SetLineColor(colors[i % 5]);
		rateCoefficientList[i].graph->SetMarkerStyle(21);

		rateCoefficientList[i].graph->GetXaxis()->SetTitle("E_d [eV]");
		rateCoefficientList[i].graph->GetYaxis()->SetTitle("alpha [m^3/s]");

		legend->AddEntry(rateCoefficientList[i].graph, rateCoefficientList[i].label.c_str(), "l");

		if (i == 0)
		{
			rateCoefficientList[i].graph->Draw("ALP");
		}
		else
		{
			rateCoefficientList[i].graph->Draw("ALP SAME");
		}
	}
	legend->Draw();
}

