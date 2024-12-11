#include "pch.h"
#include "CrossSectionManager.h"


CrossSectionManager::CrossSectionManager()
	: CrossSectionDeconvolutionModule("Cross Section")
{
}


//double CrossSectionManager::FitFunction(double* x, double* params)
//{
//	if (!x) std::cout << "x is nullptr\n";
//
//	double detuningEnergy = x[0];
//	double sum = 0;
//
//	// find correct distribution
//	//std::cout << "Ed: " << detuningEnergy << "\n";
//	EnergyDistribution* distribution = EnergyDistribution::FindByEd(detuningEnergy);
//	if (!distribution) return 0.0;
//
//	for (int i = 0; i < crossSectionFit->GetNbinsX(); i++)
//	{
//		//std::cout << i << " ";
//		//std::cout << distribution.GetPsis()[i] << " ";
//		//std::cout << params[i] << "\n";
//		sum += distribution->psi[i] * params[i];
//	}
//	//std::cout << "sum " << sum << "\n";
//	return sum;
//}

//void CrossSectionManager::FillFitPlots(double* crossSectionParamater)
//{
//	for (int i = 1; i <= crossSectionFit->GetNbinsX(); i++)
//	{
//		crossSectionFit->SetBinContent(i, crossSectionParamater[i - 1]);
//		binValuesFit.push_back(crossSectionParamater[i - 1]);
//		binCentersFit.push_back(crossSectionFit->GetBinCenter(i));
//	}
//
//	// create rate coefficient fit
//	EnergyDistributionManager* model = (EnergyDistributionManager*)Window::Get("Energy Distribution Manager");
//	std::vector<EnergyDistribution>& energyDistributionList = model->GetEnergyDistributions();
//
//	delete rateCoefficientsFit;
//	rateCoefficientsFit = new TGraph();
//	for (EnergyDistribution& eDist : energyDistributionList)
//	{
//		double x[1] = { eDist.eBeamParameter.detuningEnergy };
//		rateCoefficientsFit->AddPoint(x[0], FitFunction(x, crossSectionParamater));
//	}
//}

void CrossSectionManager::ShowUI()
{
	if (ImPlot::BeginPlot("cross section"))
	{
		ImPlot::SetupAxis(ImAxis_X1, "collision energy [eV]");
		ImPlot::SetupAxis(ImAxis_Y1, "sigma(E)"); 
		if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
		if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
		ImPlot::SetupLegend(ImPlotLocation_NorthEast);

		for (const CrossSection& cs : crossSectionList)
		{
			ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
			ImPlot::PlotErrorBars("##errors", cs.energy.data(), cs.value.data(), cs.error.data(), cs.value.size());
			ImPlot::PlotLine(cs.label.c_str(), cs.energy.data(), cs.value.data(), cs.value.size());
		}

		ImPlot::EndPlot();
	}
	ImGui::Checkbox("log X", &logX);
	ImGui::SameLine();
	ImGui::Checkbox("log Y", &logY);


	//if (ImGui::Button("calculate rate Coefficients"))
	//{
	//	CalculateRateCoefficients();
	//	PlotRateCoefficients();
	//}

	//ImGui::SetNextItemWidth(150.0f);
	//ImGui::Combo("binning options", &currentOption, binningOptions, IM_ARRAYSIZE(binningOptions));
	//
	//if (currentOption == ConstantBinning || currentOption == FactorBinning || currentOption == PaperFactorMix)
	//{
	//	ImGui::SameLine();
	//	ImGui::SetNextItemWidth(100.0f);
	//	ImGui::InputInt("number bins", &numberBins);
	//}
	//if (currentOption == FactorBinning)
	//{
	//	ImGui::SameLine();
	//	ImGui::Checkbox("limit bin size", &limitBinSize);
	//	ImGui::SameLine();
	//	ImGui::BeginDisabled(!limitBinSize);
	//	ImGui::SetNextItemWidth(100.0f);
	//	ImGui::InputDouble("min bin size", &minBinSize, 0.0, 0.0, "%.1e");
	//	ImGui::EndDisabled();
	//}
	//if (currentOption == PaperBinning)
	//{
	//	ImGui::SameLine();
	//	ImGui::InputDouble("factor", &binFactor);
	//}
	//if (ImGui::Button("fit cross section"))
	//{
	//	FitCrossSectionHistogram();
	//}
	//ImGui::SameLine();
	//if (ImGui::Button("clear fit data"))
	//{
	//	initialGuess.clear();
	//	binCentersFit.clear();
	//	binValuesFit.clear();
	//	crossSectionFit->Clear();
	//}
	//ImGui::SameLine();
	//if (ImGui::Button("test"))
	//{
	//	test();
	//}
	//ImGui::SameLine();
	//if (ImGui::InputInt2("fix parameter", &fixParamStart))
	//{
	//	fixParamStop = std::min(fixParamStop, (int)initialGuess.size());
	//}
	//ImGui::SameLine();
	//ImGui::Checkbox("limit params", &limitParamRange);
	//if (ImGui::Button("Fit with SVD"))
	//{
	//	FitWithSVD();
	//}
	//ImGui::SameLine();
	//if (ImGui::Button("Fit with Eigen SVD"))
	//{
	//	FitWithEigenSVD();
	//}
	//ImGui::SameLine();
	//if (ImGui::Button("Fit with GD"))
	//{
	//	FitWithEigenGD();
	//}
	//ImGui::SameLine();
	//ImGui::InputDouble("lr", &learningRate, 0, 0, "%e");
	//ImGui::SameLine();
	//ImGui::InputDouble("lambda", &lambda);
	//ImGui::SameLine();
	//ImGui::InputInt("iter", &iterations);
	//
	//if (ImGui::Button("Fit with Torch"))
	//{
	//	FitWithTorch();
	//}
	//ImGui::SameLine();
	//ImGui::InputDouble("lr ##torch", &torchLearningRate, 0, 0, "%e");
	//ImGui::SameLine();
	//ImGui::InputInt("epochs", &nEpochs);
	//ImGui::SameLine();
	//ImGui::InputDouble("smooth ##torch", &smoothRegularisation);
	//ImGui::SameLine();
	//ImGui::InputDouble("l2 ##torch", &l2regularisation);
}
