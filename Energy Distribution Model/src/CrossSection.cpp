#include "pch.h"

#include <TDecompSVD.h>

#include "CrossSection.h"
#include "EnergyDistributionManager.h"
#include "Constants.h"


CrossSection::CrossSection()
	: Module("Cross Section")
{
	SetupTrueCrossSection();
}

void CrossSection::test()
{
	SetupFitCrossSectionHist();

	CalculatePsis();

	binValuesFit.clear();
	binCentersFit.clear();
	binValuesFit.reserve(crossSectionFit->GetNbinsX());
	binCentersFit.reserve(crossSectionFit->GetNbinsX());

	std::vector<double> params;
	params.reserve(crossSectionFit->GetNbinsX());

	for (int i = 1; i <= crossSectionFit->GetNbinsX(); i++)
	{
		params.push_back(1 / crossSectionFit->GetBinCenter(i));
	}

	FillFitPlots(params.data());
}

void CrossSection::SetupFitCrossSectionHist()
{
	EnergyDistributionManager* model = (EnergyDistributionManager*)Module::Get("Energy Distribution Manager");
	std::vector<EnergyDistribution*>& energyDistributions = model->GetEnergyDistributions();

	std::vector<double> binEdges;
	double maxEnergy = 100; //just a gues, not fixed
	double minEnergy = energyDistributions.back()->eBeamParameter.detuningEnergy / 10;
	
	// binning like in the paper 
	if (currentOption == PaperBinning)
	{
		EnergyDistribution* representativeEnergyDist = energyDistributions[energyDistributions.size() - 1];
		double kT_trans = representativeEnergyDist->eBeamParameter.transverse_kT;
		double kT_long = 0.00047; // representativeEnergyDist->eBeamParameter.longitudinal_kT; //eBeam->GetLongitudinal_kT(labEnergiesParameter.centerLabEnergy);

		double factor = 1;

		binEdges.push_back(0);
		binEdges.push_back(factor * kT_trans / 20);

		for (int i = 1; binEdges[i] < kT_trans; i++)
		{
			binEdges.push_back(factor * 2 * binEdges[i]);
			//std::cout << binEdges[i + 1] << "\n";
		}
		while (binEdges[binEdges.size() - 1] < maxEnergy)
		{
			double previousEdge = binEdges[binEdges.size() - 1];
			double delta_E = sqrt(pow((kT_trans * log(2)), 2) + 16 * log(2) * kT_long * previousEdge);
			binEdges.push_back(previousEdge + factor * delta_E);

			//std::cout << previousEdge << " " << delta_E << "\n";
			//std::cout << binEdges[binEdges.size()] << "\n";
			//std::cout << (binEdges[binEdges.size()] < maxEnergy) << "\n";
		}
	}
	// constant binning
	if (currentOption == ConstantBinning)
	{
		binEdges.reserve(numberBins + 1);
		binEdges.push_back(0);
		double binWidth = maxEnergy / numberBins;
		for (int i = 0; i < numberBins; i++)
		{
			binEdges.push_back(binWidth * (i + 1));
		}
	}
	// bin width increses by a constant factor
	if (currentOption == FactorBinning)
	{
		if (limitBinSize)
		{
			double min = minEnergy; // std::max(energyRange[0], 1e-9f);
			double factor = TMath::Power((maxEnergy / min), (1.0 / numberBins));

			binEdges.reserve(numberBins + 1);
			binEdges.push_back(0);
			for (int i = 0; binEdges[binEdges.size() - 1] < maxEnergy; i++)
			{
				double nextBin = binEdges[i] * factor;
				double difference = std::max(nextBin - binEdges[i], minBinSize);
				binEdges.push_back(binEdges[i] + difference);
			}
		}
		else
		{
			float min = minEnergy;
			double factor = TMath::Power((maxEnergy / min), (1.0 / numberBins));

			binEdges.reserve(numberBins + 1);
			binEdges.push_back(min);
			for (int i = 0; i < numberBins; i++)
			{
				binEdges.push_back(binEdges[i] * factor);
			}
		}
		
	}
	if (currentOption == PaperFactorMix)
	{
		EnergyDistribution* representativeEnergyDist = energyDistributions[energyDistributions.size() - 1];
		double kT_trans = representativeEnergyDist->eBeamParameter.transverse_kT;
		//double kT_long = 0.00047; // representativeEnergyDist->eBeamParameter.longitudinal_kT; //eBeam->GetLongitudinal_kT(labEnergiesParameter.centerLabEnergy);

		binEdges.push_back(0);
		binEdges.push_back(kT_trans / 20);

		for (int i = 1; binEdges[i] < kT_trans; i++)
		{
			binEdges.push_back(2 * binEdges[i]);
			//std::cout << binEdges[i + 1] << "\n";
		}
		
		double factor = TMath::Power((maxEnergy / binEdges.back()), (1.0 / numberBins));

		for (int i = 0; i < numberBins; i++)
		{
			binEdges.push_back(binEdges.back() * factor);
		}
	}
	//for (double edge : binEdges)
	//{
	//	std::cout << edge << "\n";
	//}
	std::cout << "number cross section bins: " << binEdges.size() - 1 << "\n";

	crossSectionFit = new TH1D("cross section fit", "cross section fit", binEdges.size() - 1, binEdges.data());
}

void CrossSection::CalculateRateCoefficients()
{
	EnergyDistributionManager* model = (EnergyDistributionManager*)Module::Get("Energy Distribution Manager");
	std::vector<EnergyDistribution*>& energyDistributions = model->GetEnergyDistributions();
	
	//if (useSigmaHist)
	//{
	//	delete rateCoefficients1; 
	//	rateCoefficients1 = new TGraph();
	//
	//	for (auto& eDist : energyDistributions)
	//	{
	//		if (eDist->distribution && crossSection)
	//		{
	//			TH1D* temp = (TH1D*)eDist->distribution->Clone("temp");
	//			temp->Reset();
	//
	//			for (int i = 1; i <= temp->GetNbinsX(); i++)
	//			{
	//				double collisionEnergyProbability = eDist->distribution->GetBinContent(i);
	//				double collisionEnergy = eDist->distribution->GetBinCenter(i);
	//				double collosionVelocity = TMath::Sqrt(2 * collisionEnergy * TMath::Qe() / PhysicalConstants::electronMass);
	//				double crossSectionValue = crossSection->Interpolate(collisionEnergy);
	//				//std::cout << "hist value: " << crossSectionValue << " theo: " << 1 / collisionEnergy << "\n";
	//
	//				double value = collisionEnergyProbability * collosionVelocity * crossSectionValue;
	//				temp->SetBinContent(i, value);
	//			}
	//			double E_d = pow(sqrt(eDist->labEnergiesParameter.centerLabEnergy) - sqrt(eDist->eBeamParameter.coolingEnergy), 2);
	//			rateCoefficients1->AddPoint(E_d, temp->Integral());
	//			delete temp;
	//		}
	//	}
	//}
	//else
	{
		delete rateCoefficients;
		rateCoefficients = new TGraph();

		for (EnergyDistribution* eDist : energyDistributions)
		{
			if (eDist->collisionEnergies.empty()) continue;

			for (double collisionEnergy : eDist->collisionEnergies)
			{
				double crossSectionValue = 1 / collisionEnergy;
				double collosionVelocity = TMath::Sqrt(2 * collisionEnergy * TMath::Qe() / PhysicalConstants::electronMass);
				eDist->rateCoefficient += crossSectionValue * collosionVelocity;
			}
			eDist->rateCoefficient /= eDist->collisionEnergies.size();

			double E_d = eDist->eBeamParameter.detuningEnergy; //pow(sqrt(eDist->labEnergiesParameter.centerLabEnergy) - sqrt(eDist->eBeamParameter.coolingEnergy), 2);
			rateCoefficients->AddPoint(E_d, eDist->rateCoefficient);
		}
	}
}

void CrossSection::CalculatePsis()
{
	EnergyDistributionManager* model = (EnergyDistributionManager*)Module::Get("Energy Distribution Manager");
	std::vector<EnergyDistribution*>& energyDistributions = model->GetEnergyDistributions();

	for (EnergyDistribution* distribution : energyDistributions)
	{
		distribution->psi.clear();
		int nBins = crossSectionFit->GetNbinsX();
		std::cout << nBins << std::endl;
		//distribution->psi.reserve(nBins);
		distribution->psi.resize(nBins);

		std::cout << "distribution: " << distribution->index << std::endl;
		for (double energy : distribution->collisionEnergies)
		{
			int bin = crossSectionFit->FindBin(energy);
			if (bin == 0)
			{
				std::cout << "crossection hist too small, got underflow when putting in " << energy << std::endl;
				continue;
			}
			//if (bin < 10) std::cout << "bin " << bin << ": " << energy << std::endl;
			double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
			if (bin - 1 >= distribution->psi.size())
			{
				std::cout << "want to access " << bin - 1 << " but size is " << distribution->psi.size() << std::endl;
				std::cout << energy << std::endl;
			}
			distribution->psi[bin - 1] += velocity;
		}
		for (int i = 0; i < distribution->psi.size(); i++)
		{
			distribution->psi[i] /= distribution->collisionEnergies.size();
			std::cout << "Psi_" << i << ": " << distribution->psi[i] << "\t" << crossSectionFit->GetBinLowEdge(i+1)
				 << " - " << crossSectionFit->GetBinLowEdge(i+2) << "\n";
		}
	}
}

void CrossSection::SetupInitialGuess()
{
	//for (int i = 1; i <= crossSectionFit->GetNbinsX(); i++)
	//{
	//	initialGuess.push_back(2 / crossSectionFit->GetBinCenter(i));
	//}
	initialGuess.clear();

	for (int i = 1; i <= crossSectionFit->GetNbinsX(); i++)
	{
		double energy = crossSectionFit->GetBinCenter(i);
		double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
		double alpha = rateCoefficients->Eval(energy);
		initialGuess.push_back(alpha / velocity);
	}
}

void CrossSection::FitCrossSectionHistogram()
{
	EnergyDistributionManager* model = (EnergyDistributionManager*)Module::Get("Energy Distribution Manager");
	if (model->GetEnergyDistributions().empty())
	{
		std::cout << "no energy distributions\n";
		return;
	}

	if (initialGuess.empty())
	{
		// create Fit cross section
		SetupFitCrossSectionHist();
		CalculatePsis();
		SetupInitialGuess();
	}
	else
	{
		initialGuess.clear();
		initialGuess = binValuesFit;
		for (double& value : initialGuess) 
		{
			value = std::abs(value);
		}
	}
	TF1* fitFunction = new TF1("fit function", this, &CrossSection::FitFunction, 0, 99, crossSectionFit->GetNbinsX());

	fitFunction->SetParameters(initialGuess.data());
	for (int i = fixParamStart; i < fixParamStop; i++)
	{
		std::cout << i << " " << initialGuess[i] << std::endl;
		fitFunction->FixParameter(i, initialGuess[i]);
	}

	binValuesFit.clear();
	binCentersFit.clear();
	binValuesFit.reserve(crossSectionFit->GetNbinsX());
	binCentersFit.reserve(crossSectionFit->GetNbinsX());
	
	rateCoefficients->Fit(fitFunction, "RN");

	double* parameter = fitFunction->GetParameters();
	FillFitPlots(parameter);

	fitFunction->Delete();
}

void CrossSection::FitWithSVD()
{
	// solve Ax = b with SVD

	EnergyDistributionManager* model = (EnergyDistributionManager*)Module::Get("Energy Distribution Manager");
	std::vector<EnergyDistribution*>& energyDistributions = model->GetEnergyDistributions();

	if (model->GetEnergyDistributions().empty())
	{
		std::cout << "no energy distributions\n";
		return;
	}

	SetupFitCrossSectionHist();
	CalculatePsis();

	binValuesFit.clear();
	binCentersFit.clear();
	binValuesFit.reserve(crossSectionFit->GetNbinsX());
	binCentersFit.reserve(crossSectionFit->GetNbinsX());

	int n = energyDistributions.size();
	int p = crossSectionFit->GetNbinsX();

	// matrix A is all the p Psis for all the n distributions, n >= p is required so the rest is filled with 0
	TMatrixD PsiMatrix(std::max(n, p), p);
	// vector b with all the rate coeffictions
	TVectorD alphaVector(std::max(n, p));

	// fill matrix and vector
	for (int i = 0; i < std::max(n, p); i++)
	{
		
		// all fitted data points
		for (int j = 0; j < p; j++)
		{
			if (i >= n)
			{
				PsiMatrix[i][j] = 0;
			}
			else
			{
				PsiMatrix[i][j] = energyDistributions[i]->psi[j];
			}			
		}
		if (i >= n)
		{
			alphaVector[i] = 0;
		}
		else
		{
			alphaVector[i] = energyDistributions[i]->rateCoefficient;
		}
	}

	alphaVector.Print();

	TDecompSVD decomp(PsiMatrix);
	// solve for x, result is put into input vector
	decomp.Solve(alphaVector);

	alphaVector.Print();
	FillFitPlots(alphaVector.GetMatrixArray());
}

double CrossSection::FitFunction(double* x, double* params)
{
	if (!x) std::cout << "x is nullptr\n";

	double detuningEnergy = x[0];
	double sum = 0;

	// find correct distribution
	//std::cout << "Ed: " << detuningEnergy << "\n";
	EnergyDistribution* distribution = EnergyDistribution::FindByEd(detuningEnergy);
	if (!distribution) return 0.0;

	for (int i = 0; i < crossSectionFit->GetNbinsX(); i++)
	{
		//std::cout << i << " ";
		//std::cout << distribution->GetPsis()[i] << " ";
		//std::cout << params[i] << "\n";
		sum += distribution->psi[i] * params[i];
	}
	//std::cout << "sum " << sum << "\n";
	return sum;
}

void CrossSection::FillFitPlots(double* crossSectionParamater)
{
	for (int i = 1; i <= crossSectionFit->GetNbinsX(); i++)
	{
		crossSectionFit->SetBinContent(i, crossSectionParamater[i - 1]);
		binValuesFit.push_back(crossSectionParamater[i - 1]);
		binCentersFit.push_back(crossSectionFit->GetBinCenter(i));
	}

	// create rate coefficient fit
	EnergyDistributionManager* model = (EnergyDistributionManager*)Module::Get("Energy Distribution Manager");
	std::vector<EnergyDistribution*>& energyDistributions = model->GetEnergyDistributions();

	delete rateCoefficientsFit;
	rateCoefficientsFit = new TGraph();
	for (EnergyDistribution* eDist : energyDistributions)
	{
		double x[1] = { eDist->eBeamParameter.detuningEnergy};
		rateCoefficientsFit->AddPoint(x[0], FitFunction(x, crossSectionParamater));
	}
}

void CrossSection::PlotRateCoefficients()
{
	m_mainCanvas->cd(1);

	//std::cout << rateCoefficients1->GetN() << "\n";
	//std::cout << rateCoefficients2->GetN() << "\n";

	if (rateCoefficients->GetN())
	{
		
		rateCoefficients->SetLineColor(kBlue);
		rateCoefficients->SetMarkerStyle(21);
		rateCoefficients->SetTitle("rate Coefficients");
		rateCoefficients->GetXaxis()->SetTitle("E_d [eV]");
		rateCoefficients->GetYaxis()->SetTitle("alpha [m^3/s]");
		rateCoefficients->Draw("ALP");
	}
	//if (rateCoefficients2->GetN())
	//{
	//	rateCoefficients2->SetLineColor(kRed);
	//	rateCoefficients2->SetMarkerStyle(21);
	//	rateCoefficients2->Draw("ALP SAME");
	//}
}

void CrossSection::ShowUI()
{
	//if (ImGui::Button("generate Cross section"))
	//{
	//	GenerateCrossSection();
	//}
	//ImGui::SameLine();
	//ImGui::SetNextItemWidth(100.0f);
	//ImGui::InputInt("number bins", &nBins, 100);

	if (ImPlot::BeginPlot("cross section"))
	{
		ImPlot::SetupAxis(ImAxis_X1, "collision energy [eV]");
		ImPlot::SetupAxis(ImAxis_Y1, "sigma(E)");
		if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
		if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

		ImPlot::SetupLegend(ImPlotLocation_NorthEast);

		ImPlot::PlotLine("cross section", binCentersTrue.data(), binValuesTrue.data(), binCentersTrue.size());
		ImPlot::PlotLine("initial guess", binCentersFit.data(), initialGuess.data(), initialGuess.size());
		ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, 3.0);
		ImPlot::PlotLine("cross section Fit", binCentersFit.data(), binValuesFit.data(), binCentersFit.size());

		ImPlot::EndPlot();
	}
	ImGui::Checkbox("log X", &logX);
	ImGui::SameLine();
	ImGui::Checkbox("log Y", &logY);
	if (ImPlot::BeginPlot("rate coefficient"))
	{
		ImPlot::SetupAxis(ImAxis_X1, "detuning energy [eV]");
		ImPlot::SetupAxis(ImAxis_Y1, "rate coefficient");
		if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
		if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

		ImPlot::SetupLegend(ImPlotLocation_NorthEast);

		ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
		ImPlot::PlotLine("rate coefficient", rateCoefficients->GetX(), rateCoefficients->GetY(), rateCoefficients->GetN());
		ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
		ImPlot::PlotLine("rate coefficient fit", rateCoefficientsFit->GetX(), rateCoefficientsFit->GetY(), rateCoefficientsFit->GetN());
		
		ImPlot::EndPlot();
	}

	if (ImGui::Button("calculate rate Coefficients"))
	{
		CalculateRateCoefficients();
		PlotRateCoefficients();
	}

	ImGui::SetNextItemWidth(150.0f);
	ImGui::Combo("binning options", &currentOption, binningOptions, IM_ARRAYSIZE(binningOptions));
	
	if (currentOption == ConstantBinning || currentOption == FactorBinning || currentOption == PaperFactorMix)
	{
		ImGui::SameLine();
		ImGui::SetNextItemWidth(100.0f);
		ImGui::InputInt("number bins", &numberBins);
	}
	if (currentOption == FactorBinning)
	{
		ImGui::SameLine();
		ImGui::Checkbox("limit bin size", &limitBinSize);
		ImGui::SameLine();
		ImGui::BeginDisabled(!limitBinSize);
		ImGui::SetNextItemWidth(100.0f);
		ImGui::InputDouble("min bin size", &minBinSize, 0.0, 0.0, "%.1e");
		ImGui::EndDisabled();
	}
	if (ImGui::Button("fit cross section"))
	{
		FitCrossSectionHistogram();
	}
	ImGui::SameLine();
	if (ImGui::Button("clear fit data"))
	{
		initialGuess.clear();
		binCentersFit.clear();
		binValuesFit.clear();
		crossSectionFit->Clear();
	}
	ImGui::SameLine();
	if (ImGui::Button("test"))
	{
		test();
	}
	ImGui::SameLine();
	if (ImGui::InputInt2("fix parameter", &fixParamStart))
	{
		fixParamStop = std::min(fixParamStop, (int)initialGuess.size());
	}
	if (ImGui::Button("Fit with SVD"))
	{
		FitWithSVD();
	}
}

void CrossSection::SetupTrueCrossSection()
{
	double min = 1e-4;
	double max = 100;
	int numberPoints = 100;
	double step = (max - min) / numberPoints;
	for (int i = 0; i <= 1000; i++)
	{
		double energy = min + i * step;
		binCentersTrue.push_back(energy);
		binValuesTrue.push_back(1 / energy);
	}
}
