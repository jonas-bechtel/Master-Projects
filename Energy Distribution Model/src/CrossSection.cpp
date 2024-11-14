#include "CrossSection.h"
#include "EnergyDistributionManager.h"
#include "PhysicalConstants.h"

#include <TF1.h>

CrossSection::CrossSection()
	: Module("Cross Section")
{
	GenerateCrossSection();
}

void CrossSection::GenerateCrossSection()
{
	if (crossSection) delete crossSection;

	crossSection = new TH1D("generated cross section", "generated cross section", nBins, 1e-5, 100);
	binCenters.clear();
	binValues.clear();
	binCenters.reserve(nBins);
	binValues.reserve(nBins);
	for (int i = 1; i <= nBins; i++)
	{
		double energyValue = crossSection->GetBinCenter(i);
		binCenters.push_back(energyValue);
		double crossSectionValue = 1 / energyValue;
		binValues.push_back(crossSectionValue);
		crossSection->SetBinContent(i, crossSectionValue);
	}
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
	
	// binning like in the paper 
	if (currentOption == 0)
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
	if (currentOption == 1)
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
	if (currentOption == 2)
	{
		if (limitBinSize)
		{
			double min = 1e-6; // std::max(energyRange[0], 1e-9f);
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
			float min = 1e-6f;
			double factor = TMath::Power((maxEnergy / min), (1.0 / numberBins));

			binEdges.reserve(numberBins + 1);
			binEdges.push_back(min);
			for (int i = 0; i < numberBins; i++)
			{
				binEdges.push_back(binEdges[i] * factor);
			}
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

			double E_d = pow(sqrt(eDist->labEnergiesParameter.centerLabEnergy) - sqrt(eDist->eBeamParameter.coolingEnergy), 2);
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
		distribution->psi.resize(crossSectionFit->GetNbinsX());

		for (double energy : distribution->collisionEnergies)
		{
			int bin = crossSectionFit->FindBin(energy);
			double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
			distribution->psi[bin - 1] += velocity;
		}
		for (int i = 0; i < distribution->psi.size(); i++)
		{
			distribution->psi[i] /= distribution->collisionEnergies.size();
			//std::cout << "GetPsis()_" << i << ": " << distribution->GetPsis()[i] << " " << crossSectionFit->GetBinLowEdge(i+1)
			//	 << " - " << crossSectionFit->GetBinLowEdge(i+2) << "\n";
		}
	}
}

void CrossSection::FitCrossSectionHistogram()
{
	// create Fit cross section
	SetupFitCrossSectionHist();

	CalculatePsis();

	binValuesFit.clear();
	binCentersFit.clear();
	binValuesFit.reserve(crossSectionFit->GetNbinsX());
	binCentersFit.reserve(crossSectionFit->GetNbinsX());

	TF1* fitFunction = new TF1("fit function", this, &CrossSection::FitFunction, 0, 99, crossSectionFit->GetNbinsX());
	//std::cout << fitFunction->GetNpar() << "\n";

	// just a random guess of all 1s
	std::vector<double> initialGuess; // = std::vector<double>(crossSectionFit->GetNbinsX(), 0.1);
	for (int i = 1; i <= crossSectionFit->GetNbinsX(); i++)
	{
		initialGuess.push_back(2 / crossSectionFit->GetBinCenter(i));
	}
	fitFunction->SetParameters(initialGuess.data());

	rateCoefficients->Fit(fitFunction, "R");

	double* parameter = fitFunction->GetParameters();
	FillFitPlots(parameter);
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
		ImPlot::SetupAxis(ImAxis_X1, "energy [eV]");
		ImPlot::SetupAxis(ImAxis_Y1, "sigma(E)");
		if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
		if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

		ImPlot::SetupLegend(ImPlotLocation_NorthEast);

		ImPlot::PlotLine("cross section", binCenters.data(), binValues.data(), binCenters.size());
		ImPlot::PlotLine("cross section Fit", binCentersFit.data(), binValuesFit.data(), binCentersFit.size());

		ImPlot::EndPlot();
	}
	ImGui::Checkbox("log X", &logX);
	ImGui::SameLine();
	ImGui::Checkbox("log Y", &logY);
	if (ImPlot::BeginPlot("rate coefficient"))
	{
		ImPlot::SetupAxis(ImAxis_X1, "energy [eV]");
		ImPlot::SetupAxis(ImAxis_Y1, "rate coefficient");
		if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
		if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

		ImPlot::SetupLegend(ImPlotLocation_NorthEast);

		ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
		ImPlot::PlotLine("rate coefficient", rateCoefficients->GetX(), rateCoefficients->GetY(), rateCoefficients->GetN());
		ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
		ImPlot::PlotLine("rate coefficient2", rateCoefficientsFit->GetX(), rateCoefficientsFit->GetY(), rateCoefficientsFit->GetN());
		
		ImPlot::EndPlot();
	}

	if (ImGui::Button("calculate rate Coefficients"))
	{
		CalculateRateCoefficients();
		PlotRateCoefficients();
	}

	ImGui::SetNextItemWidth(150.0f);
	ImGui::Combo("binning options", &currentOption, binningOptions, IM_ARRAYSIZE(binningOptions));
	
	if (currentOption == 1 || currentOption == 2)
	{
		ImGui::SameLine();
		ImGui::SetNextItemWidth(100.0f);
		ImGui::InputInt("number bins", &numberBins);
	}
	if (currentOption == 2)
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
	if (ImGui::Button("test"))
	{
		test();
	}
}
