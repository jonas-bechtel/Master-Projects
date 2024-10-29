#include "CrossSection.h"
#include "EnergyDistributionModel.h"
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

void CrossSection::CalculateRateCoefficients()
{
	EnergyDistributionModel* model = (EnergyDistributionModel*)Module::Get("Energy Distribution Model");
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
	EnergyDistributionModel* model = (EnergyDistributionModel*)Module::Get("Energy Distribution Model");
	std::vector<EnergyDistribution*>& energyDistributions = model->GetEnergyDistributions();

	for (EnergyDistribution* distribution : energyDistributions)
	{
		distribution->psi.clear();
		distribution->psi.resize(crossSection->GetNbinsX());

		for (double energy : distribution->collisionEnergies)
		{
			int bin = crossSection->FindBin(energy);
			double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
			distribution->psi[bin] += velocity;
		}
		for (int i = 0; i < distribution->psi.size(); i++)
		{
			distribution->psi[i] /= distribution->collisionEnergies.size();
		}
	}
}

void CrossSection::FitCrossSectionHistogram()
{
	CalculatePsis();

	TF1* fitFunction = new TF1("fit function", this, &CrossSection::FitFunction, 0, 99, crossSection->GetNbinsX());
	std::cout << fitFunction->GetNpar() << "\n";
	std::vector<double> initialGuess = std::vector<double>(crossSection->GetNbinsX(), 1);
	fitFunction->SetParameters(initialGuess.data());

	rateCoefficients->Fit(fitFunction, "R");

	double* parameter = fitFunction->GetParameters();

	// create Fit cross section
	crossSectionFit = (TH1D*)crossSection->Clone("cross section Fit");
	crossSectionFit->Reset();
	binValuesFit.clear();
	binCentersFit.clear();
	binValuesFit.reserve(crossSectionFit->GetNbinsX());
	binCentersFit.reserve(crossSectionFit->GetNbinsX());

	for (int i = 1; i <= crossSectionFit->GetNbinsX(); i++)
	{
		crossSectionFit->SetBinContent(i, parameter[i-1]);
		binValuesFit.push_back(parameter[i-1]);
		binCentersFit.push_back(crossSectionFit->GetBinCenter(i));
	}

	// create rate coefficient fit
	EnergyDistributionModel* model = (EnergyDistributionModel*)Module::Get("Energy Distribution Model");
	std::vector<EnergyDistribution*>& energyDistributions = model->GetEnergyDistributions();

	delete rateCoefficientsFit;
	rateCoefficientsFit = new TGraph();
	for (EnergyDistribution* eDist : energyDistributions)
	{
		double x[1] = {eDist->eDistParameter.detuningEnergy};
		rateCoefficientsFit->AddPoint(x[0], FitFunction(x, parameter));
	}
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

	for (int i = 0; i < crossSection->GetNbinsX(); i++)
	{
		//std::cout << i << " ";
		//std::cout << distribution->psi[i] << " ";
		//std::cout << params[i] << "\n";
		sum += distribution->psi[i] * params[i];
	}
	//std::cout << "sum " << sum << "\n";
	return sum;
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
	if (ImGui::Button("generate Cross section"))
	{
		GenerateCrossSection();
	}
	ImGui::SameLine();
	ImGui::SetNextItemWidth(100.0f);
	ImGui::InputInt("number bins", &nBins, 100);

	if (ImPlot::BeginPlot("cross section"))
	{
		ImPlot::SetupAxis(ImAxis_X1, "energy [eV]");
		ImPlot::SetupAxis(ImAxis_Y1, "sigma(E)", ImPlotAxisFlags_AutoFit);
		//ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);

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
	ImGui::SameLine();
	ImGui::Checkbox("use binned cross section", &useSigmaHist);
	if (ImGui::Button("fit cross section"))
	{
		FitCrossSectionHistogram();
	}
}
