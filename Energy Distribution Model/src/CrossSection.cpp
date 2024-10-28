#include "CrossSection.h"
#include "EnergyDistributionModel.h"
#include "PhysicalConstants.h"

CrossSection::CrossSection()
	: Module("Cross Section")
{
	GenerateCrossSection();
}

void CrossSection::GenerateCrossSection()
{
	if (crossSection) delete crossSection;

	crossSection = new TH1D("generated cross section", "generated cross section", nBins, 1e-8, 100);
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
	std::vector<EnergyDistribution> energyDistributions = model->GetEnergyDistributions();
	
	if (useSigmaHist)
	{
		delete rateCoefficients1; 
		rateCoefficients1 = new TGraph();

		for (auto& eDist : energyDistributions)
		{
			if (eDist.distribution && crossSection)
			{
				TH1D* temp = (TH1D*)eDist.distribution->Clone("temp");
				temp->Reset();

				for (int i = 1; i <= temp->GetNbinsX(); i++)
				{
					double collisionEnergyProbability = eDist.distribution->GetBinContent(i);
					double collisionEnergy = eDist.distribution->GetBinCenter(i);
					double collosionVelocity = TMath::Sqrt(2 * collisionEnergy * TMath::Qe() / PhysicalConstants::electronMass);
					double crossSectionValue = crossSection->Interpolate(collisionEnergy);
					//std::cout << "hist value: " << crossSectionValue << " theo: " << 1 / collisionEnergy << "\n";

					double value = collisionEnergyProbability * collosionVelocity * crossSectionValue;
					temp->SetBinContent(i, value);
				}
				double E_d = pow(sqrt(eDist.labEnergiesParameter.centerLabEnergy) - sqrt(eDist.eBeamParameter.coolingEnergy), 2);
				rateCoefficients1->AddPoint(E_d, temp->Integral());
				delete temp;
			}
		}
	}
	else
	{
		delete rateCoefficients2;
		rateCoefficients2 = new TGraph();

		for (auto& eDist : energyDistributions)
		{
			if (eDist.collisionEnergies.empty()) continue;

			for (double collisionEnergy : eDist.collisionEnergies)
			{
				double crossSectionValue = 1 / collisionEnergy;
				double collosionVelocity = TMath::Sqrt(2 * collisionEnergy * TMath::Qe() / PhysicalConstants::electronMass);
				eDist.rateCoefficient += crossSectionValue * collosionVelocity;
			}
			eDist.rateCoefficient /= eDist.collisionEnergies.size();

			double E_d = pow(sqrt(eDist.labEnergiesParameter.centerLabEnergy) - sqrt(eDist.eBeamParameter.coolingEnergy), 2);
			rateCoefficients2->AddPoint(E_d, eDist.rateCoefficient);
		}
	}
}

void CrossSection::PlotRateCoefficients()
{
	m_mainCanvas->cd(1);

	std::cout << rateCoefficients1->GetN() << "\n";
	std::cout << rateCoefficients2->GetN() << "\n";

	if (rateCoefficients1->GetN())
	{
		
		rateCoefficients1->SetLineColor(kBlue);
		rateCoefficients2->SetMarkerStyle(21);
		rateCoefficients1->SetTitle("rate Coefficients");
		rateCoefficients1->GetXaxis()->SetTitle("E_d [eV]");
		rateCoefficients1->GetYaxis()->SetTitle("alpha [m^3/s]");
		rateCoefficients1->Draw("ALP");
	}
	if (rateCoefficients2->GetN())
	{
		rateCoefficients2->SetLineColor(kRed);
		rateCoefficients2->SetMarkerStyle(21);
		rateCoefficients2->Draw("ALP SAME");
	}
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
		ImPlot::PlotLine("rate coefficient1", rateCoefficients1->GetX(), rateCoefficients1->GetY(), rateCoefficients1->GetN());
		ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
		ImPlot::PlotLine("rate coefficient2", rateCoefficients2->GetX(), rateCoefficients2->GetY(), rateCoefficients2->GetN());
		
		ImPlot::EndPlot();
	}

	if (ImGui::Button("calculate rate Coefficients"))
	{
		CalculateRateCoefficients();
		PlotRateCoefficients();
	}
	ImGui::SameLine();
	ImGui::Checkbox("use binned cross section", &useSigmaHist);
}
