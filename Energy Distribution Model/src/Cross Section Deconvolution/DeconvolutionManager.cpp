#include "pch.h"
#include "DeconvolutionManager.h"
#include "RateCoefficient.h"
#include "FileHandler.h"
#include "CrossSectionManager.h"

DeconvolutionManager::DeconvolutionManager()
	: CrossSectionDeconvolutionModule("rate coefficient")
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
				rc.input = true;
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
				ImGui::PopID();

			}
			ImGui::EndListBox();
		}
		ImGui::End();
	}
}

void DeconvolutionManager::ShowUI()
{
	ImGui::SeparatorText("rate coefficient fitting inputs");
	ImGui::Text("energy distribution set: %s", energyDistributionSets.at(currentSetIndex).Label().c_str());
	ImGui::Text("target rate coefficient: %s", 
		rateCoefficientList.empty() ? "" : rateCoefficientList.at(currentRateCoefficientIndex).label.c_str());
	ImGui::Text("cross section binning: %s", binningOptions[currentSettings.scheme]);
	if (ImGui::Button("Deconvolve Cross Section"))
	{
		EnergyDistributionSet& currentSet = energyDistributionSets.at(currentSetIndex);
		RateCoefficient& currentRC = rateCoefficientList.at(currentRateCoefficientIndex);
		CrossSection cs = Deconvolve(currentRC, currentSet);
		std::cout << "cs energy size: " << cs.energies.size() << std::endl;

		RateCoefficient rc = Convolve(cs, currentSet);
		AddCrossSectionToList(cs);
		AddRateCoefficientToList(rc);
	}

	if (ImPlot::BeginPlot("rate coefficient"))
	{
		ImPlot::SetupAxis(ImAxis_X1, "detuning energy [eV]");
		ImPlot::SetupAxis(ImAxis_Y1, "rate coefficient");
		if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
		if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
		ImPlot::SetupLegend(ImPlotLocation_NorthEast);

		for (const RateCoefficient& rc : rateCoefficientList)
		{
			ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
			ImPlot::PlotErrorBars("##errors", rc.detuningEnergies.data(), rc.value.data(), rc.error.data(), rc.error.size());
			ImPlot::PlotLine(rc.label.c_str(), rc.detuningEnergies.data(), rc.value.data(), rc.value.size());
		}
		
		ImPlot::EndPlot();
	}

	ImGui::Checkbox("log X", &logX);
	ImGui::SameLine();
	ImGui::Checkbox("log Y", &logY);

	if (ImPlot::BeginPlot("cross section"))
	{
		ImPlot::SetupAxis(ImAxis_X1, "collision energy [eV]");
		ImPlot::SetupAxis(ImAxis_Y1, "sigma(E)");
		if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
		if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
		ImPlot::SetupLegend(ImPlotLocation_NorthEast);

		for (const CrossSection& cs : crossSectionList)
		{
			//ImGui::PushID();
			ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
			ImPlot::PlotErrorBars("##errors", cs.energies.data(), cs.values.data(), cs.errors.data(), cs.errors.size());
			ImPlot::PlotLine(cs.label.c_str(), cs.energies.data(), cs.values.data(), cs.values.size());
			ImPlot::PlotLine((cs.label + "init").c_str(), cs.energies.data(), cs.initialGuess.data(), cs.initialGuess.size());
			//ImGui::PopID();
		}

		ImPlot::EndPlot();
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

	cs.SetupBinning(currentSettings);

	set.CalculatePsisFromBinning(cs.hist);
	cs.SetupInitialGuess(rc);
	DeconvolveInPlace(rc, set, cs);

	return cs;
}

void DeconvolutionManager::DeconvolveInPlace(const RateCoefficient& rc, const EnergyDistributionSet& set, CrossSection& cs)
{
	TF1* fitFunction = new TF1("fit function", this, &DeconvolutionManager::ConvolveFit, 0, 100, cs.hist->GetNbinsX());

	fitFunction->SetParameters(cs.values.data());

	rc.graph->Fit(fitFunction, "RN");

	double* parameter = fitFunction->GetParameters();
	cs.SetValues(parameter);
}

RateCoefficient DeconvolutionManager::Convolve(CrossSection& cs, EnergyDistributionSet& set)
{
	RateCoefficient rc = RateCoefficient();
	rc.label = "bla";
	rc.input = false;

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
	for (const EnergyDistribution& eDist : set.distributions)
	{
		rc.value.push_back(0);
		for (int i = 0; i < eDist.psi.size(); i++)
		{
			rc.value.back() += eDist.psi[i] * cs.values[i];
		}
	}
}

double DeconvolutionManager::ConvolveFit(double* x, double* params)
{
	if (!x) std::cout << "x is nullptr\n";

	EnergyDistributionSet& set = energyDistributionSets.at(currentSetIndex);
	RateCoefficient& rc = rateCoefficientList.at(currentRateCoefficientIndex);

	double detuningEnergy = x[0];
	double sum = 0;
	
	// find correct distribution
	//std::cout << "Ed: " << detuningEnergy << "\n";
	int index = rc.GetIndexOfDetuningEnergy(detuningEnergy);
	if (index) return 0.0;

	EnergyDistribution& distribution = set.distributions.at(index);
	
	for (int i = 0; i < distribution.psi.size(); i++)
	{
		//std::cout << i << " ";
		//std::cout << distribution.GetPsis()[i] << " ";
		//std::cout << params[i] << "\n";
		sum += distribution.psi[i] * params[i];
	}
	//std::cout << "sum " << sum << "\n";
	return sum;
}

void DeconvolutionManager::AddRateCoefficientToList(RateCoefficient& rc)
{
	rateCoefficientList.emplace_back(std::move(rc));
}

void DeconvolutionManager::RemoveRateCoefficient(int index)
{
	rateCoefficientList.erase(rateCoefficientList.begin() + index);
}

void DeconvolutionManager::AddCrossSectionToList(CrossSection& cs)
{
	std::cout << "cs energy size: " << cs.energies.size() << std::endl;
	crossSectionList.emplace_back(std::move(cs));
}

void DeconvolutionManager::RemoveCrossSection(int index)
{
	crossSectionList.erase(crossSectionList.begin() + index);
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

