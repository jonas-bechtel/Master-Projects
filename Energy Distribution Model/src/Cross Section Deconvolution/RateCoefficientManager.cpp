#include "pch.h"
#include "RateCoefficientManager.h"
#include "RateCoefficient.h"
#include "FileHandler.h"

RateCoefficientManager::RateCoefficientManager()
	: CrossSectionDeconvolutionModule("rate coefficient")
{
}

void RateCoefficientManager::ShowUI()
{
	if (ImGui::Button("load rate coefficients"))
	{
		std::filesystem::path file = FileHandler::GetInstance().SelectFile("data\\RateCoefficients\\", { "*.dat" });
		RateCoefficient rc = FileHandler::GetInstance().LoadRateCoefficients(file);
		AddRateCoefficientToList(rc);

		PlotRateCoefficient();
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
			ImPlot::PlotErrorBars("##errors", rc.detuningEnergies.data(), rc.value.data(), rc.error.data(), rc.value.size());
			ImPlot::PlotLine(rc.label.c_str(), rc.detuningEnergies.data(), rc.value.data(), rc.value.size());
		}
		
		ImPlot::EndPlot();
	}

	ImGui::Checkbox("log X", &logX);
	ImGui::SameLine();
	ImGui::Checkbox("log Y", &logY);

}

void RateCoefficientManager::AddRateCoefficientToList(RateCoefficient& rc)
{
	rateCoefficientList.emplace_back(std::move(rc));
}

void RateCoefficientManager::RemoveRateCoefficient(int index)
{
	rateCoefficientList.erase(rateCoefficientList.begin() + index);
}

void RateCoefficientManager::PlotRateCoefficient()
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

