#include "pch.h"
#include "RateCoefficient.h"
#include "FileHandler.h"

RateCoefficient::RateCoefficient()
	: Module("rate coefficient")
{
}

void RateCoefficient::ShowUI()
{
	if (ImGui::Button("load rate coefficients"))
	{
		std::filesystem::path file = FileHandler::GetInstance().SelectFile("data\\RateCoefficients\\", {"*.dat"});
		std::vector<mbrcData> data = FileHandler::GetInstance().LoadRateCoefficients(file);
		FillVectors(data);

		PlotRateCoefficient();
	}

	if (ImPlot::BeginPlot("rate coefficient"))
	{
		ImPlot::SetupAxis(ImAxis_X1, "detuning energy [eV]");
		ImPlot::SetupAxis(ImAxis_Y1, "rate coefficient");
		if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
		if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

		ImPlot::SetupLegend(ImPlotLocation_NorthEast);

		ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
		ImPlot::PlotErrorBars("##errors", detuningEnergies.data(), rateCoefficients.data(), rateCoefficientsError.data(), rateCoefficients.size());
		ImPlot::PlotLine("rate coefficient", detuningEnergies.data(), rateCoefficients.data(), rateCoefficients.size());
		
		ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
		ImPlot::PlotLine("rate coefficient fit", detuningEnergies.data(), rateCoefficientsFit.data(), rateCoefficientsFit.size());

		ImPlot::EndPlot();
	}

	ImGui::Checkbox("log X", &logX);
	ImGui::SameLine();
	ImGui::Checkbox("log Y", &logY);

}

void RateCoefficient::FillVectors(std::vector<mbrcData>& data)
{
	detuningEnergies.clear();
	rateCoefficients.clear();
	rateCoefficientsError.clear();

	detuningEnergies.reserve(data.size());
	rateCoefficients.reserve(data.size());
	rateCoefficientsError.reserve(data.size());

	for (const mbrcData& thing : data)
	{
		detuningEnergies.push_back(thing.E_det);
		rateCoefficients.push_back(thing.rateCoef);
		rateCoefficientsError.push_back(thing.RateCoefError);
		rateCoefficientsGraph->AddPoint(thing.E_det, thing.rateCoef);
		
	}
}

void RateCoefficient::PlotRateCoefficient()
{
	m_mainCanvas->cd(1);

	rateCoefficientsGraph->SetLineColor(kBlue);
	rateCoefficientsGraph->SetMarkerStyle(21);
	rateCoefficientsGraph->SetTitle("rate Coefficients");
	rateCoefficientsGraph->GetXaxis()->SetTitle("E_d [eV]");
	rateCoefficientsGraph->GetYaxis()->SetTitle("alpha [m^3/s]");
	rateCoefficientsGraph->Draw("ALP");
}
