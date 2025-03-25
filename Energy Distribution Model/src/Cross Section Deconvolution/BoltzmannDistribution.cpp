#include "pch.h"
#include "BoltzmannDistribution.h"
#include "Constants.h"

namespace BoltzmannDistribution
{
	static float temperature = 100.0f;
	static float energyRange[2] = { 1e-6f, 100.0f };
	static float energies[2000];
	static float values[2000];
	static float valuesMultiplied[2000];
	static float color[3] = { 1.0f, 0.0f, 0.0f };

	double Function(double energy, double temperature)
	{
		double kB_T = ((TMath::K() / TMath::Qe()) * temperature);

		return 2.0 * sqrt(energy / TMath::Pi()) * pow(1.0 / kB_T, 1.5) * exp(-energy / kB_T);
	}

	void Update(CrossSection* currentCS)
	{
		float step = (energyRange[1] - energyRange[0]) / 1999;
		for (int i = 0; i < 2000; i++)
		{
			energies[i] = energyRange[0] + i * step;
			values[i] = Function(energies[i], temperature);
			float velocity = TMath::Sqrt(2 * energies[i] * TMath::Qe() / PhysicalConstants::electronMass);
			float sigma = 0;
			if (currentCS)
			{
				sigma = currentCS->GetHist()->Interpolate(energies[i]);
			}
			valuesMultiplied[i] = values[i] * velocity * sigma;
		}
	}

	void Plot(bool showMarkers)
	{
		if (showMarkers) ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);

		ImPlot::PushStyleColor(ImPlotCol_Line, { color[0], color[1], color[2], 1.0 });
		ImPlot::PushStyleVar(ImPlotStyleVar_FillAlpha, 0.25f);
		ImPlot::PlotLine("Maxwell Boltzmann distribution", energies, values, 2000);
		ImPlot::PlotLine("sigma * v * f_pl", energies, valuesMultiplied, 2000, ImPlotLineFlags_Shaded);
		ImPlot::PopStyleVar();
		ImPlot::PopStyleColor();
	}

	void ShowWindow(bool& show, CrossSection* currentCS)
	{
		if (!show)
		{
			return;
		}
		if (ImGui::Begin("plasma rate convolution extra window", &show, ImGuiWindowFlags_NoDocking))
		{
			bool changed = false;

			changed |= ImGui::SliderFloat("Temperature", &temperature, 1.0f, 5000.0f, "%.3f", ImGuiSliderFlags_Logarithmic);
			changed |= ImGui::SliderFloat2("range", energyRange, 1e-6f, 100.0f, "%.6f", ImGuiSliderFlags_Logarithmic);

			ImGui::ColorEdit3("color", color);

			if (changed)
			{
				Update(currentCS);
			}
		}
		ImGui::End();
	}
}