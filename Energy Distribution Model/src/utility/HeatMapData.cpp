#include "pch.h"
#include "HeatMapData.h"

void HeatMapData::FromTH3D(TH3D* hist, float zSliceValue)
{
	if (!hist)
		return;

	nRows = hist->GetNbinsY();
	nCols = hist->GetNbinsX();
	//std::cout << nRows << ", " << nCols << std::endl;
	values.clear();
	values.reserve(nRows * nCols);

	int z_bin = hist->GetZaxis()->FindBin(zSliceValue);

	for (int y_bin = nRows; y_bin >= 1; y_bin--)
	{
		for (int x_bin = 1; x_bin <= nCols; x_bin++)
		{
			double content = hist->GetBinContent(x_bin, y_bin, z_bin);
			values.push_back(content);
		}
	}

	minValue = *std::min_element(values.begin(), values.end());
	maxValue = *std::max_element(values.begin(), values.end());

	bottomLeft = { hist->GetXaxis()->GetBinLowEdge(1), hist->GetYaxis()->GetBinLowEdge(1) };
	topRight = { hist->GetXaxis()->GetXmax(), hist->GetYaxis()->GetXmax() };
}

void HeatMapData::Plot(std::string label, std::string title) const
{
	ImPlot::PushColormap(9);
	
	if (ImPlot::BeginPlot(title.c_str(), ImVec2(350, 350)))
	{
		ImPlot::SetupAxes("x", "y");
		ImPlot::PlotHeatmap(label.c_str(), values.data(), nRows, nCols, minValue, maxValue, nullptr, bottomLeft, topRight);
		if(ImPlot::IsPlotHovered())
			ShowColorValueTooltip();
		ImPlot::EndPlot();
	}
	ImGui::SameLine();
	ImPlot::ColormapScale("##HeatScale", minValue, maxValue, ImVec2(0, 350));

	ImPlot::PopColormap();
}

void HeatMapData::ShowColorValueTooltip() const
{
	ImPlotPoint mouse = ImPlot::GetPlotMousePos();
	
	// Convert mouse X, Y to array indices
	int i = static_cast<int>((mouse.x - bottomLeft.x) / ((topRight.x - bottomLeft.x) / nCols));
	int j = static_cast<int>((mouse.y - bottomLeft.y) / ((topRight.y - bottomLeft.y) / nRows));

	// Flip Y-axis indexing (ImPlot uses bottom-left origin)
	j = (nRows - 1) - j;

	// Check if within valid range
	if (i >= 0 && i < nCols && j >= 0 && j < nRows) {
		float z_value = values[j * nCols + i]; // Correct indexing

		// Show tooltip with X, Y, and Z values
		ImGui::BeginTooltip();
		ImGui::Text("X: %.3f", mouse.x);
		ImGui::Text("Y: %.3f", mouse.y);
		ImGui::Text("Z: %.4e", z_value);
		ImGui::EndTooltip();
	}
}
