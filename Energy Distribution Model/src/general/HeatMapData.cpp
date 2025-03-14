#include "pch.h"
#include "HeatMapData.h"

void HeatMapData::Plot(std::string title) const
{
	ImPlot::PushColormap(9);
	
	if (ImPlot::BeginPlot("XY Slice", ImVec2(400, 400)))
	{
		ImPlot::SetupAxes("x", "y");
		ImPlot::PlotHeatmap(title.c_str(), values.data(), nRows, nCols, minValue, maxValue, nullptr, bottomLeft, topRight);
		ImPlot::EndPlot();
	}
	ImGui::SameLine();
	ImPlot::ColormapScale("##HeatScale", minValue, maxValue, ImVec2(0, 400));

	ImPlot::PopColormap();
}
