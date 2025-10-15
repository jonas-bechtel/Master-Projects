#pragma once

class HeatMapData
{
public:
	void FromTH3D(TH3D* hist, float zSliceValue);
	void Plot(std::string label, std::string title) const;

private:
	void ShowColorValueTooltip() const;

private:
	std::vector<double> values;
	int nRows = 0;
	int nCols = 0;
	double minValue = 0;
	double maxValue = 1;
	ImPlotPoint bottomLeft = { 0,0 };
	ImPlotPoint topRight = { 1,1 };
};
