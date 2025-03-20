#pragma once

class HeatMapData
{
public:
	void FromTH3D(TH3D* hist, float zSliceValue);
	void Plot(std::string title) const;

private:
	std::vector<double> values;
	int nRows;
	int nCols;
	double minValue;
	double maxValue;
	ImPlotPoint bottomLeft;
	ImPlotPoint topRight;
};
