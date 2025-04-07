#pragma once

#include "HeatMapData.h"
#include "ROOTCanvas.h"

class PlotBeamData
{
private:
	TH3D* fullHistogram = nullptr;
	TH3D* fullHistogramSmall = nullptr;

	std::vector<double> xAxis;
	std::vector<double> yAxis;
	std::vector<double> zAxis;

	std::vector<double> projectionValuesX;
	std::vector<double> projectionValuesY;
	std::vector<double> projectionValuesZ;

	std::vector<double> centerValue;
	std::vector<double> outsideValue;

	HeatMapData slice;

	std::string label;

	static int rebinningFactors[3];

public:
	PlotBeamData() {};
	PlotBeamData(TH3D* hist);
	PlotBeamData(const PlotBeamData& other) = delete;
	PlotBeamData& operator=(const PlotBeamData& other) = delete;
	PlotBeamData(PlotBeamData&& other) noexcept;
	PlotBeamData& operator=(PlotBeamData&& other) noexcept;
	~PlotBeamData();

	TH3D* GetHist();
	std::string GetLabel();
	void SetLabel(std::string str);
	void UpdateSlice(float zValue);
	void UpdateData();

	void Plot3D(ROOTCanvas* canvas, int pos);
	void PlotSlice() const;
	void PlotProjectionX(ImPlotLineFlags_ flags = ImPlotLineFlags_None) const;
	void PlotProjectionY(ImPlotLineFlags_ flags = ImPlotLineFlags_None) const;
	void PlotProjectionZ(ImPlotLineFlags_ flags = ImPlotLineFlags_None) const;
	void PlotInsideOutsideValue() const;

	static bool ShowRebinningFactorsInput();
};

