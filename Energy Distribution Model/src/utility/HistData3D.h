#pragma once

#include "HeatMapData.h"
#include "ROOTCanvas.h"

class HistData3D
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
	HistData3D() {};
	HistData3D(TH3D* hist);
	HistData3D(const HistData3D& other) = delete;
	HistData3D& operator=(const HistData3D& other) = delete;
	HistData3D(HistData3D&& other) noexcept;
	HistData3D& operator=(HistData3D&& other) noexcept;
	~HistData3D();

	void Setup(TH3D* hist);

	TH3D* GetHist() const;
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

