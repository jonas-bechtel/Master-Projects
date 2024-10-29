#pragma once
#include "Module.h"

#include <TGraph.h>

class CrossSection : public Module
{
public:
	CrossSection();

private:
	void ShowUI() override;
	void GenerateCrossSection();

	void CalculateRateCoefficients();
	void CalculatePsis();
	void FitCrossSectionHistogram();
	double FitFunction(double* x, double* params);

	void PlotRateCoefficients();

private:
	TH1D* crossSection = nullptr;
	std::vector<double> binCenters;
	std::vector<double> binValues;
	int nBins = 10;

	TH1D* crossSectionFit = nullptr;
	std::vector<double> binCentersFit;
	std::vector<double> binValuesFit;

	TGraph* rateCoefficients = new TGraph();
	TGraph* rateCoefficientsFit = new TGraph();

	bool useSigmaHist = false;

	// plot parameters
	bool logX = true;
	bool logY = true;
};

