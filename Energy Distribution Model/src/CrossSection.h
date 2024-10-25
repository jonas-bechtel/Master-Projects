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

	void PlotRateCoefficients();

private:
	TH1D* crossSection = nullptr;
	std::vector<double> binCenters;
	std::vector<double> binValues;
	int nBins = 100;

	TGraph* rateCoefficients1 = nullptr;
	TGraph* rateCoefficients2 = nullptr;

	bool useSigmaHist = false;
};

