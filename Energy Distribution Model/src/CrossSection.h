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

	void test();

	void SetupFitCrossSectionHist();
	void CalculateRateCoefficients();
	void CalculatePsis();
	void FitCrossSectionHistogram();
	double FitFunction(double* x, double* params);

	void FillFitPlots(double* crossSectionParamater);

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

	const char* binningOptions[3] = {"paper binning" , "const binning", "factor binning"};
	//bool optionUsed[3] = { true, false ,false };
	int currentOption = 0;

	int numberBins = 100;
	bool limitBinSize = false;
	double minBinSize = 1e-5;

	// plot parameters
	bool logX = true;
	bool logY = true;
};

