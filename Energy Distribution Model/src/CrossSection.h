#pragma once

#include "Module.h"

class CrossSection : public Module
{
	enum Binning {PaperBinning, ConstantBinning, FactorBinning, PaperFactorMix, FWHM};

public:
	CrossSection();

private:
	void ShowUI() override;
	void SetupTrueCrossSection();

	void test();

	void SetupFitCrossSectionHist();
	void CalculateRateCoefficients();
	void CalculatePsis();
	void SetupInitialGuess();
	void FitCrossSectionHistogram();
	void FitWithSVD();
	void FitWithEigenSVD();
	void FitWithEigenGD();
	double FitFunction(double* x, double* params);

	void FillFitPlots(double* crossSectionParamater);

	void PlotRateCoefficients();

private:
	std::vector<double> binCentersTrue;
	std::vector<double> binValuesTrue;

	std::vector<double> initialGuess;
	TH1D* crossSectionFit = nullptr;
	std::vector<double> binCentersFit;
	std::vector<double> binValuesFit;

	TGraph* rateCoefficients = new TGraph();
	TGraph* rateCoefficientsFit = new TGraph();

	const char* binningOptions[4] = {"paper binning" , "const binning", "factor binning", "paper/factor mix"};
	//bool optionUsed[3] = { true, false ,false };
	int currentOption = PaperFactorMix;

	int numberBins = 50;
	bool limitBinSize = false;
	double minBinSize = 1e-5;
	int fixParamStart = 0;
	int fixParamStop = 0;

	// plot parameters
	bool logX = true;
	bool logY = true;

	// test 
	double learningRate = 1e-1;
	int iterations = 100000;
	double lambda = 1;
};

