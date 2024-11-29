#pragma once

#include "Module.h"

class CrossSection : public Module
{
	enum Binning {PaperBinning, ConstantBinning, FactorBinning, PaperFactorMix, Paper_FWHM};

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
	void FitWithTorch();
	torch::Tensor custom_loss(const torch::Tensor& x, const torch::Tensor& A, const torch::Tensor& b);
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

	const char* binningOptions[5] = {"paper binning" , "const binning", "factor binning", "paper/factor mix", "paper/FWHM"};
	//bool optionUsed[3] = { true, false ,false };
	int currentOption = PaperFactorMix;

	double binFactor = 1;

	int numberBins = 50;
	bool limitBinSize = false;
	double minBinSize = 1e-5;

	// plot parameters
	bool logX = true;
	bool logY = true;

	//fitting things
	bool limitParamRange = false;

	// test 
	double learningRate = 1e-1;
	int iterations = 100000;
	double lambda = 1;

	// torch parameter
	double torchLearningRate = 1e-14;
	int nEpochs = 10000;
	double l2regularisation = 0.01;
	double smoothRegularisation = 0.01;
};

