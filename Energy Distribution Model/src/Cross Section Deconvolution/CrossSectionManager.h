#pragma once
#include "Module.h"
#include "CrossSection.h"

class CrossSectionManager : public CrossSectionDeconvolutionModule
{
	enum Binning { PaperBinning, ConstantBinning, FactorBinning, PaperFactorMix, Paper_FWHM };

public:
	CrossSectionManager();
	int GetCurrentBinningOption();
	std::string&& GetCurrentBinningString();

private:
	//void ShowUI() override;

	RateCoefficient ConvolveCrossSection(CrossSection& cs, EnergyDistributionSet& set);
	void ConvolveCrossSection(CrossSection& cs, EnergyDistributionSet& set, RateCoefficient& rc);


	void test();

	void SetupFitCrossSectionHist();
	void CalculateRateCoefficients();

	void SetupInitialGuess();
	void FitCrossSectionHistogram();
	void FitWithSVD();
	void FitWithEigenSVD();
	void FitWithEigenGD();
	void FitWithTorch();
	torch::Tensor custom_loss(const torch::Tensor& x, const torch::Tensor& A, const torch::Tensor& b);
	double FitFunction(double* x, double* params);

	void FillFitPlots(double* crossSectionParamater);

private:
	std::vector<double> initialGuess;

	const char* binningOptions[5] = { "paper binning" , "const binning", "factor binning", "paper/factor mix", "paper/FWHM" };
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

