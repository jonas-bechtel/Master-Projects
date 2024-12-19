#pragma once
#include "Module.h"
#include "CrossSection.h"
#include "RateCoefficient.h"

class DeconvolutionManager : public CrossSectionDeconvolutionModule
{
	//enum Binning { PaperBinning, ConstantBinning, FactorBinning, PaperFactorMix, Paper_FWHM };

public:
	DeconvolutionManager();

	void ShowRateCoefficientListWindow();
	void ShowCrossSectionListWindow();
	void ShowPlasmaRateListWindow();

private:
	void ShowUI() override;
	void ShowSettings();
	void ShowPlots();
	void SetupFitOptionsPopup();
	void SetupBinningOptionsPopup();

	CrossSection Deconvolve(const RateCoefficient& rc, EnergyDistributionSet& set);
	void DeconvolveInPlace(const RateCoefficient& rc, const EnergyDistributionSet& set, CrossSection& cs);
	double* DeconvolveWithSVD(const RateCoefficient& rc, const EnergyDistributionSet& set, const CrossSection& cs);
	double* DeconvolveWithROOT(const RateCoefficient& rc, const EnergyDistributionSet& set, const CrossSection& cs);
	double* DeconvolveWithGradientDescent(const RateCoefficient& rc, const EnergyDistributionSet& set, const CrossSection& cs);

	RateCoefficient Convolve(const CrossSection& cs, EnergyDistributionSet& set);
	void ConvolveInPlace(const CrossSection& cs, const EnergyDistributionSet& set, RateCoefficient& rc);

	double ConvolveFit(double* x, double* param);

	double MaxwellBoltzmannDistribution(double energy, double temperature);
	PlasmaRateCoefficient ConvolveIntoPlasmaRate(const CrossSection& cs);

	void AddRateCoefficientToList(RateCoefficient& rc);
	void RemoveRateCoefficient(int index);

	void AddCrossSectionToList(CrossSection& cs);
	void RemoveCrossSection(int index);

	void AddPlasmaRateToList(PlasmaRateCoefficient& prc);
	void RemovePlasmaRate(int index);

	void PlotRateCoefficient();
	void PlotCrossSections();

private:

	const char* binningOptions[4] = { "paper binning", "factor binning", "paper/factor mix", "paper/FWHM" };
	CrossSectionBinningSettings currentSettings;

	// fit options
	bool ROOT_fit = true;
	bool SVD_fit = false;
	bool limitROOTparameterRange = false;
	bool GD_fit = false;
	int iterations = 100000;
	double learningRate = 1;

	// plot parameters
	bool logX = true;
	bool logY = true;
	bool showMarkers = false;

	ImVec4 inputColor = ImVec4(0.6, 0.2, 0.1, 1.0);

	char CSnameInput[64] = "cross section name";
	char RCnameInput[64] = "rate coefficient name";
};

