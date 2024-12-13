#pragma once
#include "Module.h"
#include "CrossSection.h"
#include "RateCoefficient.h"

class DeconvolutionManager : public CrossSectionDeconvolutionModule
{
	enum Binning { PaperBinning, ConstantBinning, FactorBinning, PaperFactorMix, Paper_FWHM };

public:
	DeconvolutionManager();

	void ShowRateCoefficientListWindow();
	void ShowCrossSectionListWindow();

private:
	void ShowUI() override;

	CrossSection Deconvolve(const RateCoefficient& rc, EnergyDistributionSet& set);
	void DeconvolveInPlace(const RateCoefficient& rc, const EnergyDistributionSet& set, CrossSection& cs);

	RateCoefficient Convolve(CrossSection& cs, EnergyDistributionSet& set);
	void ConvolveInPlace(const CrossSection& cs, const EnergyDistributionSet& set, RateCoefficient& rc);

	double ConvolveFit(double* x, double* param);

	void AddRateCoefficientToList(RateCoefficient& rc);
	void RemoveRateCoefficient(int index);

	void AddCrossSectionToList(CrossSection& cs);
	void RemoveCrossSection(int index);

	void PlotRateCoefficient();

private:

	const char* binningOptions[4] = { "paper binning", "factor binning", "paper/factor mix", "paper/FWHM" };
	CrossSectionBinningSettings currentSettings;

	// plot parameters
	bool logX = true;
	bool logY = true;
};

