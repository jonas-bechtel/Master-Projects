#pragma once

struct RateCoefficient;

enum CrossSectionBinningScheme { PaperBinning, FactorBinning, PaperFactorMix, Paper_FWHM };


struct CrossSectionBinningSettings
{
	int numberBins = 100;
	int maxRatio = 10;

	CrossSectionBinningScheme scheme = CrossSectionBinningScheme::PaperBinning;
};


struct CrossSection
{
	CrossSection();
	CrossSection(const CrossSection& other) = delete;
	CrossSection& operator=(const CrossSection& other) = delete;
	CrossSection(CrossSection&& other) = default;
	CrossSection& operator=(CrossSection&& other) = default;

	void SetValues(double* newValues, bool square = false);
	void SetupBinning(CrossSectionBinningSettings binSettings, const RateCoefficient& rc);
	void SetupInitialGuess(const RateCoefficient& rc, bool squareRoot = false);
	void FillWithOneOverE(int scale = 1);

public:
	// main data
	TH1D* hist;
	std::vector<double> energies;
	std::vector<double> values;
	std::vector<double> errors;

	std::vector<double> initialGuess;

	// labelling things
	std::string label = "cs";
	std::filesystem::path file;
	std::filesystem::path energyDistriubtionSetFolder;
	std::filesystem::path mergedBeamRateCoefficientFile;
};

