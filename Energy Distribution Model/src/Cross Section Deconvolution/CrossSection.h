#pragma once

struct RateCoefficient;

enum CrossSectionBinningScheme { PaperBinning, FactorBinning, PaperFactorMix, Paper_FWHM };


struct CrossSectionBinningSettings
{
	int numberBins = 100;

	CrossSectionBinningScheme scheme = CrossSectionBinningScheme::PaperFactorMix;
};


struct CrossSection
{
	CrossSection();
	CrossSection(const CrossSection& other) = delete;
	CrossSection& operator=(const CrossSection& other) = delete;
	CrossSection(CrossSection&& other) = default;
	CrossSection& operator=(CrossSection&& other) = default;

	void SetValues(double* newValues);
	void SetupBinning(CrossSectionBinningSettings binSettings);
	void SetupInitialGuess(const RateCoefficient& rc);

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

