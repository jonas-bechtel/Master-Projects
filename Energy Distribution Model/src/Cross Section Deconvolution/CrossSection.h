#pragma once

class EnergyDistributionSet;
class RateCoefficient;

enum CrossSectionBinningScheme { PaperBinning, FactorBinning, PaperFactorMix, Paper_FWHM };

struct CrossSectionBinningSettings
{
	const char* binningOptions[4] = { "paper binning", "factor binning", "paper/factor mix", "paper/FWHM" };

	int numberBins = 100;
	int maxRatio = 10;

	CrossSectionBinningScheme scheme = CrossSectionBinningScheme::PaperBinning;

	void ShowWindow(bool& show);
};

struct FittingOptions
{
	// fit options
	bool ROOT_fit = true;
	bool SVD_fit = false;
	bool GD_fit = false;
	int iterations = 100000;
	double learningRate = 1;

	void ShowWindow(bool& show);
};


class CrossSection
{
public:
	CrossSection();
	CrossSection(const CrossSection& other) = delete;
	CrossSection& operator=(const CrossSection& other) = delete;
	CrossSection(CrossSection&& other) = default;
	CrossSection& operator=(CrossSection&& other) = default;

	TH1D* GetHist();
	std::string GetLabel() const;

	void SetLabel(std::string label);
	void FillWithOneOverE(int scale = 1);
	void Deconvolve(const RateCoefficient& rc, EnergyDistributionSet& set, const FittingOptions& fitSettings, const CrossSectionBinningSettings& binSettings);

	void Plot(bool showMarkers) const;

	void Save() const;
	void Load(std::filesystem::path& filename);

private:
	void SetupBinning(const CrossSectionBinningSettings& binSettings, const RateCoefficient& rc);
	void SetupInitialGuess(const RateCoefficient& rc);

	void FitWithSVD(const RateCoefficient& rc, const EnergyDistributionSet& set);
	void FitWithROOT(const RateCoefficient& rc, const EnergyDistributionSet& set);

private:
	// main data
	TH1D* hist;
	std::vector<double> energies;
	std::vector<double> values;
	std::vector<double> errors;

	//std::vector<double> initialGuess;

	// labelling things
	std::string label = "cs";
	std::filesystem::path energyDistriubtionSetFolder;
	std::filesystem::path mergedBeamRateCoefficientFile;

	friend class RateCoefficient;
	friend class PlasmaRateCoefficient ;
};

