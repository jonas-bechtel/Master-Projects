#pragma once

class EnergyDistributionSet;
class EnergyDistribution;
class RateCoefficient;

enum CrossSectionBinningScheme { PaperBinning, FactorBinning, PaperFactorMix };

struct CrossSectionBinningSettings
{
	const char* binningOptions[3] = { "paper binning", "factor binning", "paper/factor mix" };

	int numberBins = 7;
	int maxRatio = 10;
	double boundaryEnergy = 0.1;
	double binFactor = 1.5;

	CrossSectionBinningScheme scheme = CrossSectionBinningScheme::PaperFactorMix;

	void ShowWindow(bool& show);
};

struct FittingOptions
{
	// fit options
	bool ROOT_fit = true;
	bool SVD_fit = false;
	bool EigenNNLS_fit = false;
	bool NNLS_ROOT_fit = false;

	int errorIterations = 1; // number of iterations to calculate errors
	int fit_iterations = 1;
	int maxIterations = 1000;
	double tolerance = 1e-6;
	double learningRate = 1;
	bool fixParameters = false;
	float fixedParameterRange[2] = { 1,100 };

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
	void FillWithOneOverE(double scale = 1.0);
	void SetupBinning(const CrossSectionBinningSettings& binSettings, const RateCoefficient& rc);
	void SetInitialGuessValues(const RateCoefficient& rc);
	void Deconvolve(RateCoefficient& rc, EnergyDistributionSet& set, const FittingOptions& fitSettings, const CrossSectionBinningSettings& binSettings);

	void Plot(bool showMarkers) const;

	void Clear();
	void Save() const;
	void Load(std::filesystem::path& filename);

private:
	double ConvolveFit(double Ed, double* csBins, const EnergyDistributionSet& set,
		bool squareCS = true,
		std::unordered_map<double, EnergyDistribution*>& map = std::unordered_map<double, EnergyDistribution*>()) const;
	void FitWithSVD(const RateCoefficient& rc, const EnergyDistributionSet& set);
	void FitWithROOT(const RateCoefficient& rc, const EnergyDistributionSet& set, const FittingOptions& fitSettings);
	void FitWithEigenNNLS(const RateCoefficient& rc, const EnergyDistributionSet& set, const FittingOptions& fitSettings);

	void ResetNonFixedParameters(const RateCoefficient& rc, const FittingOptions& fitSettings);
private:
	// main data
	TH1D* hist = nullptr;

	std::vector<double> energies;
	std::vector<double> values;
	std::vector<double> errors;

	// array with values from error iterations
	std::vector<double> valueArray;

	// labelling things
	std::string label = "cs";
	std::filesystem::path energyDistriubtionSetFolder;
	std::filesystem::path mergedBeamRateCoefficientFile;

	friend class RateCoefficient;
	friend class PlasmaRateCoefficient ;
};

