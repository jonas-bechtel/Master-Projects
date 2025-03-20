#pragma once

#include "ParameterImplementations.h"
#include "CoolingForceData.h"

using RNG_engine = std::mersenne_twister_engine<std::uint_fast64_t,
	64, 312, 156, 31,
	0xb5026f5aa96619e9, 29,
	0x5555555555555555, 17,
	0x71d67fffeda60000, 37,
	0xfff7eee000000000, 43,
	6364136223846793005>;

struct BinningSettings
{
	float energyRange[2] = { 1e-6f, 100.0f };

	bool constantBinSize = false;
	double normalStepSize = 0.5;
	double peakStepSize = 0.01;

	bool factorBinning = true;
	int binsPerDecade = 20;
	int binsAtPeak = 50;

	bool increasePeakResolution = true;

	void ShowWindow(bool& show);
};

struct PeakFitSettings
{
	double initialRange[2] = {0.0, 100.0};

	int fitRounds = 4;
	bool freeDetuningEnergy[4] = { true, false, true, true };
	bool freekT_long[4] = { false, true, true, true };
	bool freeKT_trans[4] = { false, false, false, false };
	bool adjustRange[4] = { true, true, true, true };
	bool showLimitedFit = true;

	void ShowWindow(bool& show);
};

class EnergyDistribution : public TH1D
{
public:
	EnergyDistribution();
	~EnergyDistribution();

	EnergyDistribution(const EnergyDistribution& other) = delete;
	EnergyDistribution& operator=(const EnergyDistribution& other) = delete;
	EnergyDistribution(EnergyDistribution&& other) noexcept;
	EnergyDistribution& operator=(EnergyDistribution&& other) noexcept;

	void Generate(std::filesystem::path descriptionFile, int index, const BinningSettings& binSettings, const PeakFitSettings& fitSettings);
	void CalculatePsisFromBinning(TH1D* crossSection);
	void Plot(bool showMarkers, bool showFit) const;

	void SetPlot(bool plot);
	void SetNormalised(bool normalised);

	double GetDetuningEnergy();

	void ShowListItem();

	void SaveSamples(std::filesystem::path folder) const;
	void SaveHist(std::filesystem::path folder) const;
	void Load(std::filesystem::path& file, bool loadSamples);
	
private:
	void CopyParameters();
	void ResetDefaultValues();
	void SetupLabellingThings();
	void SetupBinning(const BinningSettings& binSettings);
	void FillVectorsFromHist();
	void RemoveEdgeZeros();
	void CalculateFWHM();
	void FitAnalyticalToPeak(const PeakFitSettings& settings);

	std::string HeaderString() const;
	std::string Filename() const;

private:
	// actual distribution data
	std::vector<double> collisionEnergies;
	std::vector<double> binCenters;
	std::vector<double> binValues;
	std::vector<double> binValuesNormalised;
	std::vector<double> fitX;
	std::vector<double> fitY;

	// all the parameters used to create it
	MCMC_Parameters mcmcParameter;
	ElectronBeamParameters eBeamParameter;
	IonBeamParameters ionBeamParameter;
	LabEnergyParameters labEnergiesParameter;

	// output stuff
	OutputParameters outputParameter;

	// additional labelling things
	std::string label = "";
	std::string tags = "";
	int index = 0;

	// fitting things done by the Cross Section class
	std::vector<double> psi;

	// cooling force
	CoolingForceData cfData;

	// plot parameters
	bool showPlot = false;
	bool showNormalisedByWidth = true;

private:
	// random number generation things
	static std::mersenne_twister_engine<std::uint_fast64_t,
		64, 312, 156, 31,
		0xb5026f5aa96619e9, 29,
		0x5555555555555555, 17,
		0x71d67fffeda60000, 37,
		0xfff7eee000000000, 43,
		6364136223846793005> generator;

	static std::normal_distribution<double> longitudinalNormalDistribution;
	static std::normal_distribution<double> transverseNormalDistribution;

	friend struct SetInformation;
	friend class CrossSection;
	friend class RateCoefficient;
};
