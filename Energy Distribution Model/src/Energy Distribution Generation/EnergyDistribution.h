#pragma once

#include "ParameterImplementations.h"
#include "CoolingForceData.h"

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
};

struct EnergyDistribution : public TH1D
{
	EnergyDistribution();
	~EnergyDistribution();

	EnergyDistribution(const EnergyDistribution& other) = delete;// { std::cout << "using illegal copy constructor" << std::endl; }// = delete;
	EnergyDistribution& operator=(const EnergyDistribution& other) = delete;
	EnergyDistribution(EnergyDistribution&& other) noexcept;
	EnergyDistribution& operator=(EnergyDistribution&& other) noexcept;
	
	void ResetDefaultValues();
	void SetupLabellingThings();
	void SetupBinning(const BinningSettings& binSettings);
	void FillVectorsFromHist();
	void RemoveEdgeZeros();
	void CalculateFWHM();
	void FitAnalyticalToPeak(const PeakFitSettings& settings);
	void CalculatePsisFromBinning(TH1D* crossSection);

	std::string String() const;
	std::string Filename() const;

public:
	// actual distribution data
	std::vector<double> collisionEnergies;
	std::vector<double> binCenters;
	std::vector<double> binValues;
	std::vector<double> binValuesNormalised;
	std::vector<double> fitX;
	std::vector<double> fitY;

	// all the things used to create it
	MCMC_Parameters mcmcParameter;
	ElectronBeamParameters eBeamParameter;
	IonBeamParameters ionBeamParameter;
	LabEnergyParameters labEnergiesParameter;
	OutputParameters outputParameter;
	SimplificationParameter simplifyParams;

	// additional labelling things
	std::string label = "";
	std::string tags = "";
	int index = 0;

	// fitting things done by the Cross Section class
	std::vector<double> psi;

	// cooling force
	CoolingForceData cfData;

	// plot parameters
	bool plotted = false;
	bool showNormalisedByWidth = true;
};

struct SetInformation
{
	std::vector<int> indeces;
	std::vector<double> centerLabEnergy;
	std::vector<double> detuningEnergy;
	std::vector<double> fitDetuningEnergy;
	std::vector<double> fitLongitudinalTemperature;
	std::vector<double> fitTransverseTemperature;
	std::vector<double> fitScalingFactor;
	std::vector<double> fitFWHM;
	std::vector<double> FWHM;
	std::vector<double> detuningVelocity;
	std::vector<double> longitudinalCoolingForce;
	//std::vector<double> effectiveLength;

	void AddDistributionValues(const EnergyDistribution& dist);
};

struct EnergyDistributionSet
{
	EnergyDistributionSet()
	{
		distributions.reserve(100);
		//std::cout << "calling EnergyDistSet default constructor" << std::endl;
	}
	EnergyDistributionSet(const EnergyDistributionSet& other) = delete;
	EnergyDistributionSet& operator=(const EnergyDistributionSet& other) = delete;

	EnergyDistributionSet(EnergyDistributionSet&& other) noexcept
	{
		distributions = std::move(other.distributions);
		EdToDistMap = std::move(other.EdToDistMap);
		info = std::move(other.info);
		folder = std::move(other.folder);
		subFolder = std::move(other.subFolder);

		//std::cout << "calling EnergyDistSet move constructor" << std::endl;
	}
	EnergyDistributionSet& operator=(EnergyDistributionSet&& other) noexcept
	{
		distributions = std::move(other.distributions);
		EdToDistMap = std::move(other.EdToDistMap);
		info = std::move(other.info);
		folder = std::move(other.folder);
		subFolder = std::move(other.subFolder);

		//std::cout << "calling EnergyDistSet move assignment op" << std::endl;
		return *this;
	}

	void AddDistribution(EnergyDistribution&& distribution);
	std::string Label()
	{
		return (folder / subFolder).string();
	}
	EnergyDistribution* FindByEd(double detuningEnergy);

	void SetAllPlotted(bool plotted);
	void SetAllShowNormalised(bool showNormalised);

	void CalculatePsisFromBinning(TH1D* crossSection);
	
public:
	std::vector<EnergyDistribution> distributions;
	std::unordered_map<double, EnergyDistribution*> EdToDistMap;
	SetInformation info;

	std::filesystem::path folder = "Test";
	std::filesystem::path subFolder = "subfolder";

	bool plotInfo = false;
};

