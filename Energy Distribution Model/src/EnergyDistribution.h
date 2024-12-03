#pragma once

#include "MCMC.h"
#include "ElectronBeam.h"
#include "IonBeam.h"
#include "LabEnergies.h"

struct BinningSettings;

struct EnergyDistribution : public TH1D
{
	EnergyDistribution();
	~EnergyDistribution();

	EnergyDistribution(const EnergyDistribution& other) = delete;
	EnergyDistribution& operator=(const EnergyDistribution& other) = delete;
	EnergyDistribution(EnergyDistribution&& other);
	EnergyDistribution& operator=(EnergyDistribution&& other);
	
	void ResetDefaultValues();
	void SetupLabellingThings();
	void SetupBinning(const BinningSettings& binSettings);
	void FillVectorsFromHist();
	void RemoveEdgeZeros();
	void FitAnalyticalToPeak();

	std::string String();
	std::string Filename();

	static EnergyDistribution* FindByEd(double detuningEnergy);
	static std::unordered_map<double, EnergyDistribution*> s_allDistributions;

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
	AnalyticalDistributionParameters analyticalParameter;
	SimplificationParameter simplifyParams;

	// additional labelling things
	std::string label = "";
	std::string tags = "";
	std::filesystem::path folder = "Test";
	std::filesystem::path subFolder = "";
	int index = 0;

	// fitting things done by the Cross Section class
	double rateCoefficient = 0;
	std::vector<double> psi;

	// plot parameters
	bool plotted = false;
	bool showNormalisedByWidth = true;
};

