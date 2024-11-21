#pragma once

#include "MCMC.h"
#include "ElectronBeam.h"
#include "IonBeam.h"
#include "LabEnergies.h"

struct BinningSettings;

struct EnergyDistributionParameters : public Parameters
{
	EnergyDistributionParameters()
	{
		setName("energy distribution parameters");
	}

	ParameterValue<bool> cutOutZValues = ParameterValue(false, "cut out z values", "%d", true);
	ParameterValue<float2> cutOutRange = ParameterValue(float2(0.0f, 0.4f), "cut out range", "%.2f, %.2f m", true);

private:
	int GetSize() override
	{
		return sizeof(*this);
	}
};

struct AnalyticalDistributionParameters : public Parameters
{
	AnalyticalDistributionParameters()
	{
		setName("analytical/fit distribution parameters");
	}

	ParameterValue<double> detuningEnergy = ParameterValue(1.0, "detuning energy", "%.3f eV");
	ParameterValue<double> longitudinalTemperature = ParameterValue(0.0005, "longitudinal kT", "%.2e eV");
	ParameterValue<double> transverseTemperature = ParameterValue(0.002, "transverse kT", "%.2e eV");
	ParameterValue<double> FWHM = ParameterValue(0.0, "FWHM", "%.3f eV");
	ParameterValue<double> effectiveLength = ParameterValue(0.0, "effective length", "%.3f m");

private:
	int GetSize() override
	{
		return sizeof(*this);
	}
};

struct EnergyDistribution : public TH1D
{
	EnergyDistribution();
	~EnergyDistribution();
	
	void SetupLabellingThings();
	void SetupBinning(const BinningSettings& binSettings);
	void RemoveEdgeZeros();

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
	std::array<double, 200> fitX;
	std::array<double, 200> fitY;

	// all the things used to create it
	MCMC_Parameters mcmcParameter;
	ElectronBeamParameters eBeamParameter;
	IonBeamParameters ionBeamParameter;
	LabEnergyParameters labEnergiesParameter;
	EnergyDistributionParameters eDistParameter;
	AnalyticalDistributionParameters analyticalParameter;

	// additional labelling things
	std::string label = "";
	std::string tags = "";
	std::filesystem::path folder = "Test";
	std::filesystem::path subFolder = "";
	int index = 0;
	bool isAnalytical = false;

	// fitting things done by the Cross Section class
	double rateCoefficient = 0;
	std::vector<double> psi;

	// plot parameters
	bool plotted = false;
	bool showNormalisedByWidth = true;
};

