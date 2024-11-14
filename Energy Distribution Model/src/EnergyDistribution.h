#pragma once

#include <TH1D.h>

#include "MCMC.h"
#include "ElectronBeam.h"
#include "IonBeam.h"
#include "LabEnergies.h"

struct EnergyDistributionParameters : public Parameters
{
	EnergyDistributionParameters()
	{
		setName("energy distribution parameters");
	}

	ParameterValue<bool> limitBinSize = ParameterValue(false, "limit bin size", "%d");
	ParameterValue<double> minBinSize = ParameterValue(5.0, "min bin size", "%.1e eV");
	ParameterValue<bool> cutOutZValues = ParameterValue(false, "cut out z values", "%d");
	ParameterValue<float2> cutOutRange = ParameterValue(float2(0.0f, 0.4f), "cut out range", "%.2f, %.2f m");

private:
	int GetSize() override
	{
		return sizeof(*this);
	}
};

struct EnergyDistribution : public TH1D
{
	EnergyDistribution();
	
	void SetupFromCurrentEnvironment();
	void SetupFromHeader(std::string& header);
	void RemoveEdgeZeros();

	std::string String();
	std::string Filename();

	static EnergyDistribution* FindByEd(double detuningEnergy);
	static std::unordered_map<double, EnergyDistribution*> s_allDistributions;

	// actual distribution data
	std::vector<double> collisionEnergies;
	std::vector<double> binCenters;
	std::vector<double> binValues;
	std::vector<double> binValuesNormalised;

	// all the things used to create it
	MCMC_Parameters mcmcParameter;
	ElectronBeamParameters eBeamParameter;
	IonBeamParameters ionBeamParameter;
	LabEnergyParameters labEnergiesParameter;
	EnergyDistributionParameters eDistParameter;

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
	bool showNormalisedByWidth = false;
};

