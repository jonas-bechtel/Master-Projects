#pragma once

#include "EnergyDistributionManager.h"

#include <TH1D.h>

class EnergyDistribution : public TH1D
{
public:
	EnergyDistribution();
	std::vector<double>& GetCollisionEnergies();
	std::vector<double>& GetBinCenters();
	std::vector<double>& GetBinValues();
	std::vector<double>& GetNormalisedBinValues();

	MCMC_Parameters& GetMCMC_Parameter();
	ElectronBeamParameters& GetElectronBeamParameter();
	IonBeamParameters& GetIonBeamParameter();
	LabEnergiesParameters& GetLabEnergyParameter();
	EnergyDistributionParameters& GetEnergyDistributionParameter();

	void SetLabel(std::string label);
	double& GetRateCoefficient();
	std::vector<double>& GetPsis();
	std::string& GetLabel();
	std::string& GetTags();
	std::filesystem::path& GetFolder();
	std::filesystem::path& GetSubFolder();

	bool& IsPlotted();
	bool& IsPlottedNormalised();
	
	void SetupFromCurrentEnvironment();
	void RemoveEdgeZeros();

	std::string String();
	std::string Filename();

	static EnergyDistribution* FindByEd(double detuningEnergy);
	static std::unordered_map<double, EnergyDistribution*> s_allDistributions;

private:

private:
	// actual distribution
	std::vector<double> collisionEnergies;
	std::vector<double> binCenters;
	std::vector<double> binValues;
	std::vector<double> binValuesNormalised;

	// all the things used to create it
	MCMC_Parameters mcmcParameter;
	ElectronBeamParameters eBeamParameter;
	IonBeamParameters ionBeamParameter;
	LabEnergiesParameters labEnergiesParameter;
	EnergyDistributionParameters eDistParameter;

	// additional labelling things
	std::string m_label = "";
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

