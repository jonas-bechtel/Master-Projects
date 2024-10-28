#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <filesystem>

#include <TH3D.h>

#include "MCMC.h"
#include "ElectronBeam.h"
#include "IonBeam.h"
#include "LabEnergies.h"
#include "Point3D.h"

struct EnergyDistributionParameters
{
	double driftTubeVoltage = 0;
	bool normaliseByBinWidth = true;
	bool energyDefinedBinning = true;

	// parameters for simpler test 
	bool cutOutZValues = false;
	float cutOutRange[2] = { 0, 0.35 };

	std::string String();
};

struct EnergyDistribution
{
	// resulting values
	double rateCoefficient = 0;

	// actual distribution
	TH1D* distribution = nullptr;
	std::vector<double> collisionEnergies;
	std::vector<double> binCenters;
	std::vector<double> binValues;

	// all the things used to create it
	MCMC_Parameters mcmcParameter;
	ElectronBeamParameters eBeamParameter;
	IonBeamParameters ionBeamParameter;
	LabEnergiesParameters labEnergiesParameter;
	EnergyDistributionParameters eDistParameter;

	// additional labelling things
	std::string label = "";
	std::string tags = "";
	//std::filesystem::path descriptionFile;
	std::filesystem::path folder = "Test";
	int index = 0;

	//double cathodeVoltage;

	bool plotted = false;

	TH1D* operator->()
	{
		return distribution;
	}

	std::string String();
	std::string Filename();

	//EnergyDistribution()
	//{
	//	std::cout << "creating energy distribution\n";
	//}
	//~EnergyDistribution()
	//{
	//	std::cout << "deleting energy distribution\n";
	//	delete distribution;
	//}
};

class EnergyDistributionModel : public Module
{
public:
	EnergyDistributionModel();
	std::vector<EnergyDistribution>& GetEnergyDistributions();
	void GenerateEnergyDistribution();
	void GenerateEnergyDistributionsFromFile(std::filesystem::path file);

private:
	void ShowUI() override;
	void ShowEnergyDistributionList();

	void SetupEnergyDistribution();

	void PlotEnergyDistributions();
	void PLotZweightByEnergy();
	void PlotLongkTDistribution();
	void PlotLongVelAddition();
	void PlotRateCoefficients();

	void ClearDistributionList();

private:
	std::vector<EnergyDistribution> energyDistributions;
	EnergyDistribution currentDistribution;
	EnergyDistributionParameters parameter;
	
	// graphs and plots
	TGraph* rateCoefficients = nullptr;

	TH1D* zPositions = nullptr;
	TH1D* zWeightByEnergy = nullptr;
	TH1D* long_ktDistribution = nullptr;
	TH1D* long_VelAddition = nullptr;

	// currently loaded files
	std::filesystem::path currentDescriptionFile = std::filesystem::path("data\\C60\\100x100x100_Ie0.012_Ucath47.0_RelTol-1e-3_sort_energies.asc");
	int maxIndex = 0;

	// start/end index in description file to generate distribution for
	int startIndex = 1;
	int endIndex = 1;
	bool doAll = false;

	// parameters for energy distribution generation
	float energyRange[2] = { 1e-8, 100 };
	int binsPerDecade = 20000;

	// plot parameters
	bool logX = true;
	bool logY = true;

	// random number generation things
	std::mersenne_twister_engine<std::uint_fast64_t,
		64, 312, 156, 31,
		0xb5026f5aa96619e9, 29,
		0x5555555555555555, 17,
		0x71d67fffeda60000, 37,
		0xfff7eee000000000, 43,
		6364136223846793005> generator;

	std::normal_distribution<double> longitudinalNormalDistribution;
	std::normal_distribution<double> transverseNormalDistribution;
	std::uniform_real_distribution<double> angleDistribution = std::uniform_real_distribution<double>(0.0, 2 * TMath::Pi());
};


