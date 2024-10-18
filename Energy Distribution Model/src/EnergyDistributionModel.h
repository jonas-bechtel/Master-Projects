#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <filesystem>

#include <TH3D.h>

#include "MCMC.h"
#include "ElectronBeam.h"
#include "IonBeam.h"
#include "Point3D.h"

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

	std::filesystem::path densityFile;
	std::filesystem::path energyFile;

	// additional labelling things
	std::string outputFileName = "";
	std::string label = "";
	//std::filesystem::path descriptionFile;
	// energy range ??
	std::filesystem::path folder;
	int index = 0;
	double driftTubeVoltage = 0;
	double centerLabEnergy = 0;

	//double cathodeVoltage;

	TH1D* operator->()
	{
		return distribution;
	}

	std::string String();

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
	void LoadLabEnergyFile(std::filesystem::path file);
	void GenerateEnergyDistribution();
	void GenerateEnergyDistributionsFromFile(std::filesystem::path file);

private:
	void ShowUI() override;

	void SetupEnergyDistribution();
	void PlotEnergyDistributions();
	void PlotRateCoefficients();
	void PlotCurrentEnergyDistribution();

	void PlotLabEnergyProjections();
	void PLotZweightByEnergy();
	void ClearDistributionList();

private:
	std::vector<EnergyDistribution> energyDistributions;
	EnergyDistribution currentDistribution;
	std::vector<double> currentEnergies;
	
	TGraph* rateCoefficients;
	TH1D* labEnergyProjectionX;
	TH1D* labEnergyProjectionY;
	TH1D* labEnergyProjectionZ;

	TH1D* zPositions;
	TH1D* zWeightByEnergy;

	std::filesystem::path loadedEnergyFile;
	std::filesystem::path currentDescriptionFile;

	// start/end index in description file to generate distribution for
	int startIndex = 1;
	int endIndex = 2;
	bool doAll = false;

	float energyRange[2] = { 6e-0, 30 };
	bool normalise = true;

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


