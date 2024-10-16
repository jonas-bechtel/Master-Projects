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
	// actual distribution
	TH1D* distribution;

	// all the things used to create it
	MCMC_Parameters mcmcParameter;
	ElectronBeamParameters eBeamParameter;
	IonBeamParameters ionBeamParameter;

	std::filesystem::path densityFile;
	std::filesystem::path energyFile;

	// additional labelling things
	std::string name = "";
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
	void ShowUI();

	void SetupEnergyDistribution();
	void PlotEnergyDistributions();
	void PlotCurrentEnergyDistribution();

	void PlotLabEnergyProjections();
	void ClearDistributionList();

private:
	std::vector<EnergyDistribution> energyDistributions;
	EnergyDistribution currentDistribution;
	
	TH1D* labEnergyProjectionX;
	TH1D* labEnergyProjectionY;
	TH1D* labEnergyProjectionZ;

	std::filesystem::path loadedEnergyFile;

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


