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

class EnergyDistributionModel : public Module
{
public:
	EnergyDistributionModel();

private:
	void ShowUI();

	TH3D* MultiplyElectronAndIonDensities(TH3D* electronDensities, TH3D* ionDensities);
	void LoadLabEnergyFile(std::filesystem::path file);
	void CreateEnergyDistributionHistogram();
	void GenerateEnergyDistribution();
	void PlotEnergyDistribution();

	void PlotLabEnergyProjections();

private:
	ElectronBeam electronBeam;
	IonBeam ionBeam;
	MCMC mcmcSampler;

	TH1D* energyDistribution;
	
	TH1D* labEnergyProjectionX;
	TH1D* labEnergyProjectionY;
	TH1D* labEnergyProjectionZ;

	// sigma of gaussian used for iona beam, in [m]
	double ionBeamRadius = 0.006;
	float energyRange[2] = { 1e-7, 10 };

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


