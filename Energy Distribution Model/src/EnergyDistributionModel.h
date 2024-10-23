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

struct EnergyDistributionParameters
{
	double driftTubeVoltage = 0;
	double centerLabEnergy = 0;

	std::filesystem::path energyFile;

	// parameters for simpler test 
	bool useUniformEnergies = false;
	bool useOnlySliceXY = false;
	float sliceToFill = 0.5;
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
	void LoadLabEnergyFile(std::filesystem::path file);
	void GenerateEnergyDistribution();
	void GenerateEnergyDistributionsFromFile(std::filesystem::path file);

private:
	void ShowUI() override;
	void ShowPlots() override;
	void ShowEnergyDistributionList();

	void GenerateUniformLabEnergy(float energy);
	void FillEnergiesWithXY_Slice();

	void SetupEnergyDistribution();

	void PlotLabEnergySlice();
	void PlotEnergyDistributions();
	void PlotCurrentEnergyDistribution();
	void PlotLabEnergyProjections();
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
	TH1D* labEnergyProjectionX = nullptr;
	TH1D* labEnergyProjectionY = nullptr;
	TH1D* labEnergyProjectionZ = nullptr;
	TH2D* labEnergySliceXY = nullptr;

	TH1D* zPositions = nullptr;
	TH1D* zWeightByEnergy = nullptr;
	TH1D* long_ktDistribution = nullptr;
	TH1D* long_VelAddition = nullptr;

	// currently loaded files
	std::filesystem::path currentDescriptionFile;
	int maxIndex = 0;

	// start/end index in description file to generate distribution for
	int startIndex = 1;
	int endIndex = 1;
	bool doAll = false;

	// parameters for energy distribution generation
	float energyRange[2] = { 6e-0, 100 };
	bool normalise = true;

	// list things
	//int currentItem = -1;

	// plot parameters
	bool logX = true;
	bool logY = true;

	// z value for the xy slice of the lab energies
	float SliceZ = 0.0f;

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


