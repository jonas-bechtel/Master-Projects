#pragma once

#include "MCMC.h"
#include "ElectronBeam.h"
#include "IonBeam.h"
#include "LabEnergies.h"
#include "EnergyDistribution.h"

struct BinningSettings
{
	float energyRange[2] = { 1e-5, 100 };

	bool constantBinSize = true;
	double normalStepSize = 0.5;
	double peakStepSize = 0.03;

	bool factorBinning = false;
	int binsPerDecade = 200;

	bool increasePeakResolution = true;
};

class EnergyDistributionManager : public Module
{
public:
	EnergyDistributionManager();
	EnergyDistributionParameters& GetParameter();
	std::vector<EnergyDistribution*>& GetEnergyDistributions();
	void GenerateEnergyDistribution();
	void GenerateEnergyDistributionsFromFile(std::filesystem::path file);

private:
	void ShowUI() override;
	void ShowSettings();
	void ShowEnergyDistributionList();
	void ShowEnergyDistributionPlot();

	void SetupEnergyDistribution(EnergyDistribution* distribution);
	void AddDistributionToList(EnergyDistribution* distribution);
	void RemoveDistributionFromList(int index);

	// functions for the analytical model
	void FitAnalyticalToGeneratedDistribution(EnergyDistribution* distribution);
	void GenerateAnalyticalDistribution();
	double AnalyticalEnergyDistributionFit(double* x, double* params);
	double AnalyticalEnergyDistribution(double* x, double* params);
	double AnalyticalEnergyDistribution(double Ecm, double Ed, double Ttr, double Tlong);
	double ComplexErrorFunction(double* x, double* par);
	double DawsonIntegral(double* x, double* par);
	double ExpDiff(double* x, double* par);

	void PlotEnergyDistributions();
	void PLotZweightByEnergy();
	void PlotLongkTDistribution();
	void PlotLongVelAddition();

	void ClearDistributionList();

private:
	std::vector<EnergyDistribution*> energyDistributions;
	EnergyDistribution* currentDistribution;// = new EnergyDistribution();
	EnergyDistributionParameters parameter;
	
	// graphs and plots
	TH1D* zPositions = nullptr;
	TH1D* zWeightByEnergy = nullptr;
	TH1D* long_ktDistribution = nullptr;
	TH1D* long_VelAddition = nullptr;

	// currently loaded files
	std::filesystem::path currentDescriptionFile = std::filesystem::path("data\\C60 (2)\\100x100x100_Ie11.3_Ucath44.2_RelTol0_sort_energies.asc");
	int maxIndex = 0;

	// start/end index in description file to generate distribution for
	int startIndex = 1;
	int endIndex = 1;
	bool doAll = false;

	// save all sampled values in a file to load it again
	bool saveAsHist = false;
	bool saveSamplesToFile = false;
	bool loadSamples = true;

	// binning related parameters
	BinningSettings binSettings;
	
	// parameters for analytical energy distribution generation
	AnalyticalDistributionParameters analyticalParameter;
	float analyticalEnergyRange[2] = { 1e-1, 10 };
	int analyticalNumberBins = 1000;
	
	// plot parameters
	bool logX = true;
	bool logY = true;
	bool showMarkers = false;

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


