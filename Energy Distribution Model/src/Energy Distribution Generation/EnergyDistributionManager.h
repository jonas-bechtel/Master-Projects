#pragma once

#include "MCMC.h"
#include "ElectronBeam.h"
#include "IonBeam.h"
#include "LabEnergies.h"
#include "EnergyDistribution.h"

class EnergyDistributionManager : public EnergyDistributionModule, public EnergyDistribtionSetsContainer
{
public:
	EnergyDistributionManager();
	std::vector<EnergyDistributionSet>& GetEnergyDistributionSets();
	void GenerateEnergyDistribution();
	void GenerateEnergyDistributionsFromFile(std::filesystem::path file);

	void ShowSetListWindow();
private:
	void SetupDistribution(std::filesystem::path = "") override {}

	void ShowUI() override;
	void ShowSettings();
	void ShowTabsWithSets();
	void ShowEnergyDistributionSet(int setIndex);
	void ShowEnergyDistributionPlot();
	void ShowSetInformationWindow();
	void ShowAnalyticalParameterWindow();

	void CreateNewSet();
	void RemoveSet(int setIndex);
	void AddDistributionToSet(EnergyDistribution&& distribution, int setIndex);
	void RemoveDistributionFromSet(int index, int setIndex);
	void PrepareCurrentSet(std::filesystem::path folder, std::filesystem::path subfolder = "");

	void SetupSecondaryPlots();
	void PlotEnergyDistributions();
	void PLotZweightByEnergy();
	void PlotLongkTDistribution();
	void PlotLongVelAddition();

	void UpdateAnalytical();

	void ClearDistributionsInSet(int setIndex);

private:	
	// graphs and plots
	TH1D* zPositions = nullptr;
	TH1D* zWeightByEnergy = nullptr;
	TH1D* long_ktDistribution = nullptr;
	TH1D* long_VelAddition = nullptr;

	// currently loaded file
	std::filesystem::path currentDescriptionFile = std::filesystem::path("data\\dataset1\\100x100x100_Ie0.95_Ucath44.2_RelTol0_mbrc1_energies.asc");
	int maxIndex = 0;

	// start/end index in description file to generate distribution for
	int startIndex = 1;
	int endIndex = 1;
	bool doAll = false;

	// save all sampled values in a file to load it again
	bool loadSamples = true;

	// set information window things
	bool showSetInformation = false;
	bool infoPlotsLogX = false;
	bool infoPlotsLogY = false;

	// binning related parameters
	BinningSettings binSettings;

	// fit options
	bool fixKT_trans = true;
	bool fixDetuningEnergy = false;
	bool showLimitedFit = true;
	
	// plot parameters
	bool logX = true;
	bool logY = true;
	bool showMarkers = false;
	bool showFits = true;
	bool showAnalytical = false;

	// Analytical Parameter
	float scale = 1;
	float E_d = 10;
	float kT_long = 1e-4;
	float kT_trans = 0.002;
	float energyRange[2] = { 0,1 };
	float energies[200];
	float values[200];
	float color[3] = { 0.6, 0.2, 0.4};

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


