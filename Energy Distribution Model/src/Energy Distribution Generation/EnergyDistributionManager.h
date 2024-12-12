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

private:
	void SetupDistribution(std::filesystem::path = "") override {}

	void ShowUI() override;
	void ShowSettings();
	void ShowTabsWithSets();
	void ShowEnergyDistributionSet(int setIndex);
	void ShowEnergyDistributionPlot();

	void CreateNewSet();
	void RemoveSet(int setIndex);
	void AddDistributionToSet(EnergyDistribution&& distribution, int setIndex);
	void RemoveDistributionFromSet(int index, int setIndex);

	void SetupSecondaryPlots();
	void PlotEnergyDistributions();
	void PLotZweightByEnergy();
	void PlotLongkTDistribution();
	void PlotLongVelAddition();

	void ClearDistributionsInSet(int setIndex);

private:	
	// graphs and plots
	TH1D* zPositions = nullptr;
	TH1D* zWeightByEnergy = nullptr;
	TH1D* long_ktDistribution = nullptr;
	TH1D* long_VelAddition = nullptr;

	int currentSetIndex = 0;

	// currently loaded file
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


