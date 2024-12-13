#pragma once
#include "Module.h"

struct RateCoefficient
{
	RateCoefficient();
	RateCoefficient(const RateCoefficient& other) = delete;
	RateCoefficient& operator=(const RateCoefficient& other) = delete;
	RateCoefficient(RateCoefficient&& other) = default;
	RateCoefficient& operator=(RateCoefficient&& other) = default;

	int GetIndexOfDetuningEnergy(double Ed);
	

public:
	// main data
	TGraph* graph = new TGraph();

	std::vector<double> detuningEnergies;
	std::vector<double> value;
	std::vector<double> error;

	// labelling things
	std::string label = "mbrc";
	std::filesystem::path energyDistriubtionSetFolder;
	std::filesystem::path crossSectionFile;

	bool input = true;
};

