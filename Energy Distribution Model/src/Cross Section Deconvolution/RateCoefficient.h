#pragma once
#include "Module.h"

struct RateCoefficient
{
	RateCoefficient();


public:
	// main data
	TGraph* graph = new TGraph();

	std::vector<double> detuningEnergies;
	std::vector<double> value;
	std::vector<double> error;

	// labelling things
	std::string label = "";
	std::filesystem::path energyDistriubtionSetFolder;
	std::filesystem::path crossSectionFile;
};

