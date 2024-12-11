#pragma once

struct CrossSection
{
	// main data
	std::vector<double> energy;
	std::vector<double> value;
	std::vector<double> error;

	// labelling things
	std::string label = "";
	std::filesystem::path energyDistriubtionSetFolder;
	std::filesystem::path mergedBeamRateCoefficientFile;
};

