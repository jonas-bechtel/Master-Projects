#pragma once

struct PlasmaRateCoefficient
{
	std::vector<double> temperature;
	std::vector<double> value;
	std::vector<double> error;

	std::string label = "";
	std::filesystem::path crossSectionFile;
};

