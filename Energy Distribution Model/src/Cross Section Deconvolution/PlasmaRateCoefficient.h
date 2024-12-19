#pragma once

struct PlasmaRateCoefficient
{
	std::vector<double> temperatures;
	std::vector<double> values;
	std::vector<double> errors;

	std::string label = "";
	std::filesystem::path file = "";
	std::filesystem::path crossSectionFile;
};

