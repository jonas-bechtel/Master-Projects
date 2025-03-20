#pragma once
#include "CrossSection.h"

class PlasmaRateCoefficient
{
public:
	std::string GetLabel();
	void SetLabel(std::string label);

	void Convolve(const CrossSection& cs);

	void Plot(bool showMarkers) const;

	void Save() const;
	void Load(std::filesystem::path file);

private:
	std::vector<double> temperatures;
	std::vector<double> values;
	std::vector<double> errors;

	std::string label = "plasma rate";
	std::filesystem::path crossSectionFile;
};

