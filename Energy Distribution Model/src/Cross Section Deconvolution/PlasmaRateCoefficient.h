#pragma once
#include "CrossSection.h"

class PlasmaRateCoefficient
{
public:
	std::string GetLabel();
	void SetLabel(std::string label);

	void Convolve(const CrossSection& cs);
	void ConvolveFromErrorIterationArray(const CrossSection& cs);

	void Plot(bool showMarkers) const;

	void Save() const;
	void Load(std::filesystem::path file);

	static void ShowConvolutionParamterInputs();

private:
	std::vector<double> temperatures;
	std::vector<double> values;
	std::vector<double> errors;

	std::string label = "plasma rate";
	std::filesystem::path crossSectionFile;

	static inline int numberValues = 1000;
	static inline double startTemperature = 10.0; 
	static inline double endTemperature = 50000.0;
};

