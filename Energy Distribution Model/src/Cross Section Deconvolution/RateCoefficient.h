#pragma once

class CrossSection;
class EnergyDistributionSet;

class RateCoefficient
{
public:
	RateCoefficient();
	RateCoefficient(const RateCoefficient& other) = delete;
	RateCoefficient& operator=(const RateCoefficient& other) = delete;
	RateCoefficient(RateCoefficient&& other) = default;
	RateCoefficient& operator=(RateCoefficient&& other) = default;

	void SetLabel(std::string label);
	std::string GetLabel();

	void Convolve(const CrossSection& cs, EnergyDistributionSet& set);
	double ConvolveFit(double Ed, double* csBins, const EnergyDistributionSet& set) const;

	void Plot(bool showMarkers) const;
	void PlotSubfunctions() const;

	void Load(std::filesystem::path& file);
	void Save() const;

private:
	int GetIndexOfDetuningEnergy(double Ed) const;

private:
	// main data
	TGraph* graph = new TGraph();

	std::vector<double> detuningEnergies;
	std::vector<double> value;
	std::vector<double> error;
	std::vector<std::vector<double>> psiSubfunctions;

	// labelling things
	std::string label = "mbrc";
	std::filesystem::path energyDistriubtionSetFolder;
	std::filesystem::path crossSectionFile;

	bool measured = true;

	friend class CrossSection;
};

