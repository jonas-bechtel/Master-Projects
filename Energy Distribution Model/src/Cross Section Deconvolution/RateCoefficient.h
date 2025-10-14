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

	void VaryGraphValues();
	void ResetGraphValues();

	void SetLabel(std::string label);
	std::string GetLabel();

	void Convolve(const CrossSection& cs, EnergyDistributionSet& set);

	void Plot(bool showMarkers) const;
	void PlotSubfunctions() const;

	void Clear();
	void Load(std::filesystem::path& file);
	void Save() const;

private:
	int GetIndexOfDetuningEnergy(double Ed) const;
	void SortValuesByDetuningEnergy();

private:
	// main data
	TGraphErrors* graph = new TGraphErrors();

	// ordered in asceding detuning energy
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

