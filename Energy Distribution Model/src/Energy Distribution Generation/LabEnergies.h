#pragma once

#include "Module.h"


struct HeatMapData
{
	std::vector<double> values;
	int nRows;
	int nCols;
	double minValue;
	double maxValue;
	ImPlotPoint bottomLeft;
	ImPlotPoint topRight;

	void FromTH3D(TH3D* hist, float zSliceValue);
};

struct LabEnergy
{
	TH3D* fullHistogram;
	std::vector<double> xAxis;
	std::vector<double> yAxis;
	std::vector<double> zAxis;

	std::vector<double> projectionValuesX;
	std::vector<double> projectionValuesY;
	std::vector<double> projectionValuesZ;

	HeatMapData slice;

	std::vector<double> labEnergyInside;
	std::vector<double> labEnergyOutside;

	std::string label;

	void FillData(const ElectronBeam* eBeam);

	LabEnergy() {}
	LabEnergy(const LabEnergy& other) = delete;
	LabEnergy& operator=(const LabEnergy& other) = delete;
	LabEnergy(LabEnergy&& other);
	LabEnergy& operator=(LabEnergy&& other);
	~LabEnergy()
	{
		//std::cout << "deleting " << fullHistogram << std::endl;
		delete fullHistogram;
	}
};

class LabEnergyWindow : public EnergyDistributionModule
{
public:
	LabEnergyWindow();
	double Get(double x, double y, double z);
	void SetupDistribution(std::filesystem::path energyfile = "") override;

private:
	void ShowUI() override;
	void ShowList();
	void ShowSettings();
	void ShowLabEnergyPlots();

	void LoadLabEnergyFile(std::filesystem::path file);
	void LoadToLookAt(std::filesystem::path file);

	void GenerateUniformLabEnergy();
	void FillEnergiesWithXY_Slice();

private:
	LabEnergyParameters& m_parameters;

	std::vector<LabEnergy> labEnergiesToLookAt;
	int selectedIndex = -1;

	// z value for the xy slice of the lab energies
	float SliceZ = 0.0f;

};

