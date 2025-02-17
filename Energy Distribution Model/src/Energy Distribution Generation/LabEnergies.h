#pragma once

#include "Module.h"
#include "HeatMapData.h"

struct LabEnergy
{
	TH3D* fullHistogram = nullptr;
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

	void FillData(const ElectronBeamWindow* eBeam);

	LabEnergy() {}
	LabEnergy(const LabEnergy& other) = delete;
	LabEnergy& operator=(const LabEnergy& other) = delete;
	LabEnergy(LabEnergy&& other) noexcept;
	LabEnergy& operator=(LabEnergy&& other) noexcept;
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

	void RemoveBeamFromList(int index);

	void GenerateUniformLabEnergy();
	void FillEnergiesWithXY_Slice();

private:
	LabEnergyParameters& m_parameters;

	std::vector<LabEnergy> labEnergiesToLookAt;
	int selectedIndex = -1;

	// z value for the xy slice of the lab energies
	float SliceZ = 0.0f;

};

