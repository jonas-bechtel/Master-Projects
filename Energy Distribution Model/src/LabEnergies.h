#pragma once

#include "Module.h"

class LabEnergies : public EnergyDistributionModule
{
public:
	LabEnergies();
	double Get(double x, double y, double z);
	void SetupDistribution(std::filesystem::path energyfile = "") override;

private:
	void ShowUI() override;
	void LoadLabEnergyFile(std::filesystem::path file);

	void GenerateUniformLabEnergy();
	void FillEnergiesWithXY_Slice();

	void PlotLabEnergySlice();
	void PlotLabEnergyProjections();
	void PlotOutInsideEnergyOnZ();

private:
	LabEnergyParameters& m_parameters;

	TH1D* labEnergyProjectionX = nullptr;
	TH1D* labEnergyProjectionY = nullptr;
	TH1D* labEnergyProjectionZ = nullptr;
	TH2D* labEnergySliceXY = nullptr;
	TGraph* labEnergyInside = nullptr;
	TGraph* labEnergyOutside = nullptr;

	std::vector<double> zValues;
	std::vector<std::vector<double>> energyValuesInside;
	std::vector<std::vector<double>> energyValuesOutside;
	int counter = 0;

	// z value for the xy slice of the lab energies
	float SliceZ = 0.0f;

};

