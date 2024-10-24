#pragma once
#include "Module.h"

#include <filesystem>

struct LabEnergiesParameters
{
	double centerLabEnergy = 0;
	std::filesystem::path energyFile;

	// parameters for simpler model
	bool useUniformEnergies = false;
	bool useOnlySliceXY = false;
	float sliceToFill = 0.5;

	std::string String();
};

class LabEnergies : public Distribution3D
{
public:
	LabEnergies();
	double Get(double x, double y, double z);
	LabEnergiesParameters GetParameter();
	void SetCenterLabEnergy(double energy);
	void SetupDistribution(std::filesystem::path energyfile = "") override;

private:
	void ShowUI() override;
	void LoadLabEnergyFile(std::filesystem::path file);

	void GenerateUniformLabEnergy();
	void FillEnergiesWithXY_Slice();

	void PlotLabEnergySlice();
	void PlotLabEnergyProjections();

private:
	LabEnergiesParameters parameter;

	TH1D* labEnergyProjectionX = nullptr;
	TH1D* labEnergyProjectionY = nullptr;
	TH1D* labEnergyProjectionZ = nullptr;
	TH2D* labEnergySliceXY = nullptr;

	// z value for the xy slice of the lab energies
	float SliceZ = 0.0f;
};

