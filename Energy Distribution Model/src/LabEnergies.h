#pragma once

#include "Module.h"
#include "Parameter.h"

struct LabEnergyParameters : public Parameters
{
	LabEnergyParameters() { setName("lab energy parameters"); }

	ParameterValue<double> centerLabEnergy = ParameterValue(0.0, "lab energy in center", "%.3e eV");
	ParameterValue<double> driftTubeVoltage = ParameterValue(0.0, "drift tube voltage", "%e V");
	ParameterValue<Path> energyFile = ParameterValue(Path(""), "energy file", "%s");

	ParameterValue<bool> useUniformEnergies = ParameterValue(false, "use uniform energy", "%d", true);
	ParameterValue<bool> useOnlySliceXY = ParameterValue(false, "use slice of energies", "%d", true);
	ParameterValue<double> sliceToFill = ParameterValue(0.5, "z value of slice", "%.3f", true);

private:
	int GetSize() override
	{
		return sizeof(*this);
	}
};

class LabEnergies : public Distribution3D
{
public:
	LabEnergies();
	double Get(double x, double y, double z);
	LabEnergyParameters& GetParameter();
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
	LabEnergyParameters m_parameters;

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

