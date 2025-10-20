#pragma once
#include "Point3D.h"
#include "HeatMapData.h"
#include "ROOTCanvas.h"
#include "HistData3D.h"
#include "ParameterImplementations.h"

namespace LabEnergy
{
	void Init();

	TH3D* Get();
	double GetValue(double x, double y, double z);
	LabEnergyParameters GetParameters();
	void SetParameters(const LabEnergyParameters& params);

	double GetCenterLabEnergy();
	double GetOutsideDriftTubeEnergy();
	void SetDriftTubeVoltage(double voltage);
	void SetCenterEnergy(double energy);

	std::string GetTags();

	void SetupDistribution(std::filesystem::path energyfile);

	TH3D* LoadLabEnergyFile(std::filesystem::path file);
	TH3D* GenerateUniformLabEnergies();

	void SelectedItemChanged();
	void AddBeamToList(HistData3D& beamData);
	void RemoveBeamFromList(int index);

	void ShowWindow();
	void ShowList();
	void ShowParameterControls();
	void ShowPlots();
}

