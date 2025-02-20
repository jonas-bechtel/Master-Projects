#pragma once

#include "Module.h"
#include "Point3D.h"
#include "HeatMapData.h"

class IonBeamWindow : public EnergyDistributionModule
{
public:
	IonBeamWindow();
	void SetupDistribution(std::filesystem::path file = "") override;
	TH3D* MultiplyWithElectronDensities();

private:
	void ShowUI() override;
	void ShowSettings();
	void ShowPlots();

	void FillHistogram(TH3D* hist);
	void UpdateDataToLookAt();

private:
	IonBeamParameters& m_parameters;
	bool useSecondGaus = false;

	// ion beam to look at
	std::vector<double> xAxis;
	std::vector<double> yAxis;
	std::vector<double> zAxis;

	std::vector<double> projectionValuesX;
	std::vector<double> projectionValuesY;
	std::vector<double> projectionValuesZ;

	HeatMapData slice;
	float SliceZ = 0.0f;
};

