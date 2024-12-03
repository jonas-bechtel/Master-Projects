#pragma once

#include "Module.h"
#include "Point3D.h"

class IonBeam : public Distribution3D
{
public:
	IonBeam();
	void SetupDistribution(std::filesystem::path file = "") override;
	TH3D* MultiplyWithElectronDensities();

	void PlotIonBeamProjections();

private:
	void ShowUI() override;

private:
	IonBeamParameters& m_parameters;

	TH1D* ionBeamProjectionX = nullptr;
	TH1D* ionBeamProjectionY = nullptr;
	TH1D* ionBeamProjectionZ = nullptr;
};

