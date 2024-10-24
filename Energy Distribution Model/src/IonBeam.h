#pragma once

#include <TCanvas.h>
#include <TH3D.h>

#include "Module.h"
#include "Point3D.h"

struct IonBeamParameters
{
	// sigma of gaussian shape in [m]
	float radius = 0.0010f;
	float shift[2] = { 0, 0 };

	std::string String();
};

class IonBeam : public Distribution3D
{
public:
	IonBeam();
	float Radius();
	IonBeamParameters GetParameter();
	void SetParameter(IonBeamParameters params);
	TH3D* MultiplyWithElectronDensities(TH3D* electronDensities);

private:
	void ShowUI() override;

	void CreateIonBeam(TH3D* referenceDensity);
private:
	IonBeamParameters parameter;
};

