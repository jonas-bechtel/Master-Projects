#pragma once

#include <TCanvas.h>
#include <TH3D.h>

#include "Module.h"
#include "Point3D.h"

class IonBeam : public Module
{
public:
	IonBeam();
	bool HasDensityChanged();
	void CreateIonBeam(TH3D* referenceDensity);

private:
	void ShowUI() override;

private:
	// sigma of gaussian shape in [m]
	float size = 0.006;

	bool densityChanged = false;
};

