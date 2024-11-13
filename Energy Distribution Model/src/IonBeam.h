#pragma once

#include <TCanvas.h>
#include <TH3D.h>

#include "Module.h"
#include "Point3D.h"

struct IonBeamParameters : public Parameters
{
	IonBeamParameters()
	{
		setName("ion beam parameters");
	}

	ParameterValue<double> radius = ParameterValue(0.0010, "radius", "%.4f m");
	ParameterValue<float2> shift = ParameterValue(float2(0.0f, 0.0f), "shift in x and y", "%.4f, %.4f m");

private:
	int GetSize() override
	{
		return sizeof(*this);
	}
};

class IonBeam : public Distribution3D
{
public:
	IonBeam();
	double GetRadius();
	void SetupDistribution(std::filesystem::path file = "") override;
	IonBeamParameters GetParameter();
	TH3D* MultiplyWithElectronDensities();

private:
	void ShowUI() override;

private:
	IonBeamParameters m_parameters;
};

