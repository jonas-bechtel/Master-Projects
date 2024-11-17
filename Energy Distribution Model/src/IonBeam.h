#pragma once

#include "Module.h"
#include "Point3D.h"

struct IonBeamParameters : public Parameters
{
	IonBeamParameters()
	{
		setName("ion beam parameters");
	}

	ParameterValue<float2> shift = ParameterValue(float2(0.0f, 0.0f), "shift in x and y", "%.4f, %.4f m");
	ParameterValue<double> amplitude1 = ParameterValue(10.1, "amplitude 1", "%.4f");
	ParameterValue<double> amplitude2 = ParameterValue(8.1, "amplitude 2", "%.4f");
	ParameterValue<float2> shape1 = ParameterValue(float2(9.5e-3f, 5.7e-3f), "shape 1 (x,y)", "%.4f, %.4f m");
	ParameterValue<float2> shape2 = ParameterValue(float2(1.39e-3f, 2.15e-3f), "shape 2 (x,y)", "%.4f, %.4f m");

	ParameterValue<bool> useSingleGaussian = ParameterValue(false, "using single gaussian", "%d", true);
	ParameterValue<double> radius = ParameterValue(0.0010, "radius", "%.4f m", true);

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
	void SetupDistribution(std::filesystem::path file = "") override;
	IonBeamParameters& GetParameter();
	TH3D* MultiplyWithElectronDensities();

	void PlotIonBeamProjections();

private:
	void ShowUI() override;

private:
	IonBeamParameters m_parameters;

	TH1D* ionBeamProjectionX = nullptr;
	TH1D* ionBeamProjectionY = nullptr;
	TH1D* ionBeamProjectionZ = nullptr;
};

