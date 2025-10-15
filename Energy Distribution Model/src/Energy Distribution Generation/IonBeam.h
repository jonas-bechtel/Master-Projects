#pragma once

#include "Point3D.h"
#include "HeatMapData.h"
#include "ParameterImplementations.h"

namespace IonBeam
{
	IonBeamParameters GetParameters();
	void SetParameters(const IonBeamParameters& params);

	void Init();

	void CreateFromReference(TH3D* reference);
	void UpdateHistData();

	TVector3 GetDirection();
	TVector3 GetVelocity();
	double GetValue(double x, double y, double z);
	double GetVelocityMagnitude();
	float GetSigmaX();
	float GetSigmaY();
	float* GetLimitedRange();
	bool IsRangeLimited();
	TH3D* Get();
	std::string GetTags();

	void ShowWindow();
	void ShowPlots();
	void ShowParameterControls();
	void ShowCoolingForceParameterControls();
}
