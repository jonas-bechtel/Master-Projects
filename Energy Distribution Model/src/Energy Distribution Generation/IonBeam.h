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
	void UpdateMainData();
	void UpdatePlotData();

	TVector3 GetDirection();
	TVector3 GetVelocity();
	double GetValue(const Point3D& point);
	double GetVelocityMagnitude();
	int GetCharge();
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

	std::vector<Point3D> GeneratePositions();
}
