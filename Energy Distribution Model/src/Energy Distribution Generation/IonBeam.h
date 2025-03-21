#pragma once

#include "Point3D.h"
#include "HeatMapData.h"
#include "ParameterImplementations.h"

namespace IonBeam
{
	IonBeamParameters GetParameters();
	void Init();

	void CreateFromReference(TH3D* reference);
	void UpdateMainData();
	void UpdatePlotData();

	TVector3 GetDirection();
	int GetCharge();
	TH3D* Get();
	std::string GetTags();

	void ShowWindow();
	void ShowPlots();
	void ShowParameterControls();
	void ShowCoolingForceParameterControls();

	std::vector<Point3D> GeneratePositions();
}
