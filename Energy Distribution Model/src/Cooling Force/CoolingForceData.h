#pragma once

#include "HeatMapData.h"

struct CoolingForceData
{
	// x,y,z component of cooling force at each position
	TH3D* forceX = nullptr;
	TH3D* forceY = nullptr;
	TH3D* forceZ = nullptr;

	TH3D* positionSamples = nullptr;

	// xy slice of each component
	HeatMapData forceXSlice;
	HeatMapData forceYSlice;
	HeatMapData forceZSlice;

	// sum in z direction of average of all xy planes, represents total energy lost
	double forceXIntegral = 0.0;
	double forceYIntegral = 0.0;
	double forceZIntegral = 0.0;

	// actual resulting value of the force to be compared to measurements (integral divided by length)
	double forceXValue = 0.0;
	double forceYValue = 0.0;
	double forceZValue = 0.0;

	// axes for plotting
	std::vector<double> xAxis;
	std::vector<double> yAxis;
	std::vector<double> zAxis;

	// projections of longitudinal component (z)
	std::vector<double> forceZProjectionX;
	std::vector<double> forceZProjectionY;
	std::vector<double> forceZProjectionZ;

	void SetupHistograms(TH3D* reference);
	void FillData();
	double CalculateIntegral(TH3D* hist);

};

