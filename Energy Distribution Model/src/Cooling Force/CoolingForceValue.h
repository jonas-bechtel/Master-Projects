#pragma once

#include "HeatMapData.h"
#include "ElectronBeam.h"
#include "IonBeam.h"
#include "LabEnergies.h"

using RNG_engine = std::mersenne_twister_engine<std::uint_fast64_t,
	64, 312, 156, 31,
	0xb5026f5aa96619e9, 29,
	0x5555555555555555, 17,
	0x71d67fffeda60000, 37,
	0xfff7eee000000000, 43,
	6364136223846793005>;

struct CoolingForceValue
{
public:
	CoolingForceValue();
	~CoolingForceValue();

	CoolingForceValue(const CoolingForceValue& other) = delete;
	CoolingForceValue& operator=(const CoolingForceValue& other) = delete;
	CoolingForceValue(CoolingForceValue&& other) noexcept;
	CoolingForceValue& operator=(CoolingForceValue&& other) noexcept;

	void Calculate(std::filesystem::path descriptionFile, int index);

	bool ShowListItem(bool selected) const;

	void SetupHistograms(TH3D* reference);
	void FillData();
	double CalculateIntegral(TH3D* hist);

private:
	void CopyParameters();
	void SetupLabellingThings();
	void ResetDefaultValues();

	std::string GetHeaderString() const;
	std::string Filename() const;

private:
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

	// all the parameters used to create it
	ElectronBeamParameters eBeamParameter;
	IonBeamParameters ionBeamParameter;
	LabEnergyParameters labEnergiesParameter;

	// additional labelling things
	std::string label = "";
	std::string tags = "";
	int index = 0;

	// random number generation things
	static std::mersenne_twister_engine<std::uint_fast64_t,
		64, 312, 156, 31,
		0xb5026f5aa96619e9, 29,
		0x5555555555555555, 17,
		0x71d67fffeda60000, 37,
		0xfff7eee000000000, 43,
		6364136223846793005> generator;

	static std::normal_distribution<double> longitudinalNormalDistribution;
	static std::normal_distribution<double> transverseNormalDistribution;

	friend class CoolingForceCurve;
};

