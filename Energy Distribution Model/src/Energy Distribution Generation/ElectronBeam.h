#pragma once

#include "Module.h"
#include "Point3D.h"

class ElectronBeam : public EnergyDistributionModule
{
public:
	ElectronBeam();
	void SetupDistribution(std::filesystem::path densityfile = "") override;
	TH3D* GetDistribution() override;

	void CalculateDetuningEnergy();
	
	void LoadDensityFile(std::filesystem::path file);
	void GenerateElectronBeamDensity();

	double Trajectory(double z) const;

	TVector3 GetDirection(double z);
	TVector3 GetDirection(Point3D point);
	double GetLongitudinal_kT(double labEnergy);
	double GetTransverse_kT();

private:
	void ShowUI() override;

	void CutZerosFromDistribution();
	void CreateLargeDistribution();

	TVector3 GetNormal(double z);
	
	double Derivative(double z) const;
	double DistancePointToTrajectoryOfZ(double z, Point3D point);

	void PlotDensitySlice();
	void PlotGeneratedDensities();
	void PlotTrajectory();
	void PlotProjections();

private:
	ElectronBeamParameters& m_parameters;
	TProfile2D* bla = nullptr;

	TH3D* generatedBeamDensity = nullptr;
	TH3D* generatedBeamDensitySmall = nullptr;
	TH2D* densitySliceXY = nullptr;

	TH1D* electronBeamProjectionX = nullptr;
	TH1D* electronBeamProjectionY = nullptr;
	TH1D* electronBeamProjectionZ = nullptr;

	TH1D* generatedBeamProjectionX = nullptr;
	TH1D* generatedBeamProjectionY = nullptr;
	TH1D* generatedBeamProjectionZ = nullptr;

	// parameters to increase histogram resolution by interpolation
	bool increaseHist = false;
	int factor = 3;

	float sliderZ = 0;
	float sliderY = 0;

	float SliceZ = 0.0f;
};



