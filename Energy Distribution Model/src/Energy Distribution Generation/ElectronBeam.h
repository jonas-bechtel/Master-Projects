#pragma once

#include "Point3D.h"

#include "HistData3D.h"
#include "ParameterImplementations.h"

namespace ElectronBeam
{
	void Init();

	void SetupDistribution(std::filesystem::path densityfile);

	TH3D* Get();
	ElectronBeamParameters GetParameters();
	void SetParameters(const ElectronBeamParameters& params);

	TVector3 GetDirection(double z);
	TVector3 GetVelocity(double z, double energy);
	double GetVelocityMagnitude(double energy);
	double GetLongitudinal_kT(double labEnergy);
	double GetTransverse_kT();
	double GetCoolingEnergy();
	double GetDensity(const Point3D& point);

	void SetElectronCurrent(double current);

	std::string GetTags();

	void CalculateDetuningEnergy();
	void CalculateDetuningVelocity();
	void CalculateEstimateLongkT();

	void ShowWindow();
	void ShowList();
	void ShowPlots();
	void ShowParameterControls();
	void PlotTrajectory();

	TH3D* LoadDensityFile(std::filesystem::path file);
	TH3D* GenerateElectronBeamDensity();

	HistData3D* GetSelected();
	std::string GetSelectedBeamLabel();
	void SelectedItemChanged();
	void AddBeamToList(HistData3D& eBeam);
	void RemoveBeamFromList(int index);

	TH3D* CutZerosFromDistribution(TH3D* input);
	TH3D* MirrorDistributionAtZ(TH3D* input);
	TH3D* CreateLargeDistribution(TH3D* input);

	double Trajectory(double z);
	TVector3 GetNormal(double z);
	double Derivative(double z);
	//double DistancePointToTrajectoryOfZ(double z, Point3D point);
}




