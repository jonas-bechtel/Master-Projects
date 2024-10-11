#pragma once
#include <filesystem>

#include <TCanvas.h>
#include <TPolyLine3D.h>
#include <TF1.h>
#include <TVector3.h>
#include <TH3D.h>

#include "Module.h"
#include "Point3D.h"


struct ElectronBeamParameters
{
	double transverse_kT = 2.0e-3;
	//double longitudinal_kT = -3e-4;		// negative value => calculates kTlong automatically from Ie, including LLR and accel. energy spread. Abs(kTlong) still used for resolution estimates - CS binding

	double coolingEnergy = 1; // 33.9342;		// cooling energy [eV]
	double cathodeRadius = 0.0012955;
	double expansionFactor = 30;

	double electronCurrent = 1.2e-08;	// electron current [A]
	double cathodeTemperature = 300;	// cathode temperature [K]
	//double Ecath = 60.528;				// beam energy outside drift tubes
	double LLR = 1.9;					// llr scaling factor
	double sigmaLabEnergy = 0.0;		// sigma of gausian spread of the acceleration voltage (Elab) [eV]
	double extractionEnergy = 31.26;	// extraction voltage for LLR (<=0 -> use Elab)
};

class ElectronBeam : public Module
{
public:
	ElectronBeam();
	ElectronBeamParameters& GetParamters();
	void LoadDensityFile(std::filesystem::path file);
	bool HasDensityChanged();

	TVector3 GetDirection(double z);
	TVector3 GetDirection(Point3D point);
	TVector3 GetNormal(double z);
	double GetLongitudinal_kT(double labEnergy);
	double GetTransverse_kT();

private:
	void ShowUI();

	void PlotTrajectory();
	double Trajectory(double z);
	double Derivative(double z);

	double DistancePointToTrajectoryOfZ(double z, Point3D point);

private:
	ElectronBeamParameters parameters;

	float sliderZ = 0;
	float sliderY = 0;
	bool densityChanged = false;
};



