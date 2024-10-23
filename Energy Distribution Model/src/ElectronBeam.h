#pragma once
#include <filesystem>

#include <TCanvas.h>
#include <TF1.h>
#include <TVector3.h>
#include <TH3D.h>
#include <TProfile2D.h>

#include "Module.h"
#include "Point3D.h"


struct ElectronBeamParameters
{
	double transverse_kT = 2.0e-3;
	//double longitudinal_kT = -3e-4;		// negative value => calculates kTlong automatically from Ie, including LLR and accel. energy spread. Abs(kTlong) still used for resolution estimates - CS binding

	double coolingEnergy = 0.15; 		// cooling energy [eV]
	double cathodeRadius = 0.0012955;
	double expansionFactor = 30;

	double electronCurrent = 1.2e-08;	// electron current [A]
	double cathodeTemperature = 300;	// cathode temperature [K]
	//double Ecath = 60.528;				// beam energy outside drift tubes
	double LLR = 1.9;					// llr scaling factor
	double sigmaLabEnergy = 0.0;		// sigma of gausian spread of the acceleration voltage (Elab) [eV]
	double extractionEnergy = 31.26;	// extraction voltage for LLR (<=0 -> use Elab)

	std::filesystem::path densityFile;

	// parameters for analytical electron beam shape to test model
	bool hasGaussianShape = false;
	bool hasNoBending = false;
	bool hasFixedLongitudinalTemperature = false;
	double radius = 0.003;
	double longitudinal_kT = 3e-4;

	std::string String();
};

class ElectronBeam : public Module
{
public:
	ElectronBeam();
	ElectronBeamParameters GetParameter();
	TH3D* GetDistribution() override;
	void SetParameter(ElectronBeamParameters params);
	void SetCurrent(double current);

	void LoadDensityFile(std::filesystem::path file);
	void GenerateElectronBeamDensity();

	TVector3 GetDirection(double z);
	TVector3 GetDirection(Point3D point);
	double GetLongitudinal_kT(double labEnergy);
	double GetTransverse_kT();

private:
	void ShowUI() override;

	void CutZerosFromDistribution();
	void CreateLargeDistribution();

	TVector3 GetNormal(double z);
	double Trajectory(double z);
	double Derivative(double z);
	double DistancePointToTrajectoryOfZ(double z, Point3D point);

	void PlotDensitySlice();
	void PlotGeneratedDensities();
	void PlotTrajectory();
	void PlotProjections();

private:
	ElectronBeamParameters parameters;

	TH3D* generatedBeamDensity = nullptr;
	TH3D* generatedBeamDensitySmall = nullptr;
	TH2D* densitySliceXY = nullptr;

	TH1D* electronBeamProjectionX = nullptr;
	TH1D* electronBeamProjectionY = nullptr;
	TH1D* electronBeamProjectionZ = nullptr;

	TH1D* generatedBeamProjectionX = nullptr;
	TH1D* generatedBeamProjectionY = nullptr;
	TH1D* generatedBeamProjectionZ = nullptr;

	//TProfile2D* electronBeamProfileXY = nullptr;
	//TProfile2D* electronBeamProfileXZ = nullptr;
	//TProfile2D* electronBeamProfileYZ = nullptr;

	// parameters to increase histogram resolution by interpolation
	bool increaseHist = false;
	int factor = 3;

	float sliderZ = 0;
	float sliderY = 0;

	float SliceZ = 0.0f;
};



