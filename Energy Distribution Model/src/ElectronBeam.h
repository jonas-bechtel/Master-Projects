#pragma once
#include <filesystem>

#include <TCanvas.h>
#include <TF1.h>
#include <TVector3.h>
#include <TH3D.h>
#include <TProfile2D.h>

#include "Module.h"
#include "Point3D.h"

struct ElectronBeamParameters : public Parameters
{
	ElectronBeamParameters()
	{
		setName("electron beam parameters");
	}

	ParameterValue<double> detuningEnergy = ParameterValue(10.0, "detuning energy", "%.3f eV");

	ParameterValue<double> transverse_kT = ParameterValue(2.0e-3, "transverse kT", "%.2e eV");
	ParameterValue<double> coolingEnergy = ParameterValue(0.15, "cooling energy", "%.3f eV");
	ParameterValue<double> cathodeRadius = ParameterValue(0.0012955, "cathode radius", "%.3e m");
	ParameterValue<double> expansionFactor = ParameterValue(30.0, "expansion factor", "%.1f");

	ParameterValue<double> electronCurrent = ParameterValue(1.2e-08, "electron current", "%.2e A");
	ParameterValue<double> cathodeTemperature = ParameterValue(300.0, "cathode tempperature", "%.1f K");
	//ParameterValue<double> Ecath = ParameterValue(0.0, "cathode energy", "%.3e eV");
	ParameterValue<double> LLR = ParameterValue(1.9, "LLR", "%.1f");
	// sigma of gaussian spread of the acceleration voltage (Elab) [eV]
	ParameterValue<double> sigmaLabEnergy = ParameterValue(0.0, "sigma lab energy", "%.3f eV");
	ParameterValue<double> extractionEnergy = ParameterValue(31.26, "extraction energy", "%.3f eV");
	
	ParameterValue<Path> densityFile = ParameterValue(Path(""), "density file", "%s");

	// parameters for analytical electron beam shape to test model
	ParameterValue<bool> hasGaussianShape = ParameterValue(false, "using gaussian beam shape", "%d");
	ParameterValue<bool> hasNoBending = ParameterValue(false, "exclude bend", "%d");
	ParameterValue<bool> hasFixedLongitudinalTemperature = ParameterValue(false, "using fixed longitudial Temperature", "%d");
	ParameterValue<double> radius = ParameterValue(0.003, "radius", "%.4f m");
	ParameterValue<double> longitudinal_kT = ParameterValue(0.0, "longitudinal kT", "%.2e eV");

private:
	int GetSize() override
	{
		return sizeof(*this);
	}
};

class ElectronBeam : public Distribution3D
{
public:
	ElectronBeam();
	void SetupDistribution(std::filesystem::path densityfile = "") override;
	ElectronBeamParameters& GetParameter();
	TH3D* GetDistribution() override;
	void SetParameter(ElectronBeamParameters params);
	void SetCurrent(double current);
	void SetLong_kTFromCenterLabEnergy(double centerLabEnergy);

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
	ElectronBeamParameters m_parameters;

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



