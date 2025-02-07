#pragma once

#include "Module.h"
#include "Point3D.h"
#include "HeatMapData.h"

struct ElectronBeam
{
	TH3D* fullHistogram;
	std::vector<double> xAxis;
	std::vector<double> yAxis;
	std::vector<double> zAxis;

	std::vector<double> projectionValuesX;
	std::vector<double> projectionValuesY;
	std::vector<double> projectionValuesZ;

	HeatMapData slice;

	std::string label;

	void FillData();

	ElectronBeam() {}
	ElectronBeam(const ElectronBeam& other) = delete;
	ElectronBeam& operator=(const ElectronBeam& other) = delete;
	ElectronBeam(ElectronBeam&& other);
	ElectronBeam& operator=(ElectronBeam&& other);
	~ElectronBeam()
	{
		//std::cout << "deleting " << fullHistogram << std::endl;
		delete fullHistogram;
	}
};

class ElectronBeamWindow : public EnergyDistributionModule
{
public:
	ElectronBeamWindow();
	void SetupDistribution(std::filesystem::path densityfile = "") override;
	//TH3D* GetDistribution() override;

	void CalculateDetuningEnergy();

	double Trajectory(double z) const;

	TVector3 GetDirection(double z);
	TVector3 GetDirection(Point3D point);
	double GetLongitudinal_kT(double labEnergy);
	double GetTransverse_kT();

private:
	void ShowUI() override;
	void ShowList();
	void ShowSettings();
	void ShowElectronBeamPlots();

	void LoadDensityFile(std::filesystem::path file);
	void LoadToLookAt(std::filesystem::path file);

	void GenerateElectronBeamDensity();

	TH3D* CutZerosFromDistribution(TH3D* input);
	void CreateLargeDistribution();

	TVector3 GetNormal(double z);
	
	double Derivative(double z) const;
	double DistancePointToTrajectoryOfZ(double z, Point3D point);

	//void PlotDensitySlice();
	//void PlotGeneratedDensities();
	void PlotTrajectory();
	//void PlotProjections();

private:
	ElectronBeamParameters& m_parameters;

	std::vector<ElectronBeam> eBeamsToLookAt;
	int selectedIndex = -1;

	TProfile2D* bla = nullptr;

	//TH3D* generatedBeamDensity = nullptr;
	//TH3D* generatedBeamDensitySmall = nullptr;
	//TH2D* densitySliceXY = nullptr;

	//TH1D* electronBeamProjectionX = nullptr;
	//TH1D* electronBeamProjectionY = nullptr;
	//TH1D* electronBeamProjectionZ = nullptr;

	//TH1D* generatedBeamProjectionX = nullptr;
	//TH1D* generatedBeamProjectionY = nullptr;
	//TH1D* generatedBeamProjectionZ = nullptr;

	// parameters to increase histogram resolution by interpolation
	bool increaseHist = false;
	int factor = 3;

	float sliderZ = 0;
	float sliderY = 0;

	float SliceZ = 0.0f;
};



