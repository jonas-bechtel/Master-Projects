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

	void CalculateDetuningEnergy();

	double Trajectory(double z) const;

	TVector3 GetDirection(double z);
	TVector3 GetDirection(Point3D point);
	double GetLongitudinal_kT(double labEnergy);
	double GetTransverse_kT();

	ElectronBeam* GetSelected();

private:
	void ShowUI() override;
	void ShowList();
	void ShowSettings();
	void ShowElectronBeamPlots();

	void LoadDensityFile(std::filesystem::path file);
	void LoadToLookAt(std::filesystem::path file);

	void ChangeSelectedItem(int i);
	void AddBeamToList(ElectronBeam& eBeam);
	void RemoveBeamFromList(int index);

	TH3D* GenerateElectronBeamDensity();

	TH3D* CutZerosFromDistribution(TH3D* input);
	void CreateLargeDistribution();

	TVector3 GetNormal(double z);
	
	double Derivative(double z) const;
	double DistancePointToTrajectoryOfZ(double z, Point3D point);

	void PlotTrajectory();
	void Plot3DBeam(TH3D* eBeam);

private:
	ElectronBeamParameters& m_parameters;

	TH3D* m_distributionSmall = nullptr;

	std::vector<ElectronBeam> eBeamsToLookAt;
	int selectedIndex = -1;

	// parameters to increase histogram resolution by interpolation
	bool increaseHist = false;
	int factor = 3;

	float sliderZ = 0;
	float sliderY = 0;

	float SliceZ = 0.0f;
};



