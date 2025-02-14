#pragma once

#include "Module.h"
#include "Point3D.h"
#include "HeatMapData.h"

using RNG_engine = std::mersenne_twister_engine<std::uint_fast64_t,
	64, 312, 156, 31,
	0xb5026f5aa96619e9, 29,
	0x5555555555555555, 17,
	0x71d67fffeda60000, 37,
	0xfff7eee000000000, 43,
	6364136223846793005>;

struct MCMC_Data
{
	TH3D* targetDist;
	TH3D* generatedDist;

	std::vector<double> xAxis;
	std::vector<double> yAxis;
	std::vector<double> zAxis;

	std::vector<double> targetProjectionValuesX;
	std::vector<double> targetProjectionValuesY;
	std::vector<double> targetProjectionValuesZ;

	std::vector<double> generatedProjectionValuesX;
	std::vector<double> generatedProjectionValuesY;
	std::vector<double> generatedProjectionValuesZ;

	HeatMapData targetSlice;
	HeatMapData generatedSlice;

	std::string label;

	void FillData();

	MCMC_Data() {}
	MCMC_Data(const MCMC_Data& other) = delete;
	MCMC_Data& operator=(const MCMC_Data& other) = delete;
	MCMC_Data(MCMC_Data&& other);
	MCMC_Data& operator=(MCMC_Data&& other);
	~MCMC_Data()
	{
		delete targetDist;
		delete generatedDist;
	}
};

class MCMC_Window : public EnergyDistributionModule
{
public:
	MCMC_Window();
	void SetupDistribution(std::filesystem::path file = "") override;
	std::vector<Point3D>& GetSamples();
	void SetParameter(MCMC_Parameters params);
	void GenerateSamples();

private:
	void ShowUI() override;
	void ShowList();
	void ShowSettings();
	void ShowMCMCDataPlots();

	float GenerateSubchain(int length, int offset, RNG_engine& generator);
	bool GenerateSingleSample(Point3D& current, double& currentValue, RNG_engine& generator);

	void AddMCMCDataToList(MCMC_Data& mcmcData);
	void RemoveMCMCDataFromList(int index);

	void PlotTargetDistribution();
	void PlotAutocorrelation();
	//void PlotProjections();

private:
	MCMC_Parameters& m_parameters;

	std::vector<Point3D> chain;

	TH3D* m_distributionSmall = nullptr;
	TH3D* targetDist = nullptr;
	TH3D* targetDistSmall = nullptr;

	std::vector<MCMC_Data> mcmcDataToLookAt;
	int selectedIndex = -1;
	float SliceZ = 0.0f;

	// graphs
	//TH3D* targetDist = nullptr;
	//TH3D* generatedDist = nullptr;
	//TH1D* projectionXTarget = nullptr;
	//TH1D* projectionYTarget = nullptr;
	//TH1D* projectionZTarget = nullptr;
	//TH1D* projectionX = nullptr;
	//TH1D* projectionY = nullptr;
	//TH1D* projectionZ = nullptr;

	float acceptanceRate = 0.0;
	bool changeSeed = true;
	bool automaticProposalStd = true;
	bool useInterpolation = false;

	bool generateAsync = true;
	int numThreads = std::thread::hardware_concurrency();
	std::vector<std::future<float>> futures;

	//std::default_random_engine generator;
	std::vector<RNG_engine> generatorList;
	//std::subtract_with_carry_engine< std::uint_fast64_t, 48, 5, 12> generator; // is slightly faster but maybe less quality
	//std::linear_congruential_engine<std::uint_fast32_t, 48271, 0, 2147483647> generator; // is even faster
	std::uniform_real_distribution<double> uniformDist = std::uniform_real_distribution<double>(0.0, 1.0);
	std::normal_distribution<double> normalDistX = std::normal_distribution<double>(0.0, m_parameters.proposalSigma.get().x);
	std::normal_distribution<double> normalDistY = std::normal_distribution<double>(0.0, m_parameters.proposalSigma.get().y);
	std::normal_distribution<double> normalDistZ = std::normal_distribution<double>(0.0, m_parameters.proposalSigma.get().z);

	double totalTime = 0.0;
	double interpolationTime = 0;
};

