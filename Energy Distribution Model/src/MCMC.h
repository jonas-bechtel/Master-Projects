#pragma once

#include "imgui.h"

#include <TCanvas.h> 
#include <TH1F.h> 
#include <TH3D.h>
#include <TGraph.h>

#include <string>
#include <iostream>
#include <vector>
#include <future>
#include <random>
#include <math.h>

#include "Module.h"
#include "Point3D.h"

class MCMC : public Module
{
public:
	MCMC();
	void SetTargetDistribution(TH3D* targetDist);
	//void GenerateTargetDistribution();
	std::vector<Point3D>& GetSamples();

private:
	void ShowUI();

	void GenerateChain();
	float GenerateSubchain(int index);
	bool GenerateSample(Point3D& current, double& currentValue, std::mersenne_twister_engine<std::uint_fast64_t,
		64, 312, 156, 31,
		0xb5026f5aa96619e9, 29,
		0x5555555555555555, 17,
		0x71d67fffeda60000, 37,
		0xfff7eee000000000, 43,
		6364136223846793005>& generator);

	void PlotTargetDistribution();
	void PlotAutocorrelation();
	void PlotProjections();

private:
	
	// graphs
	TH3D* targetDist = nullptr;
	TH3D* targetDistSmall = nullptr;
	TH1D* projectionXTarget = nullptr;
	TH1D* projectionYTarget = nullptr;
	TH1D* projectionZTarget = nullptr;
	TH1D* projectionX = nullptr;
	TH1D* projectionY = nullptr;
	TH1D* projectionZ = nullptr;

	// axis ranges of the target distribution adn the generated one
	double axisRanges[6] = { 0, 0, 0, 0, 0, 0 };

	int chainLength = 3e5;
	int burnIn = 1000;
	int lag = 30;
	float proposalSigma[3] = {0.005, 0.005, 0.2}; // { 1.5, 0.8, 0.8 };
	float acceptanceRate = 0.0;
	int seed = std::time(0);
	bool changeSeed = true;

	bool generateAsync = true;
	int numThreads = std::thread::hardware_concurrency();
	std::vector<std::future<float>> futures;

	std::vector<Point3D> chain;
	std::vector<Point3D> subchains[4];

	//std::default_random_engine generator;
	std::mersenne_twister_engine<std::uint_fast64_t,
		64, 312, 156, 31,
		0xb5026f5aa96619e9, 29,
		0x5555555555555555, 17,
		0x71d67fffeda60000, 37,
		0xfff7eee000000000, 43,
		6364136223846793005> generator; // best one so far
	//std::subtract_with_carry_engine< std::uint_fast64_t, 48, 5, 12> generator; // is slightly faster but maybe less quality
	//std::linear_congruential_engine<std::uint_fast32_t, 48271, 0, 2147483647> generator; // is even faster
	std::uniform_real_distribution<double> uniformDist = std::uniform_real_distribution<double>(0.0, 1.0);
	std::normal_distribution<double> normalDistX = std::normal_distribution<double>(0.0, proposalSigma[0]);
	std::normal_distribution<double> normalDistY = std::normal_distribution<double>(0.0, proposalSigma[1]);
	std::normal_distribution<double> normalDistZ = std::normal_distribution<double>(0.0, proposalSigma[2]);

	double totalTime = 0.0;
	double interpolationTime = 0;
};

