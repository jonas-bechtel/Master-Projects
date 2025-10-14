#pragma once

#include "Point3D.h"
#include "HeatMapData.h"
#include "ROOTCanvas.h"
#include "HistData3D.h"
#include "ParameterImplementations.h"

using RNG_engine = std::mersenne_twister_engine<std::uint_fast64_t,
	64, 312, 156, 31,
	0xb5026f5aa96619e9, 29,
	0x5555555555555555, 17,
	0x71d67fffeda60000, 37,
	0xfff7eee000000000, 43,
	6364136223846793005>;

namespace MCMC
{
	void Init();

	MCMC_Parameters GetParameters();
	void SetParameters(const MCMC_Parameters& params);
	std::vector<Point3D>& GetSamples();

	std::string GetTags();

	void SetTargetDist(TH3D* targeDist);

	void GenerateSamples();

	float GenerateSubchain(int length, int offset, RNG_engine& generator);
	bool GenerateSingleSample(Point3D& current, double& currentValue, RNG_engine& generator);

	void SelectedItemChanged();
	void AddMCMCDataToList(HistData3D& target, HistData3D& result);
	void RemoveMCMCDataFromList(int index);

	void UpdateAutocorrelationData();

	void ShowWindow();
	void ShowList();
	void ShowParameterControls();
	void ShowPlots();
	void ShowAutoCorrelationPlots();
}
