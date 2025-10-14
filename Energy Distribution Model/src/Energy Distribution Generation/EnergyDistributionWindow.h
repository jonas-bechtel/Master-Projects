#pragma once

#include "EnergyDistributionSet.h"

namespace EnergyDistributionWindow
{
	void Init();

	EnergyDistributionSet& GetCurrentSet();
	std::vector<EnergyDistributionSet>& GetSetList();
	int GetCurrentSetIndex();

	void CreateNewSet();
	void SetupSet(std::filesystem::path folder, std::filesystem::path subfolder = "");
	void LoadSet();

	void ShowWindow();
	void ShowSettings();
	void ShowTabs();
	void ShowPlot();
	void ShowAllParametersWindow();

	void ShowSetInformationWindow();
	void ShowSetList();
}



