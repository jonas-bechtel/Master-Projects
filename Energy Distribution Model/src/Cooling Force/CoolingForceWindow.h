#pragma once
namespace CoolingForceWindow
{
	void Init();

	void CreateNewCurve();
	void SetupCurve(std::filesystem::path folder, std::filesystem::path subfolder = "");

	void ShowWindow();
	void ShowSettings();
	void ShowTabs();
	void ShowPlots();
	void ShowAllParametersWindow();
}

