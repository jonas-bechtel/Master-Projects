#pragma once
namespace CoolingForceWindow
{
	void Init();

	void CreateNewCurve();
	void SetupCurve(std::filesystem::path folder, std::filesystem::path subfolder = "");

	float GetSlice();

	void ShowWindow();
	void ShowSettings();
	void ShowTabs();
	void ShowPlots();
	void ShowAllParametersWindow();
	void ShowForceDetailWindow();
}

