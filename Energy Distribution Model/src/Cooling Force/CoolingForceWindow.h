#pragma once
namespace CoolingForce
{
	void Init();

	void CreateNewCurve();
	void SetupCurve(std::filesystem::path folder, std::filesystem::path subfolder = "");

	float GetSliceValue();

	void ShowWindow();
	void ShowSettings();
	void ShowTabs();
	void ShowPlots();
	void ShowAllParametersWindow();
	void ShowForceDetailWindow();
}

