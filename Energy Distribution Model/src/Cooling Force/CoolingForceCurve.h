#pragma once
#include "CoolingForceValue.h"
#include "CoolingForceModel.h"

class CoolingForceCurve
{
public:
	CoolingForceCurve();
	CoolingForceCurve(const CoolingForceCurve& other) = delete;
	CoolingForceCurve& operator=(const CoolingForceCurve& other) = delete;

	CoolingForceCurve(CoolingForceCurve&& other) = default;
	CoolingForceCurve& operator=(CoolingForceCurve&& other) = default;

	void IntegrateNumerically(NumericalIntegrationParameter& params);
	void AddForceValue(CoolingForceValue&& value);
	void RemoveForceValue(int index);

	void SetFolder(std::filesystem::path path);
	void SetSubfolder(std::filesystem::path path);

	std::filesystem::path GetFolder() const;
	std::filesystem::path GetSubfolder() const;
	std::string GetLabel() const;
	bool Empty() const;
	bool IsNumerical() const;

	void ShowList();
	void SelectedItemChanged();

	void PlotForceX() const;
	void PlotForceY() const;
	void PlotForceZ() const;
	void PlotDetails() const;

	void UpdateSlice(float zValue);

	void Save() const;
	void Load(const std::filesystem::path& input);

private:
	std::vector<CoolingForceValue> values;
	std::vector<double> forceX;
	std::vector<double> forceY;
	std::vector<double> forceZ;
	std::vector<double> forceZscaled;
	std::vector<double> detuningVelocites;

	float scale = 1.0f;

	std::filesystem::path folder = "Test";
	std::filesystem::path subFolder = "subfolder";

	int selectedIndex = -1;
	bool numerical = false; // as opposed to mcmc generated
	NumericalIntegrationParameter numericalParams;
};

