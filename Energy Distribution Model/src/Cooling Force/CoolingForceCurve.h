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

	void ShowList();

	void PlotForceZ() const;

	//void Save() const;
	//void Load(std::filesystem::path& folder);

private:
	std::vector<CoolingForceValue> values;
	std::vector<double> forceX;
	std::vector<double> forceY;
	std::vector<double> forceZ;
	std::vector<double> detuningVelocites;

	std::filesystem::path folder = "Test";
	std::filesystem::path subFolder = "subfolder";

	int selectedIndex = -1;
	bool numerical = false; // as opposed to mcmc generated
};

