#pragma once
#include "CoolingForceValue.h"
#include "CoolingForceModel.h"

namespace CoolingForce
{
	struct MeasuredCurveParameter
	{
		int ionCharge = 1;
		double coolingEnergy = 0;
		double effectiveBunchingVoltageDirect = 0;
		double effectiveBunchingVoltageSync = 0;

		//void FromString(std::string& input);
		//void ShowValues();
	};

	class Curve
	{
	public:
		Curve();
		Curve(const Curve& other) = delete;
		Curve& operator=(const Curve& other) = delete;

		Curve(Curve&& other) = default;
		Curve& operator=(Curve&& other) = default;

		void IntegrateNumerically(const Model::Parameter& params);
		void AddForceValue(Value&& value);
		void RemoveForceValue(int index);

		void SetFolder(std::filesystem::path path);
		void SetSubfolder(std::filesystem::path path);

		std::filesystem::path GetFolder() const;
		std::filesystem::path GetSubfolder() const;
		std::string GetLabel() const;
		bool Empty() const;
		bool IsSimpleModel() const;
		bool IsMeasured() const;

		void ShowContent();
		void SelectedItemChanged();
		void UpdateMovedScaledValues();

		void PlotForce() const;
		void PlotDetails() const;

		void UpdateSlice(float zValue);

		void Save() const;
		void Load(const std::filesystem::path& input);
		void LoadMeasured(const std::filesystem::path& input);

	private:
		std::vector<Value> values;
		std::vector<double> force;
		std::vector<double> forceError;
		std::vector<double> detuningVelocites;

		// moved and scaled versions of the force, velocity lists
		std::vector<double> forceMovedScaled;
		std::vector<double> forceErrorMovedScaled;
		std::vector<double> detuningVelocitesMovedScaled;

		float scale = 1.0f;
		float moveForce = 0.0f;
		float moveVelocity = 0.0f;

		std::filesystem::path folder = "Test";
		std::filesystem::path subFolder = "subfolder";

		int selectedIndex = -1;
		bool simpleModel = false; // as opposed to 3D model
		bool measured = false;
		Model::Parameter modelParams;
		MeasuredCurveParameter measuredParams;
	};
}


