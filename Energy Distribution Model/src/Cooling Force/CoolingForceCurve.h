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
		bool IsNumerical() const;
		bool IsMeasured() const;

		void ShowContent();
		void SelectedItemChanged();
		void UpdateMovedScaledValues();

		void PlotForceX() const;
		void PlotForceY() const;
		void PlotForceZ() const;
		void PlotDetails() const;

		void UpdateSlice(float zValue);

		void Save() const;
		void Load(const std::filesystem::path& input);
		void LoadMeasured(const std::filesystem::path& input);

	private:
		std::vector<Value> values;
		std::vector<double> forceX;
		std::vector<double> forceY;
		std::vector<double> forceZ;
		std::vector<double> forceZError;
		std::vector<double> detuningVelocites;

		// moved and scaled versions of the force, velocity lists
		std::vector<double> forceZmovedScaled;
		std::vector<double> forceZErrorScaled;
		std::vector<double> detuningVelocitesMoved;

		float scale = 1.0f;
		float moveForce = 0.0f;
		float moveVelocity = 0.0f;

		std::filesystem::path folder = "Test";
		std::filesystem::path subFolder = "subfolder";

		int selectedIndex = -1;
		bool numerical = false; // as opposed to mcmc generated
		bool measured = false;
		Model::Parameter modelParams;
		MeasuredCurveParameter measuredParams;
	};
}


