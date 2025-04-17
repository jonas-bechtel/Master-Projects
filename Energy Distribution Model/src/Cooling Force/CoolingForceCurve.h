#pragma once
#include "CoolingForceValue.h"
#include "CoolingForceModel.h"

namespace CoolingForce
{
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

		void ShowContent();
		void SelectedItemChanged();

		void PlotForceX() const;
		void PlotForceY() const;
		void PlotForceZ() const;
		void PlotDetails() const;

		void UpdateSlice(float zValue);

		void Save() const;
		void Load(const std::filesystem::path& input);

	private:
		std::vector<Value> values;
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
		Model::Parameter modelParams;
	};
}


