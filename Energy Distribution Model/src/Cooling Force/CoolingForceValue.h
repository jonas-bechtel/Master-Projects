#pragma once

#include "HeatMapData.h"
#include "HistData3D.h"
#include "ElectronBeam.h"
#include "IonBeam.h"
#include "LabEnergies.h"

#include "CoolingForceModel.h"

using RNG_engine = std::mersenne_twister_engine<std::uint_fast64_t,
	64, 312, 156, 31,
	0xb5026f5aa96619e9, 29,
	0x5555555555555555, 17,
	0x71d67fffeda60000, 37,
	0xfff7eee000000000, 43,
	6364136223846793005>;

namespace CoolingForce
{
	struct Value
	{
	public:
		Value();
		~Value();

		Value(const Value& other) = delete;
		Value& operator=(const Value& other) = delete;
		Value(Value&& other) noexcept;
		Value& operator=(Value&& other) noexcept;

		void Calculate(std::filesystem::path descriptionFile, int index, Model::Parameter params);
		double Integrand(double* position, double* params);
		//double IntegrandTrans(double* position, double* params);

		bool ShowListItem(bool selected) const;
		static bool ShowParallelPrecalculationCheckbox();
		static bool ShowCalcTransForceCheckbox();

		void SetupHistogramsFromReference(TH3D* reference);

		void PlotForceSlice() const;
		void UpdateSlice(float zValue);

		ElectronBeamParameters GetElectronBeamParameters() const;
		IonBeamParameters GetIonBeamParameters() const;
		LabEnergyParameters GetLabEnergyParameters() const;

		void Save(std::filesystem::path folder) const;
		void Load(std::filesystem::path file);

	private:
		void PrepareCalculation(std::filesystem::path descriptionFile, int index);
		void CalculateForce3D(Model::Parameter params);
		void CopyParameters();
		void SetupLabel();
		void SetupTags();
		void ResetDefaultValues();

		std::string GetHeaderString() const;
		std::string Filename() const;

	private:
		HistData3D force3D;
		//HistData3D forceTrans;

		// actual resulting value of the force to be compared to measurements (integral divided by length)
		double forceValue = 0.0;
		//double forceTransValue = 0.0;

		// all the parameters used to create it
		ElectronBeamParameters eBeamParameter;
		IonBeamParameters ionBeamParameter;
		LabEnergyParameters labEnergiesParameter;

		// additional labelling things
		std::string label = "";
		std::string tags = "";
		int index = 0;

		static inline bool parallelForcePrecalculation = true;
		static inline bool calculateTransverseForce = false;

		// random number generation things
		//static std::mersenne_twister_engine<std::uint_fast64_t,
		//	64, 312, 156, 31,
		//	0xb5026f5aa96619e9, 29,
		//	0x5555555555555555, 17,
		//	0x71d67fffeda60000, 37,
		//	0xfff7eee000000000, 43,
		//	6364136223846793005> generator;
		//
		//static std::normal_distribution<double> longitudinalNormalDistribution;
		//static std::normal_distribution<double> transverseNormalDistribution;

		friend class Curve;
	};
}


