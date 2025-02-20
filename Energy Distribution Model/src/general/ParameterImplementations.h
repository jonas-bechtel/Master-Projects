#pragma once
#include "Parameter.h"

struct MCMC_Parameters : public Parameters
{
	MCMC_Parameters() { setName("mcmc sampling parameters"); }

	ParameterValue<int> numberSamples = ParameterValue((int)3e5, "number of samples", "%d");
	ParameterValue<int> burnIn = ParameterValue(1000, "burn in", "%d");
	ParameterValue<int> lag = ParameterValue(30, "lag", "%d");
	ParameterValue<float3> proposalSigma = ParameterValue(float3(0.005f, 0.005f, 0.2f), "proposal sigmas", "%.4f, %.4f, %.3f m");
	ParameterValue<int> seed = ParameterValue((int)std::time(0), "seed", "%d");

private:
	int GetSize() const override
	{
		return sizeof(*this);
	}
};

struct ElectronBeamParameters : public Parameters
{
	ElectronBeamParameters()
	{
		setName("electron beam parameters");
	}

	ParameterValue<double> detuningEnergy = ParameterValue(10.0, "detuning energy", "%.6e eV");

	ParameterValue<double> transverse_kT = ParameterValue(2.0e-3, "transverse kT", "%.2e eV");
	ParameterValue<double> longitudinal_kT_estimate = ParameterValue(0.0, "estimated longitudinal kT", "%.2e eV");
	ParameterValue<double> coolingEnergy = ParameterValue(0.15263, "cooling energy", "%.6e eV");
	ParameterValue<double> cathodeRadius = ParameterValue(0.0012955, "cathode radius", "%.3e m");
	ParameterValue<double> expansionFactor = ParameterValue(30.0, "expansion factor", "%.1f");

	ParameterValue<double> electronCurrent = ParameterValue(1.2e-08, "electron current", "%.2e A");
	ParameterValue<double> cathodeTemperature = ParameterValue(300.0, "cathode tempperature", "%.1f K");
	//ParameterValue<double> Ecath = ParameterValue(0.0, "cathode energy", "%.3e eV");
	ParameterValue<double> LLR = ParameterValue(1.9, "LLR", "%.1f");
	// sigma of gaussian spread of the acceleration voltage (Elab) [eV]
	ParameterValue<double> sigmaLabEnergy = ParameterValue(0.0, "sigma lab energy", "%.3f eV");
	ParameterValue<double> extractionEnergy = ParameterValue(31.26, "extraction energy", "%.3f eV");

	ParameterValue<Path> densityFile = ParameterValue(Path(""), "density file", "%s");

	
private:
	int GetSize() const override
	{
		return sizeof(*this);
	}
};

struct LabEnergyParameters : public Parameters
{
	LabEnergyParameters() { setName("lab energy parameters"); }

	ParameterValue<double> centerLabEnergy = ParameterValue(0.0, "lab energy in center", "%e eV");
	ParameterValue<double> driftTubeVoltage = ParameterValue(0.0, "drift tube voltage", "%e V");
	ParameterValue<Path> energyFile = ParameterValue(Path(""), "energy file", "%s");

private:
	int GetSize() const override
	{
		return sizeof(*this);
	}
};

struct IonBeamParameters : public Parameters
{
	IonBeamParameters()
	{
		setName("ion beam parameters");
	}

	ParameterValue<float2> shift = ParameterValue(float2(0.0f, 0.0f), "shift in x and y", "%.4f, %.4f m");
	ParameterValue<float2> angles = ParameterValue(float2(0.0f, 0.0f), "horizontal, vertical angle", "%.4f, %.4f rad");

	// always one gaussian
	ParameterValue<double> amplitude = ParameterValue(10.1, "amplitude", "%.4f");
	ParameterValue<float2> sigma = ParameterValue(float2(9.5e-3f, 5.7e-3f), "sigmas (x,y)", "%.4f, %.4f m");

	// maybe a second one
	ParameterValue<double> amplitude2 = ParameterValue(8.1, "amplitude 2", "%.4f");
	ParameterValue<float2> sigma2 = ParameterValue(float2(1.39e-3f, 2.15e-3f), "sigmas 2 (x,y)", "%.4f, %.4f m");

private:
	int GetSize() const override
	{
		return sizeof(*this);
	}
};

struct OutputParameters : public Parameters
{
	OutputParameters()
	{
		setName("analytical/fit distribution parameters");
	}

	ParameterValue<float2> fitRange = ParameterValue(float2(0.0f, 1.0f), "fit range", "%.4e, %.4e eV");
	ParameterValue<double> fitDetuningEnergy = ParameterValue(1.0, "fit detuning energy", "%.6e eV");
	ParameterValue<double> fitLongitudinalTemperature = ParameterValue(0.0005, "fit longitudinal kT", "%.2e eV");
	ParameterValue<double> fitTransverseTemperature = ParameterValue(0.002, "fit transverse kT", "%.2e eV");
	ParameterValue<double> fitFWHM = ParameterValue(0.0, "fit FWHM", "%.4f eV");
	ParameterValue<double> fitScalingFactor = ParameterValue(0.0, "fit scaling factor", "%.3f");
	ParameterValue<double> effectiveLength = ParameterValue(0.0, "effective length", "%.3f m");

	ParameterValue<double> FWHM = ParameterValue(0.0, "FWHM", "%.4f eV");

private:
	int GetSize() const override
	{
		return sizeof(*this);
	}
};

struct SimplificationParameter : public Parameters
{
	SimplificationParameter()
	{
		setName("simplification parameter");
	}

	// ion beam
	//ParameterValue<bool> singleGaussianIonBeam = ParameterValue(false, "using single gaussian ion beam", "%d");
	//ParameterValue<double> ionBeamRadius = ParameterValue(0.0010, "ion beam radius", "%.4f m");

	// lab energies
	ParameterValue<bool> uniformLabEnergies = ParameterValue(false, "use uniform lab energy", "%d");
	ParameterValue<bool> sliceLabEnergies = ParameterValue(false, "use slice of lab energies", "%d");
	ParameterValue<double> sliceToFill = ParameterValue(0.5, "z value of slice", "%.3f");
	
	// electron beam
	ParameterValue<bool> gaussianElectronBeam = ParameterValue(false, "using gaussian e-beam shape", "%d");
	ParameterValue<bool> cylindricalElectronBeam = ParameterValue(false, "using cylindrical e-beam shape", "%d");
	ParameterValue<bool> noElectronBeamBend = ParameterValue(false, "exclude bend", "%d");
	ParameterValue<bool> fixedLongitudinalTemperature = ParameterValue(false, "using fixed longitudial Temperature", "%d");
	ParameterValue<double> electronBeamRadius = ParameterValue(0.003, "e-beam radius", "%.4f m");

	// during generation
	ParameterValue<bool> cutOutZValues = ParameterValue(false, "cut out z values", "%d");
	ParameterValue<float2> cutOutRange = ParameterValue(float2(0.0f, 0.4f), "cut out range", "%.2f, %.2f m");


private:
	int GetSize() const override
	{
		return sizeof(*this);
	}
};
