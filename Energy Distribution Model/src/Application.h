#pragma once

#include "MCMC.h"
#include "EnergyDistributionManager.h"
#include "DeconvolutionManager.h"

class Application
{
public:
	void Run();
	Application();
	
private:
	void ShowWindows();

private:
	TApplication app = TApplication("app", nullptr, nullptr);

	ElectronBeamWindow electronBeam;
	IonBeamWindow ionBeam;
	MCMC_Window mcmc;
	LabEnergyWindow labEnergies;
	EnergyDistributionManager energyDistributionManager;

	DeconvolutionManager deconvolutionManager;

};

