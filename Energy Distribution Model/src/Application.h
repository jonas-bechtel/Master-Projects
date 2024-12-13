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

	ElectronBeam electronBeam;
	IonBeam ionBeam;
	MCMC mcmc;
	LabEnergies labEnergies;
	EnergyDistributionManager energyDistributionManager;

	DeconvolutionManager deconvolutionManager;

};

