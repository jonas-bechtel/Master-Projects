#pragma once

#include "MCMC.h"
#include "EnergyDistributionManager.h"
#include "CrossSection.h"

#include <TApplication.h>  // For initializing the ROOT application
#include <TSystem.h>


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
	EnergyDistributionManager model;

	CrossSection crossSection;
};

