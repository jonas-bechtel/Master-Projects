#pragma once

#include "MCMC.h"
#include "EnergyDistributionModel.h"
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
	EnergyDistributionModel model;

	CrossSection crossSection;
};

