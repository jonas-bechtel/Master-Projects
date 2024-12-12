#pragma once

#include "ParameterImplementations.h"

struct EnergyDistributionSet;
struct EnergyDistribution;
class MCMC;
class ElectronBeam;
class IonBeam;
class LabEnergies;
class EnergyDistributionManager;

struct CrossSection;
struct RateCoefficient;
struct PlasmaRateCoefficient;
class RateCoefficientManager;
class CrossSectionManager;

class Window
{
public:
	Window(std::string name);
	virtual ~Window();

	void ShowWindow();

private:
	bool IsCanvasShown(TCanvas* canvas);
	void ShowCanvas(TCanvas* canvas);
	void HideCanvas(TCanvas* canvas);
	void ShowHideCanvasButton(TCanvas* canvas);
	void UpdateCanvas();

	virtual void ShowUI() = 0;

protected:	
	std::string m_name;
	TCanvas* m_mainCanvas;
	TCanvas* m_secondCanvas;
};

class EnergyDistribtionSetsContainer
{
protected:
	static std::vector<EnergyDistributionSet> energyDistributionSets;
};

class EnergyDistributionModule : public Window
{
public:
	EnergyDistributionModule(std::string name);
	virtual void SetupDistribution(std::filesystem::path = "") = 0;
	virtual TH3D* GetDistribution();
	void PlotDistribution();
	 
	virtual ~EnergyDistributionModule();

protected:
	bool RebinningFactorInput();

protected:
	TH3D* m_distribution;
	TH3D* m_distributionSmall;

	// shared working energy distribution
	static EnergyDistribution activeDist;

	// all the components affecting it
	static MCMC* mcmc;
	static ElectronBeam* eBeam;
	static IonBeam* ionBeam;
	static LabEnergies* labEnergies;
	static EnergyDistributionManager* manager;

	static int s_rebinningFactors[3];
};

class CrossSectionDeconvolutionModule : public EnergyDistribtionSetsContainer, public Window
{
public:
	CrossSectionDeconvolutionModule(std::string name);

	virtual ~CrossSectionDeconvolutionModule();

protected:

protected:
	// shared lists of data	
	static std::vector<CrossSection> crossSectionList;
	static std::vector<RateCoefficient> rateCoefficientList;
	static std::vector<PlasmaRateCoefficient> plasmaRateCoefficientList;

	// components working on the data
	static CrossSectionManager* CSmanager;
	static RateCoefficientManager* RCmanager;
};

