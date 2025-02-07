#pragma once

#include "ParameterImplementations.h"

struct EnergyDistributionSet;
struct EnergyDistribution;
class MCMC;
class ElectronBeamWindow;
class IonBeam;
class LabEnergyWindow;
class EnergyDistributionManager;

struct CrossSection;
struct RateCoefficient;
struct PlasmaRateCoefficient;
class DeconvolutionManager;

class Window
{
public:
	Window(std::string name, int numberCanvasses = 2);
	virtual ~Window();

	void ShowWindow();

protected:
	void ShowCanvasButtons(); 

private:
	bool IsCanvasShown(TCanvas* canvas);
	void ShowCanvas(TCanvas* canvas);
	void HideCanvas(TCanvas* canvas);
	void ShowHideCanvasButton(TCanvas* canvas);
	void UpdateCanvas();

	virtual void ShowUI() = 0;

protected:	
	std::string m_name;
	TCanvas* m_mainCanvas = nullptr;
	TCanvas* m_secondCanvas = nullptr;
};

class EnergyDistribtionSetsContainer
{
protected:
	static std::vector<EnergyDistributionSet> energyDistributionSets;
	static int currentSetIndex;
};

class EnergyDistributionModule : public Window
{
public:
	EnergyDistributionModule(std::string name, int numberCanvasses = 2);
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
	static ElectronBeamWindow* eBeam;
	static IonBeam* ionBeam;
	static LabEnergyWindow* labEnergies;
	static EnergyDistributionManager* manager;

	static int s_rebinningFactors[3];
};

class CrossSectionDeconvolutionModule : public EnergyDistribtionSetsContainer, public Window
{
public:
	CrossSectionDeconvolutionModule(std::string name, int numberCanvasses = 2);

	virtual ~CrossSectionDeconvolutionModule();

protected:

protected:
	// shared lists of data	
	static std::vector<CrossSection> crossSectionList;
	static std::vector<RateCoefficient> rateCoefficientList;
	static std::vector<PlasmaRateCoefficient> plasmaRateCoefficientList;
	static int currentCrossSectionIndex;
	static int currentRateCoefficientIndex;
	static int currentPlasmaRateCoefficientIndex;

	// components working on the data
	//static CrossSectionManager* CSmanager;
	static DeconvolutionManager* deconvolutionManager;
};

