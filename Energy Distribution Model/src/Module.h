#pragma once

#include "ParameterImplementations.h"

struct EnergyDistribution;

class Module
{
public:
	Module(std::string name);
	virtual ~Module();

	static Module* Get(std::string name);
	static std::unordered_map<std::string, Module*>& GetModuleMap();

	void ShowWindow();

private:
	bool IsCanvasShown(TCanvas* canvas);
	void ShowCanvas(TCanvas* canvas);
	void HideCanvas(TCanvas* canvas);
	void ShowHideCanvasButton(TCanvas* canvas);
	void UpdateCanvas();

	virtual void ShowUI() = 0;

protected:
	static std::unordered_map<std::string, Module*> s_moduleMap;
	
	std::string m_name;
	TCanvas* m_mainCanvas;
	TCanvas* m_secondCanvas;
};

class EnergyDistributionModule : public Module
{
public:
	EnergyDistributionModule(std::string name);
	virtual void SetupDistribution(std::filesystem::path = "") = 0;
	virtual TH3D* GetDistribution();
	void PlotDistribution();
	 
	virtual ~EnergyDistributionModule() override;

protected:
	bool RebinningFactorInput();

protected:
	static EnergyDistribution activeDist;

	TH3D* m_distribution;
	TH3D* m_distributionSmall;

	static int s_rebinningFactors[3];

};

