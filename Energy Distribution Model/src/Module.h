#pragma once

#include "ParameterImplementations.h"

struct EnergyDistribution;

class EnergyDistributionModule
{
public:
	EnergyDistributionModule(std::string name);
	virtual ~EnergyDistributionModule();

	static EnergyDistributionModule* Get(std::string name);
	static std::unordered_map<std::string, EnergyDistributionModule*>& GetModuleMap();

	void ShowWindow();

private:
	bool IsCanvasShown(TCanvas* canvas);
	void ShowCanvas(TCanvas* canvas);
	void HideCanvas(TCanvas* canvas);
	void ShowHideCanvasButton(TCanvas* canvas);
	void UpdateCanvas();

	virtual void ShowUI() = 0;

protected:
	static std::unordered_map<std::string, EnergyDistributionModule*> s_moduleMap;
	
	static EnergyDistribution activeDist;
	std::string m_name;
	TCanvas* m_mainCanvas;
	TCanvas* m_secondCanvas;
};

class Distribution3D : public EnergyDistributionModule
{
public:
	Distribution3D(std::string name);
	virtual void SetupDistribution(std::filesystem::path = "") = 0;
	virtual TH3D* GetDistribution();
	void PlotDistribution();
	 
	virtual ~Distribution3D() override;

protected:
	bool RebinningFactorInput();

protected:
	TH3D* m_distribution;
	TH3D* m_distributionSmall;

	static int s_rebinningFactors[3];

};

