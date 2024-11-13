#pragma once
#include "imgui.h"
#include "implot.h"

#include <string>
#include <unordered_map>
#include <filesystem>

#include <TCanvas.h>
#include <TH3D.h>

#include "Parameter.h"

class Module
{
public:
	Module(std::string name);
	virtual ~Module();

	static Module* Get(std::string name);
	static std::unordered_map<std::string, Module*>& GetModuleMap();

	//Parameters GetParameter();

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
	
	//Parameters m_parameters;

	std::string m_name;
	TCanvas* m_mainCanvas;
	TCanvas* m_secondCanvas;
};

class Distribution3D : public Module
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

