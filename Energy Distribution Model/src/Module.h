#pragma once
#include "imgui.h"
#include "implot.h"

#include <string>
#include <unordered_map>

#include <TCanvas.h>
#include <TH3D.h>

class Module
{
public:
	Module(std::string name);
	static Module* Get(std::string name);
	static std::unordered_map<std::string, Module*>& GetModuleMap();

	void ShowWindow();
	virtual TH3D* GetDistribution();
	void PlotDistribution();

	virtual ~Module();

protected:
	virtual void PlotAllOnMainCanvas();
	virtual void PlotAllOnSecondCanvas();
	bool RebinningFactorInput();

private:
	bool IsCanvasShown(TCanvas* canvas);
	void ShowCanvas(TCanvas* canvas);
	void HideCanvas(TCanvas* canvas);
	void ShowHideCanvasButton(TCanvas* canvas);
	void UpdateCanvas();

	virtual void ShowUI() = 0;
	virtual void ShowPlots() {};

protected:
	static std::unordered_map<std::string, Module*> s_moduleMap;
	static int s_rebinningFactors[3];

	std::string m_name;
	TCanvas* m_mainCanvas;
	TCanvas* m_secondCanvas;

	TH3D* m_distribution;
	TH3D* m_distributionSmall;

	bool m_hasPlotWindow = false;
};

