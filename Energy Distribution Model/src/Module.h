#pragma once

#include <string>
#include <unordered_map>

#include <TCanvas.h>
#include <TH3D.h>

class Module
{
public:
	Module(std::string name);

	void ShowWindow();
	TH3D* GetDistribution();

	virtual ~Module();

protected:
	void PlotDistribution();
	void RebinningFactorInput();

private:
	bool IsCanvasShown(TCanvas* canvas);
	void ShowCanvas(TCanvas* canvas);
	void HideCanvas(TCanvas* canvas);
	void ShowHideCanvasButton(TCanvas* canvas);
	void UpdateCanvas();

	virtual void ShowUI() = 0;

protected:
	static std::unordered_map<std::string, Module*> s_moduleMap;
	static int s_rebinningFactors[3];

	std::string m_name;
	TCanvas* m_mainCanvas;
	TCanvas* m_secondCanvas;
	TH3D* m_distribution;
	TH3D* m_distributionSmall;

};

