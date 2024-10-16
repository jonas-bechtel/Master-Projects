#include "Module.h"

#include "imgui.h"

#include <TRootCanvas.h>

std::unordered_map<std::string, Module*> Module::s_moduleMap;
int Module::s_rebinningFactors[3] = { 10, 10, 10 };

Module::Module(std::string name)
	: m_name(name)
{
	s_moduleMap[name] = this;
	
	m_mainCanvas = new TCanvas((m_name + " main canvas").c_str(), (m_name + " main canvas").c_str(), 1200, 500);
	m_mainCanvas->Divide(2, 1);

	m_secondCanvas = new TCanvas((m_name + " analysis canvas").c_str(), (m_name + " analysis canvas").c_str(), 1500, 800);
	m_secondCanvas->Divide(3, 2);

	HideCanvas(m_mainCanvas);
	HideCanvas(m_secondCanvas);
}

Module* Module::Get(std::string name)
{
	return s_moduleMap.at(name);
}

std::unordered_map<std::string, Module*>& Module::GetModuleMap()
{
	return s_moduleMap;
}

void Module::ShowWindow()
{
	if (ImGui::Begin((m_name + " Window").c_str()))
	{
		ShowUI();
		ImGui::Separator();
		ShowHideCanvasButton(m_mainCanvas);
		ImGui::SameLine();
		ShowHideCanvasButton(m_secondCanvas);
		
		ImGui::End();
	}

	UpdateCanvas();
}

TH3D* Module::GetDistribution()
{
	return m_distribution;
}

Module::~Module()
{
	delete m_distribution;
	delete m_distributionSmall;
	delete m_mainCanvas;
	delete m_secondCanvas;
}

void Module::PlotDistribution()
{
	if (!m_distribution) return;
	if (m_distributionSmall) delete m_distributionSmall;

	m_mainCanvas->cd(2);
	m_distributionSmall = (TH3D*)m_distribution->Rebin3D(s_rebinningFactors[0],
														 s_rebinningFactors[1], 
														 s_rebinningFactors[2], "generated dist small");
	m_distributionSmall->Draw("BOX2");
}

bool Module::RebinningFactorInput()
{
	return ImGui::InputInt3("Rebinning factors", s_rebinningFactors);
}

bool Module::IsCanvasShown(TCanvas* canvas)
{
	return ((TRootCanvas*)canvas->GetCanvasImp())->IsMapped();
}

void Module::ShowCanvas(TCanvas* canvas)
{
	((TRootCanvas*)canvas->GetCanvasImp())->MapRaised();
}

void Module::HideCanvas(TCanvas* canvas)
{
	((TRootCanvas*)canvas->GetCanvasImp())->UnmapWindow();
}

void Module::ShowHideCanvasButton(TCanvas* canvas)
{
	if (ImGui::Button(("Show/Hide " + std::string(canvas->GetTitle())).c_str()))
	{
		if (IsCanvasShown(canvas))
		{
			HideCanvas(canvas);
		}
		else
		{
			ShowCanvas(canvas);
		}
	}
}

void Module::UpdateCanvas()
{
	m_mainCanvas->cd();
	m_mainCanvas->Modified();
	m_mainCanvas->Update();

	m_secondCanvas->cd();
	m_secondCanvas->Modified();
	m_secondCanvas->Update();
}
