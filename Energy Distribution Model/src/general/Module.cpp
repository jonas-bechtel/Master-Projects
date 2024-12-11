#include "pch.h"

#include "Module.h"
#include "EnergyDistribution.h"
#include "PlasmaRateCoefficient.h"
#include "RateCoefficient.h"
#include "CrossSection.h"

int EnergyDistributionModule::s_rebinningFactors[3] = { 10, 10, 10 };

EnergyDistribution EnergyDistributionModule::activeDist;
MCMC* EnergyDistributionModule::mcmc = nullptr;
ElectronBeam* EnergyDistributionModule::eBeam = nullptr;
IonBeam* EnergyDistributionModule::ionBeam = nullptr;
LabEnergies* EnergyDistributionModule::labEnergies = nullptr;
EnergyDistributionManager* EnergyDistributionModule::manager = nullptr;

std::vector<EnergyDistribution> EnergyDistribtionListContainer::energyDistributionList;

std::vector<CrossSection> CrossSectionDeconvolutionModule::crossSectionList;
std::vector<RateCoefficient> CrossSectionDeconvolutionModule::rateCoefficientList;
std::vector<PlasmaRateCoefficient> CrossSectionDeconvolutionModule::plasmaRateCoefficientList;
CrossSectionManager* CrossSectionDeconvolutionModule::CSmanager = nullptr;
RateCoefficientManager* CrossSectionDeconvolutionModule::RCmanager = nullptr;

Window::Window(std::string name)
	: m_name(name)
{
	m_mainCanvas = new TCanvas((m_name + " main canvas").c_str(), (m_name + " main canvas").c_str(), 1200, 500);
	m_mainCanvas->Divide(2, 1);

	m_secondCanvas = new TCanvas((m_name + " analysis canvas").c_str(), (m_name + " analysis canvas").c_str(), 1500, 800);
	m_secondCanvas->Divide(3, 2);

	HideCanvas(m_mainCanvas);
	HideCanvas(m_secondCanvas);
}

void Window::ShowWindow()
{
	//ImGui::begin
	if (ImGui::Begin((m_name + " Window").c_str()))
	{
		ImGui::PushItemWidth(100.0f);
		ShowUI();
		ImGui::PopItemWidth();
		ImGui::Separator();
		ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.0f, 0.2f, 0.5f, 1.0f));     
		ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.0f, 0.3f, 0.6f, 1.0f)); 
		ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.0f, 0.2f, 0.4f, 1.0f));

		ShowHideCanvasButton(m_mainCanvas);
		ShowHideCanvasButton(m_secondCanvas);
		ImGui::PopStyleColor(3);
	}
	ImGui::End();

	UpdateCanvas();
}

Window::~Window()
{
	delete m_mainCanvas;
	delete m_secondCanvas;
}

bool Window::IsCanvasShown(TCanvas* canvas)
{
	return ((TRootCanvas*)canvas->GetCanvasImp())->IsMapped();
}

void Window::ShowCanvas(TCanvas* canvas)
{
	((TRootCanvas*)canvas->GetCanvasImp())->MapRaised();

}

void Window::HideCanvas(TCanvas* canvas)
{
	((TRootCanvas*)canvas->GetCanvasImp())->UnmapWindow();
}

void Window::ShowHideCanvasButton(TCanvas* canvas)
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
			//canvas == m_mainCanvas ? PlotAllOnMainCanvas() : PlotAllOnSecondCanvas();
		}
	}
}

void Window::UpdateCanvas()
{
	m_mainCanvas->cd();
	m_mainCanvas->Modified();
	m_mainCanvas->Update();

	m_secondCanvas->cd();
	m_secondCanvas->Modified();
	m_secondCanvas->Update();
}

EnergyDistributionModule::EnergyDistributionModule(std::string name)
	: Window(name)
{

}

TH3D* EnergyDistributionModule::GetDistribution()
{
	return m_distribution;
}

void EnergyDistributionModule::PlotDistribution()
{
	if (!m_distribution) return;
	if (m_distributionSmall) delete m_distributionSmall;

	m_mainCanvas->cd(2);
	
	m_distributionSmall = (TH3D*)m_distribution->Rebin3D(s_rebinningFactors[0],
		s_rebinningFactors[1],
		s_rebinningFactors[2], "dist small");

	m_distributionSmall->GetXaxis()->SetTitle("x-axis");
	m_distributionSmall->GetYaxis()->SetTitle("y-axis");
	m_distributionSmall->GetZaxis()->SetTitle("z-axis");
	m_distributionSmall->Draw("BOX2 COLZ");
}

EnergyDistributionModule::~EnergyDistributionModule()
{
	delete m_distribution;
	delete m_distributionSmall;
}

bool EnergyDistributionModule::RebinningFactorInput()
{
	return ImGui::InputInt3("Rebinning factors", s_rebinningFactors);
}

CrossSectionDeconvolutionModule::CrossSectionDeconvolutionModule(std::string name)
	: Window(name)
{
}

CrossSectionDeconvolutionModule::~CrossSectionDeconvolutionModule()
{
}
