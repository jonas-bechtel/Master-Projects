#include "IonBeam.h"

#include "imgui.h"

#include <TRootCanvas.h>
#include "ElectronBeam.h"
#include "MCMC.h"

IonBeam::IonBeam()
	: Distribution3D("Ion Beam")
{
}

void IonBeam::SetupDistribution(std::filesystem::path file)
{
	ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");
	TH3D* electronBeam = eBeam->GetDistribution();

	if (!electronBeam) return;

	if (m_distribution != electronBeam)
	{
		m_distribution = (TH3D*)electronBeam->Clone("ion density");
		m_distribution->SetTitle("ion density");
	}
	m_distribution->Reset();
	int nXBins = electronBeam->GetXaxis()->GetNbins();
	int nYBins = electronBeam->GetYaxis()->GetNbins();
	int nZBins = electronBeam->GetZaxis()->GetNbins();

	for (int i = 1; i <= nXBins; i++) {
		for (int j = 1; j <= nYBins; j++) {
			for (int k = 1; k <= nZBins; k++) {
				// Calculate the coordinates for this bin
				double x = m_distribution->GetXaxis()->GetBinCenter(i);
				double y = m_distribution->GetYaxis()->GetBinCenter(j);
				//double z = density->GetZaxis()->GetBinCenter(k);

				// apply shift of ion beam
				x -= m_parameters.shift.get().x;
				y -= m_parameters.shift.get().y;

				double value = 0;
				// Calculate the value using a single Gaussian distribution centered at z = 0
				if (m_parameters.useSingleGaussian)
				{
					value = exp(-(x * x + y * y) / (2.0 * m_parameters.radius * m_parameters.radius));
				}
				// or use the double gaussian version
				else
				{
					value = m_parameters.amplitude1 * exp(-0.5 * ( (x * x) / pow(m_parameters.shape1.get().x, 2) + (y * y) / pow(m_parameters.shape1.get().y, 2))) + 
							m_parameters.amplitude2 * exp(-0.5 * ( (x * x) / pow(m_parameters.shape2.get().x, 2) + (y * y) / pow(m_parameters.shape2.get().y, 2)));
				}
				m_distribution->SetBinContent(i, j, k, value);
			}
		}
	}
}

IonBeamParameters& IonBeam::GetParameter()
{
	return m_parameters;
}

TH3D* IonBeam::MultiplyWithElectronDensities()
{
	ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");
	TH3D* electronDensities = eBeam->GetDistribution();

	if (!electronDensities)
	{
		std::cout << "no electron densities" << std::endl;
		return nullptr;
	}

	TH3D* result = (TH3D*)electronDensities->Clone("e-ion density");
	result->SetTitle("electron-ion density");
	result->Reset();

	double nXBins = electronDensities->GetXaxis()->GetNbins();
	double nYBins = electronDensities->GetYaxis()->GetNbins();
	double nZBins = electronDensities->GetZaxis()->GetNbins();

	for (int i = 1; i <= nXBins; i++)
	{
		for (int j = 1; j <= nYBins; j++)
		{
			for (int k = 1; k <= nZBins; k++)
			{
				double value = electronDensities->GetBinContent(i, j, k) * m_distribution->GetBinContent(i, j, k);
				result->SetBinContent(i, j, k, value);
			}
		}
	}

	return result;
}

void IonBeam::ShowUI()
{
	bool somethingChanged = false;
	ImGui::Checkbox("use single gaussian", m_parameters.useSingleGaussian);
	ImGui::BeginDisabled(!m_parameters.useSingleGaussian);
	ImGui::SetNextItemWidth(100.0f);
	somethingChanged |= ImGui::InputDouble("ion beam radius / sigma in [m]", m_parameters.radius, 0.001f, 0.001f, "%.4f", ImGuiInputTextFlags_EnterReturnsTrue);
	ImGui::EndDisabled();

	ImGui::Separator();
	ImGui::SetNextItemWidth(200.0f);
	somethingChanged |= ImGui::InputFloat2("shift in x and y [m]", m_parameters.shift, "%.4f");													ImGui::SetNextItemWidth(100.0f);
	ImGui::BeginDisabled(m_parameters.useSingleGaussian);
	somethingChanged |= ImGui::InputDouble("amplitude 1", m_parameters.amplitude1, 0.0f, 0.0f, "%.4f", ImGuiInputTextFlags_EnterReturnsTrue);			ImGui::SetNextItemWidth(100.0f);
	somethingChanged |= ImGui::InputDouble("amplitude 2", m_parameters.amplitude2, 0.0f, 0.0f, "%.4f", ImGuiInputTextFlags_EnterReturnsTrue);	ImGui::SetNextItemWidth(200.0f);
	somethingChanged |= ImGui::InputFloat2("sigmas 1 x and y [m]", m_parameters.shape1, "%.4f", ImGuiInputTextFlags_EnterReturnsTrue);			ImGui::SetNextItemWidth(200.0f);
	somethingChanged |= ImGui::InputFloat2("sigmas 2 x and y [m]", m_parameters.shape2, "%.4f", ImGuiInputTextFlags_EnterReturnsTrue);
	ImGui::EndDisabled();

	if(somethingChanged)
	{
		MCMC* mcmc = (MCMC*)Module::Get("MCMC");
		SetupDistribution();
		mcmc->SetupDistribution();
	
		PlotDistribution();
		mcmc->PlotTargetDistribution();

		PlotIonBeamProjections();
	}
}

void IonBeam::PlotIonBeamProjections()
{
	if (!m_distribution) return;

	delete ionBeamProjectionX;
	delete ionBeamProjectionY;
	delete ionBeamProjectionZ;

	m_secondCanvas->cd(1);
	ionBeamProjectionX = m_distribution->ProjectionX();
	ionBeamProjectionX->Draw();

	m_secondCanvas->cd(2);
	ionBeamProjectionY = m_distribution->ProjectionY();
	ionBeamProjectionY->Draw();

	m_secondCanvas->cd(3);
	ionBeamProjectionZ = m_distribution->ProjectionZ();
	ionBeamProjectionZ->Draw();
}

