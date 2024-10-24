#include "IonBeam.h"

#include "imgui.h"

#include <TRootCanvas.h>
#include "ElectronBeam.h"
#include "MCMC.h"

IonBeam::IonBeam()
	: Distribution3D("Ion Beam")
{
}

float IonBeam::GetRadius()
{
	return parameter.radius;
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
				x -= parameter.shift[0];
				y -= parameter.shift[1];

				// Calculate the value using the Gaussian distribution centered at z = 0
				double value = exp(-(x * x + y * y) / (2.0 * parameter.radius * parameter.radius));
				m_distribution->SetBinContent(i, j, k, value);
			}
		}
	}
}

IonBeamParameters IonBeam::GetParameter()
{
	return parameter;
}

void IonBeam::SetParameter(IonBeamParameters params)
{
	parameter = params;
}

TH3D* IonBeam::MultiplyWithElectronDensities()
{
	ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");
	TH3D* electronDensities = eBeam->GetDistribution();

	if (!electronDensities) return nullptr;

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
	ImGui::SetNextItemWidth(200.0f);
	somethingChanged = ImGui::InputFloat("ion beam radius / sigma in [m]", &parameter.radius, 0.001f, 0.001f, "%.4f");
	ImGui::SetNextItemWidth(200.0f);
	somethingChanged |= ImGui::InputFloat2("shift in x and y [m]", parameter.shift, "%.4f");

	if(somethingChanged)
	{
		MCMC* mcmc = (MCMC*)Module::Get("MCMC");
		ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");
		SetupDistribution();
		mcmc->SetupDistribution();

		PlotDistribution();
		mcmc->PlotTargetDistribution();
	}
}

std::string IonBeamParameters::String()
{
	std::string string = std::string(Form("# ion beam parameter:\n")) +
						 std::string(Form("# radius: %.4f m\n", radius)) +
						 std::string(Form("# shift in x and y: %.4f, %4f m\n", shift[0], shift[1]));

	return string;
}
