#include "IonBeam.h"

#include "imgui.h"

#include <TRootCanvas.h>
#include "ElectronBeam.h"
#include "MCMC.h"

IonBeam::IonBeam()
	: Module("Ion Beam")
{
}

float IonBeam::Radius()
{
	return parameter.radius;
}

IonBeamParameters IonBeam::GetParameter()
{
	return parameter;
}

void IonBeam::SetParameter(IonBeamParameters params)
{
	parameter = params;
}

void IonBeam::MultiplyWithElectronDensities(TH3D* electronDensities)
{
	if (!electronDensities) return;

	CreateIonBeam(electronDensities);

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

	MCMC* mcmc = (MCMC*)Module::Get("MCMC");
	mcmc->SetTargetDistribution(result);
}

void IonBeam::CreateIonBeam(TH3D* referenceDensity)
{
	if (!referenceDensity) return;

	if (m_distribution != referenceDensity)
	{
		m_distribution = (TH3D*)referenceDensity->Clone("ion density");
		m_distribution->SetTitle("ion density");
	}
	m_distribution->Reset();
	double nXBins = referenceDensity->GetXaxis()->GetNbins();
	double nYBins = referenceDensity->GetYaxis()->GetNbins();
	double nZBins = referenceDensity->GetZaxis()->GetNbins();

	for (int i = 1; i <= nXBins; i++) {
		for (int j = 1; j <= nYBins; j++) {
			for (int k = 1; k <= nZBins; k++) {
				// Calculate the coordinates for this bin
				double x = m_distribution->GetXaxis()->GetBinCenter(i);
				double y = m_distribution->GetYaxis()->GetBinCenter(j);
				//double z = density->GetZaxis()->GetBinCenter(k);

				// Calculate the value using the Gaussian distribution centered at z = 0
				//double normalisation = 1 / (2 * TMath::Pi() * ionBeamRadius * ionBeamRadius);
				double value = exp(-(x * x + y * y) / (2.0 * parameter.radius * parameter.radius));
				m_distribution->SetBinContent(i, j, k, value);
			}
		}
	}
	PlotDistribution();
}

void IonBeam::ShowUI()
{
	if (ImGui::SliderFloat("ion beam size / sigma in [m]", &parameter.radius, 0, 0.01))
	{
		ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");
		MultiplyWithElectronDensities(eBeam->GetDistribution());
	}
}

std::string IonBeamParameters::String()
{
	std::string string = std::string(Form("# radius: %.3f m\n", radius));

	return string;
}
