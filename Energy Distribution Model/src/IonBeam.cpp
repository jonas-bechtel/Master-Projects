#include "IonBeam.h"

#include "imgui.h"

#include <TRootCanvas.h>

IonBeam::IonBeam()
	: Module("Ion Beam")
{
}

bool IonBeam::HasDensityChanged()
{
	return densityChanged;
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
				double value = exp(-(x * x + y * y) / (2.0 * size * size));
				m_distribution->SetBinContent(i, j, k, value);
			}
		}
	}
	PlotDistribution();
}

void IonBeam::ShowUI()
{
	densityChanged = false;
	if (ImGui::SliderFloat("ion beam size / sigma in [m]", &size, 0, 0.01))
	{
		//CreateIonBeam(density);
		densityChanged = true;
	}
}
