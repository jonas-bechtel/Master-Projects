#include "pch.h"

#include "IonBeam.h"
#include "ElectronBeam.h"
#include "MCMC.h"
#include "EnergyDistribution.h"

IonBeamWindow::IonBeamWindow()
	: EnergyDistributionModule("Ion Beam", 0), m_parameters(activeDist.ionBeamParameter)
{
	m_distribution = new TH3D("hist to look at", "hist to look at", 100, -0.04, 0.04, 100, -0.04, 0.04, 100, 0.0, 0.7);
	UpdateDataToLookAt();
	ionBeam = this;
}

void IonBeamWindow::SetupDistribution(std::filesystem::path file)
{
	TH3D* electronBeam = eBeam->GetDistribution();

	if (!electronBeam) return;

	if (m_distribution != electronBeam)
	{
		delete m_distribution;
		m_distribution = (TH3D*)electronBeam->Clone("ion density");
		m_distribution->SetTitle("ion density");
	}
	m_distribution->Reset();
	FillHistogram(m_distribution);
}

TH3D* IonBeamWindow::MultiplyWithElectronDensities()
{
	TH3D* electronDensities = eBeam->GetDistribution();
	if (!electronDensities)
	{
		std::cout << "no electron densities" << std::endl;
		return nullptr;
	}
	if (!(electronDensities->GetNbinsX() == m_distribution->GetNbinsX()) ||
		!(electronDensities->GetNbinsY() == m_distribution->GetNbinsY()) ||
		!(electronDensities->GetNbinsZ() == m_distribution->GetNbinsZ())) 
	{
		std::cout << "Histograms must have the same binning and ranges!" << std::endl;
		return nullptr;
	}

	TH3D* result = (TH3D*)electronDensities->Clone("e*ion density");
	result->SetTitle("electron-ion density");
	result->Multiply(m_distribution);
	//result->Reset();

	//double nXBins = electronDensities->GetXaxis()->GetNbins();
	//double nYBins = electronDensities->GetYaxis()->GetNbins();
	//double nZBins = electronDensities->GetZaxis()->GetNbins();
	//
	//for (int i = 1; i <= nXBins; i++)
	//{
	//	for (int j = 1; j <= nYBins; j++)
	//	{
	//		for (int k = 1; k <= nZBins; k++)
	//		{
	//			double value = electronDensities->GetBinContent(i, j, k) * m_distribution->GetBinContent(i, j, k);
	//			result->SetBinContent(i, j, k, value);
	//		}
	//	}
	//}

	return result;
}

void IonBeamWindow::ShowUI()
{
	if (ImGui::BeginChild("settings", ImVec2(300, -1), ImGuiChildFlags_ResizeX))
	{
		ShowSettings();
		ImGui::EndChild();
	}
	ImGui::SameLine();
	ShowPlots();
}

void IonBeamWindow::ShowSettings()
{
	bool somethingChanged = false;
	somethingChanged |= ImGui::Checkbox("use single gaussian", activeDist.simplifyParams.singleGaussianIonBeam);
	ImGui::BeginDisabled(!activeDist.simplifyParams.singleGaussianIonBeam);
	somethingChanged |= ImGui::InputDouble("ion beam radius / sigma in [m]", activeDist.simplifyParams.ionBeamRadius, 0.0f, 0.0f, "%.6f", ImGuiInputTextFlags_EnterReturnsTrue);
	ImGui::EndDisabled();

	ImGui::Separator();
	ImGui::SetNextItemWidth(200.0f);
	somethingChanged |= ImGui::InputFloat2("shift in x and y [m]", m_parameters.shift, "%.4f");
	ImGui::BeginDisabled(activeDist.simplifyParams.singleGaussianIonBeam);
	somethingChanged |= ImGui::InputDouble("amplitude 1", m_parameters.amplitude1, 0.0f, 0.0f, "%.4f", ImGuiInputTextFlags_EnterReturnsTrue);
	somethingChanged |= ImGui::InputDouble("amplitude 2", m_parameters.amplitude2, 0.0f, 0.0f, "%.4f", ImGuiInputTextFlags_EnterReturnsTrue);
	ImGui::SetNextItemWidth(200.0f);
	somethingChanged |= ImGui::InputFloat2("sigmas 1 x and y [m]", m_parameters.shape1, "%.4f", ImGuiInputTextFlags_EnterReturnsTrue);
	ImGui::SetNextItemWidth(200.0f);
	somethingChanged |= ImGui::InputFloat2("sigmas 2 x and y [m]", m_parameters.shape2, "%.4f", ImGuiInputTextFlags_EnterReturnsTrue);
	ImGui::EndDisabled();

	if (ImGui::SliderFloat("slice z", &SliceZ, 0, 0.7))
	{
		slice.FromTH3D(m_distribution, SliceZ);
	}

	if (somethingChanged)
	{
		UpdateDataToLookAt();
		//SetupDistribution();
		//mcmc->SetupDistribution();
		//
		//PlotDistribution();
		//mcmc->PlotTargetDistribution();

		//PlotIonBeamProjections();
	}
}

void IonBeamWindow::ShowPlots()
{
	if (ImPlot::BeginSubplots("##labenergy subplots", 2, 3, ImVec2(-1, -1), ImPlotSubplotFlags_ShareItems))
	{
		if (ImPlot::BeginPlot("Projection X"))
		{
			ImPlot::PlotLine("", xAxis.data(), projectionValuesX.data(), xAxis.size());
			ImPlot::EndPlot();
		}

		if (ImPlot::BeginPlot("Projection Y"))
		{
			ImPlot::PlotLine("", yAxis.data(), projectionValuesY.data(), yAxis.size());
			
			ImPlot::EndPlot();
		}

		if (ImPlot::BeginPlot("Projection Z"))
		{
			
			ImPlot::PlotLine("", zAxis.data(), projectionValuesZ.data(), zAxis.size());
			
			ImPlot::EndPlot();
		}

		ImPlot::PushColormap(9);
		if (ImPlot::BeginPlot("XY Slice"))
		{
			ImPlot::SetupAxes("x", "y");
			ImPlot::PlotHeatmap("", slice.values.data(), slice.nRows,
				slice.nCols, slice.minValue, slice.maxValue, nullptr,
				slice.bottomLeft, slice.topRight);
			ImPlot::EndPlot();
		}
		ImGui::SameLine();
		ImPlot::ColormapScale("##HeatScale", slice.minValue, slice.maxValue, ImVec2(100, -1));

		ImPlot::PopColormap();

		ImPlot::EndSubplots();
	}
}

void IonBeamWindow::FillHistogram(TH3D* hist)
{
	int nXBins = hist->GetXaxis()->GetNbins();
	int nYBins = hist->GetYaxis()->GetNbins();
	int nZBins = hist->GetZaxis()->GetNbins();

	for (int i = 1; i <= nXBins; i++) {
		for (int j = 1; j <= nYBins; j++) {
			for (int k = 1; k <= nZBins; k++) {
				// Calculate the coordinates for this bin
				double x = hist->GetXaxis()->GetBinCenter(i);
				double y = hist->GetYaxis()->GetBinCenter(j);
				//double z = density->GetZaxis()->GetBinCenter(k);

				// apply shift of ion beam
				x -= m_parameters.shift.get().x;
				y -= m_parameters.shift.get().y;

				double value = 0;
				// Calculate the value using a single Gaussian distribution centered at z = 0
				if (activeDist.simplifyParams.singleGaussianIonBeam)
				{
					value = exp(-(x * x + y * y) / (2.0 * pow(activeDist.simplifyParams.ionBeamRadius, 2)));
				}
				// or use the double gaussian version
				else
				{
					value = m_parameters.amplitude1 * exp(-0.5 * ((x * x) / pow(m_parameters.shape1.get().x, 2) + (y * y) / pow(m_parameters.shape1.get().y, 2))) +
						m_parameters.amplitude2 * exp(-0.5 * ((x * x) / pow(m_parameters.shape2.get().x, 2) + (y * y) / pow(m_parameters.shape2.get().y, 2)));
				}
				hist->SetBinContent(i, j, k, value);
			}
		}
	}
}

void IonBeamWindow::UpdateDataToLookAt()
{
	xAxis.clear();
	yAxis.clear();
	zAxis.clear();

	projectionValuesX.clear();
	projectionValuesY.clear();
	projectionValuesZ.clear();

	//slice.clear();

	m_distribution->Reset();
	FillHistogram(m_distribution);

	xAxis.reserve(m_distribution->GetNbinsX());
	yAxis.reserve(m_distribution->GetNbinsY());
	zAxis.reserve(m_distribution->GetNbinsZ());

	projectionValuesX.reserve(m_distribution->GetNbinsX());
	projectionValuesY.reserve(m_distribution->GetNbinsY());
	projectionValuesZ.reserve(m_distribution->GetNbinsZ());

	for (int i = 1; i <= m_distribution->GetNbinsX(); i++)
	{
		xAxis.push_back(m_distribution->GetXaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= m_distribution->GetNbinsY(); i++)
	{
		yAxis.push_back(m_distribution->GetYaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= m_distribution->GetNbinsZ(); i++)
	{
		zAxis.push_back(m_distribution->GetZaxis()->GetBinCenter(i));
	}

	TH1D* projectionX = m_distribution->ProjectionX();
	TH1D* projectionY = m_distribution->ProjectionY();
	TH1D* projectionZ = m_distribution->ProjectionZ();

	for (int i = 1; i <= projectionX->GetNbinsX(); i++)
	{
		projectionValuesX.push_back(projectionX->GetBinContent(i));
	}
	for (int i = 1; i <= projectionY->GetNbinsX(); i++)
	{
		projectionValuesY.push_back(projectionY->GetBinContent(i));
	}
	for (int i = 1; i <= projectionZ->GetNbinsX(); i++)
	{
		projectionValuesZ.push_back(projectionZ->GetBinContent(i));
	}

	delete projectionX;
	delete projectionY;
	delete projectionZ;

	slice.FromTH3D(m_distribution, SliceZ);
}

//void IonBeamWindow::PlotIonBeamProjections()
//{
//	if (!m_distribution) return;
//
//	delete ionBeamProjectionX;
//	delete ionBeamProjectionY;
//	delete ionBeamProjectionZ;
//
//	m_secondCanvas->cd(1);
//	ionBeamProjectionX = m_distribution->ProjectionX();
//	ionBeamProjectionX->Draw();
//
//	m_secondCanvas->cd(2);
//	ionBeamProjectionY = m_distribution->ProjectionY();
//	ionBeamProjectionY->Draw();
//
//	m_secondCanvas->cd(3);
//	ionBeamProjectionZ = m_distribution->ProjectionZ();
//	ionBeamProjectionZ->Draw();
//}

