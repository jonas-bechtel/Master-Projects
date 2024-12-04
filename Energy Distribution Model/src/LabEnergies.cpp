#include "pch.h"

#include "LabEnergies.h"
#include "FileHandler.h"

LabEnergies::LabEnergies()
	: EnergyDistributionModule("Lab Energies"), m_parameters(activeDist.labEnergiesParameter)
{
	
}

double LabEnergies::Get(double x, double y, double z)
{
	int numberBinsX = m_distribution->GetXaxis()->GetNbins();
	int numberBinsY = m_distribution->GetYaxis()->GetNbins();
	int numberBinsZ = m_distribution->GetZaxis()->GetNbins();
	
	double x_modified = std::min(std::max(x, m_distribution->GetXaxis()->GetBinCenter(1)), m_distribution->GetXaxis()->GetBinCenter(numberBinsX) - 1e-4);
	double y_modified = std::min(std::max(y, m_distribution->GetYaxis()->GetBinCenter(1)), m_distribution->GetYaxis()->GetBinCenter(numberBinsY) - 1e-4);
	double z_modified = std::min(std::max(z, m_distribution->GetZaxis()->GetBinCenter(1)), m_distribution->GetZaxis()->GetBinCenter(numberBinsZ) - 1e-4);

	return m_distribution->Interpolate(x_modified, y_modified, z_modified);
}

void LabEnergies::SetupDistribution(std::filesystem::path energyfile)
{

	if (activeDist.simplifyParams.uniformLabEnergies && m_parameters.centerLabEnergy)
	{
		GenerateUniformLabEnergy();
	}
	else
	{
		LoadLabEnergyFile(energyfile);

		if (activeDist.simplifyParams.sliceLabEnergies)
		{
			FillEnergiesWithXY_Slice();
		}
	}
	
}

void LabEnergies::LoadLabEnergyFile(std::filesystem::path file)
{
	if (!file.empty())
	{
		delete m_distribution;
		m_distribution = FileHandler::GetInstance().LoadMatrixFile(file);
		m_distribution->SetTitle("lab energies");
		m_distribution->SetName("lab energies");

		m_parameters.energyFile.set(file);
	}
}

void LabEnergies::ShowUI()
{
	if (ImGui::Button("Load lab energies"))
	{
		std::filesystem::path file = FileHandler::GetInstance().SelectFile();
		LoadLabEnergyFile(file);
		PlotLabEnergyProjections();
		PlotDistribution();
		PlotLabEnergySlice();
		PlotOutInsideEnergyOnZ();
	}
	ImGui::SameLine();
	if (ImGui::InputFloat("slice z value", &SliceZ, 0.05f))
	{
		PlotLabEnergySlice();
	}
	ImGui::InputDouble("drift tube voltage", m_parameters.driftTubeVoltage, 0.0, 0.0, "%.3f");

	ImGui::BeginDisabled(activeDist.simplifyParams.sliceLabEnergies);
	ImGui::Checkbox("uniform energies", activeDist.simplifyParams.uniformLabEnergies);
	ImGui::EndDisabled();

	ImGui::BeginDisabled(activeDist.simplifyParams.uniformLabEnergies);
	ImGui::Checkbox("fill energies with slice", activeDist.simplifyParams.sliceLabEnergies);
	ImGui::EndDisabled();
	ImGui::SameLine();
	ImGui::BeginDisabled(!activeDist.simplifyParams.sliceLabEnergies);
	ImGui::InputDouble("z slice", activeDist.simplifyParams.sliceToFill, 0.05f);
	ImGui::EndDisabled();

	ImGui::Separator();
	if (ImGui::Button("clear plot"))
	{
		zValues.clear();
		energyValuesInside.clear();
		energyValuesOutside.clear();
	}
	if (ImPlot::BeginPlot("lab energy on trajectory and outside tube"))
	{
		for (int i = 0; i < energyValuesInside.size(); i++)
		{
			ImGui::PushID(i);
			ImVec4 color = ImPlot::GetColormapColor(i % ImPlot::GetColormapSize());
			ImPlot::PushStyleColor(ImPlotCol_Line, color);
			
			//ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, 2);
			ImPlot::PlotLine("inside", zValues.data(), energyValuesInside[i].data(), energyValuesInside[i].size());
			ImPlot::PlotLine("outside", zValues.data(), energyValuesOutside[i].data(), energyValuesOutside[i].size(), ImPlotLineFlags_Segments);
			
			ImPlot::PopStyleColor();
			ImGui::PopID();
		}
		
		ImPlot::EndPlot();
	}

}

void LabEnergies::GenerateUniformLabEnergy()
{
	if (m_distribution) delete m_distribution;

	m_distribution = new TH3D("uniform energies", "uniform energies", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, 0, 0.65);
	for (int x = 1; x <= m_distribution->GetNbinsX(); x++)
	{
		for (int y = 1; y <= m_distribution->GetNbinsY(); y++)
		{
			for (int z = 1; z <= m_distribution->GetNbinsZ(); z++)
			{
				m_distribution->SetBinContent(x, y, z, m_parameters.centerLabEnergy);
			}
		}
	}
}

void LabEnergies::FillEnergiesWithXY_Slice()
{
	if (!m_distribution) return;

	int z_bin = m_distribution->GetZaxis()->FindBin(activeDist.simplifyParams.sliceToFill);

	TH3D* temp = (TH3D*)m_distribution->Clone("slice filled energies");

	for (int x = 1; x <= temp->GetNbinsX(); x++)
	{
		for (int y = 1; y <= temp->GetNbinsY(); y++)
		{
			for (int z = 1; z <= temp->GetNbinsZ(); z++)
			{
				double value = m_distribution->GetBinContent(x, y, z_bin);
				temp->SetBinContent(x, y, z, value);
			}
		}
	}
	delete m_distribution;
	m_distribution = temp;
}

void LabEnergies::PlotLabEnergySlice()
{
	if (!m_distribution) return;
	if (labEnergySliceXY) delete labEnergySliceXY;

	int z_bin = m_distribution->GetZaxis()->FindBin(SliceZ);

	// Create a new TH2D histogram for the slice
	int n_bins_x = m_distribution->GetNbinsX();
	int n_bins_y = m_distribution->GetNbinsY();
	double x_min = m_distribution->GetXaxis()->GetXmin();
	double x_max = m_distribution->GetXaxis()->GetXmax();
	double y_min = m_distribution->GetYaxis()->GetXmin();
	double y_max = m_distribution->GetYaxis()->GetXmax();

	labEnergySliceXY = new TH2D("Energy Slice", Form("XY Slice at Z = %.2f", SliceZ),
		n_bins_x, x_min, x_max,
		n_bins_y, y_min, y_max);

	// Fill the 2D histogram with the contents of the corresponding Z slice
	for (int x_bin = 1; x_bin <= n_bins_x; x_bin++)
	{
		for (int y_bin = 1; y_bin <= n_bins_y; y_bin++)
		{
			double content = m_distribution->GetBinContent(x_bin, y_bin, z_bin);
			labEnergySliceXY->SetBinContent(x_bin, y_bin, content);
		}
	}
	labEnergySliceXY->SetContour(100); 

	m_mainCanvas->cd(1)->SetRightMargin(0.15);
	labEnergySliceXY->Draw("COLZ");
}

void LabEnergies::PlotLabEnergyProjections()
{
	if (!m_distribution) return;

	delete labEnergyProjectionX;
	delete labEnergyProjectionY;
	delete labEnergyProjectionZ;

	m_secondCanvas->cd(1);
	labEnergyProjectionX = m_distribution->ProjectionX();
	labEnergyProjectionX->Draw();

	m_secondCanvas->cd(2);
	labEnergyProjectionY = m_distribution->ProjectionY();
	labEnergyProjectionY->Draw();

	m_secondCanvas->cd(3);
	labEnergyProjectionZ = m_distribution->ProjectionZ();
	labEnergyProjectionZ->Draw();
}

void LabEnergies::PlotOutInsideEnergyOnZ()
{
	if (!m_distribution) return;
	if (labEnergyInside) delete labEnergyInside;
	if (labEnergyOutside) delete labEnergyOutside;

	ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");

	labEnergyInside = new TGraph(m_distribution->GetNbinsZ());
	labEnergyOutside = new TGraph(m_distribution->GetNbinsZ());

	std::vector<double> insideEnergies;
	std::vector<double> outsideEnergies;

	int binInCenterX = m_distribution->GetNbinsX() / 2;
	//int binInCenterY = m_distribution->GetNbinsY() / 2;
	for (int i = 1; i <= m_distribution->GetNbinsZ(); i++)
	{
		double zValue = m_distribution->GetZaxis()->GetBinCenter(i);
		int binInCenterY = m_distribution->GetYaxis()->FindBin(eBeam->Trajectory(zValue));
		//std::cout << binInCenterY << std::endl;
		double energyValueIn = m_distribution->GetBinContent(binInCenterX, binInCenterY, i);
		double energyValueOut = m_distribution->GetBinContent(1, 1, i);
		//std::cout << zValue << " " << energyValueIn << " " << energyValueOut << std::endl;
		labEnergyInside->SetPoint(i - 1, zValue, energyValueIn);
		labEnergyOutside->SetPoint(i - 1, zValue, energyValueOut);

		if (counter == 0)
		{
			zValues.push_back(zValue);
		}
		insideEnergies.push_back(energyValueIn);
		outsideEnergies.push_back(energyValueOut);
	}

	energyValuesInside.push_back(insideEnergies);
	energyValuesOutside.push_back(outsideEnergies);
	counter++;

	m_secondCanvas->cd(6);
	labEnergyInside->SetLineColor(kRed);
	labEnergyInside->SetTitle("lab energy on e-beam trajectory and on the outside");
	labEnergyInside->GetXaxis()->SetTitle("z");
	labEnergyInside->GetYaxis()->SetTitle("energy");
	
	labEnergyInside->Draw("ALP");
	labEnergyOutside->Draw("LP");
	
}

