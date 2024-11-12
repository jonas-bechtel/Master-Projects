#include "LabEnergies.h"
#include "FileHandler.h"

LabEnergies::LabEnergies()
	: Distribution3D("Lab Energies")
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

LabEnergiesParameters LabEnergies::GetParameter()
{
	return parameter;
}

void LabEnergies::SetCenterLabEnergy(double energy)
{
	m_parameters.Set("center lab energy", energy);
	//parameter.centerLabEnergy = energy;
}

void LabEnergies::SetupDistribution(std::filesystem::path energyfile)
{
	bool bla = m_parameters.Get("center lab energy");

	if (m_parameters.Get("use uniform energy") && (double)m_parameters.Get("center lab energy"))
	{
		GenerateUniformLabEnergy();
	}
	else
	{
		LoadLabEnergyFile(energyfile);

		if (m_parameters.Get("use xy slice"))
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

		m_parameters.Set("energy file", file);
		//parameter.energyFile = file;
	}
}

void LabEnergies::ShowUI()
{
	if (ImGui::Button("Load lab energies"))
	{
		std::filesystem::path file = FileHandler::GetInstance().OpenFileExplorer();
		LoadLabEnergyFile(file);
		PlotLabEnergyProjections();
		PlotDistribution();
		PlotLabEnergySlice();
	}
	ImGui::SameLine();
	ImGui::SetNextItemWidth(100.0f);
	if (ImGui::InputFloat("slice z value", &SliceZ, 0.05f))
	{
		PlotLabEnergySlice();
	}
	
	ImGui::BeginDisabled(m_parameters.Get("use xy slice"));
	ImGui::Checkbox("uniform energies", (bool*)&m_parameters.Get("use uniform energy"));
	ImGui::EndDisabled();

	ImGui::BeginDisabled(m_parameters.Get("use uniform energy"));
	ImGui::Checkbox("fill energies with slice", (bool*)&m_parameters.Get("use xy slice"));
	ImGui::EndDisabled();
	ImGui::SameLine();
	ImGui::BeginDisabled(!m_parameters.Get("use xy slice"));
	ImGui::SetNextItemWidth(100.0f);
	//ImGui::InputDouble("z slice", &m_parameters.Get("fill slice"), 0.05f);
	ImGui::EndDisabled();
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
				m_distribution->SetBinContent(x, y, z, m_parameters.Get("center lab energy"));
			}
		}
	}
}

void LabEnergies::FillEnergiesWithXY_Slice()
{
	if (!m_distribution) return;

	int z_bin = m_distribution->GetZaxis()->FindBin(m_parameters.Get("fill slice"));

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

	m_mainCanvas->cd(1);
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

std::string LabEnergiesParameters::String()
{
	std::string string = std::string(Form("# lab energy parameter:\n"));

	//if (m_parameters.Get<bool>("use uniform energy"))
	//{
	//	string += std::string(Form("# using a uniform energy distribution with the center lab energy\n"));
	//}
	//else
	//{
	//	string += std::string(Form("# energy file: %s\n", energyFile.filename().string().c_str()));
	//	if (useOnlySliceXY)
	//	{
	//		string += std::string(Form("# using only one xy slice to fill the rest\n"));
	//		string += std::string(Form("# z value of slice: %.3f\n", sliceToFill));
	//	}
	//}
	//string += std::string(Form("# lab energy in center: %e eV\n", centerLabEnergy));

	return string;
}
