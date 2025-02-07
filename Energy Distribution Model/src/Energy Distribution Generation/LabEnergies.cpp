#include "pch.h"

#include "LabEnergies.h"
#include "FileHandler.h"

LabEnergyWindow::LabEnergyWindow()
	: EnergyDistributionModule("Lab Energies", 0), m_parameters(activeDist.labEnergiesParameter)
{
	labEnergies = this;
}

double LabEnergyWindow::Get(double x, double y, double z)
{
	int numberBinsX = m_distribution->GetXaxis()->GetNbins();
	int numberBinsY = m_distribution->GetYaxis()->GetNbins();
	int numberBinsZ = m_distribution->GetZaxis()->GetNbins();
	
	double x_modified = std::min(std::max(x, m_distribution->GetXaxis()->GetBinCenter(1)), m_distribution->GetXaxis()->GetBinCenter(numberBinsX) - 1e-4);
	double y_modified = std::min(std::max(y, m_distribution->GetYaxis()->GetBinCenter(1)), m_distribution->GetYaxis()->GetBinCenter(numberBinsY) - 1e-4);
	double z_modified = std::min(std::max(z, m_distribution->GetZaxis()->GetBinCenter(1)), m_distribution->GetZaxis()->GetBinCenter(numberBinsZ) - 1e-4);

	return m_distribution->Interpolate(x_modified, y_modified, z_modified);
}

void LabEnergyWindow::SetupDistribution(std::filesystem::path energyfile)
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

void LabEnergyWindow::LoadLabEnergyFile(std::filesystem::path file)
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

void LabEnergyWindow::LoadToLookAt(std::filesystem::path file)
{
	if (!file.empty())
	{
		LabEnergy newLabEnergy;

		newLabEnergy.fullHistogram = FileHandler::GetInstance().LoadMatrixFile(file);
		newLabEnergy.fullHistogram->SetTitle("lab energies");
		newLabEnergy.fullHistogram->SetName("lab energies");
		newLabEnergy.FillData(eBeam);
		newLabEnergy.label = file.parent_path().parent_path().filename().string() + ": index " + std::stoi(file.filename().string().substr(0, 4));

		labEnergiesToLookAt.push_back(std::move(newLabEnergy));
		if (labEnergiesToLookAt.size() == 1)
		{
			selectedIndex = 0;
			LabEnergy& first = labEnergiesToLookAt.at(0);
			first.slice.FromTH3D(first.fullHistogram, SliceZ);
		}
	}
}

void LabEnergyWindow::ShowUI()
{
	ImGui::BeginGroup();
	ShowList();
	ShowSettings();
	ImGui::EndGroup();

	ImGui::SameLine();

	ShowLabEnergyPlots();
}

void LabEnergyWindow::ShowList()
{
	if (ImGui::BeginListBox("##lab energies", ImVec2(250.0f, 400.0f)))
	{
		for (int i = 0; i < labEnergiesToLookAt.size(); i++)
		{
			ImGui::PushID(i);
			LabEnergy& le = labEnergiesToLookAt.at(i);

			if (ImGui::Selectable(le.label.c_str(), i == selectedIndex, ImGuiSelectableFlags_AllowItemOverlap))
			{
				selectedIndex = i;
				le.slice.FromTH3D(le.fullHistogram, SliceZ);
			}
			
			ImGui::SameLine();
			if (ImGui::SmallButton("x"))
			{
				labEnergiesToLookAt.erase(labEnergiesToLookAt.begin() + i);
				selectedIndex = std::min(selectedIndex, (int)labEnergiesToLookAt.size() - 1);

				if (selectedIndex >= 0)
				{
					LabEnergy& newSelected = labEnergiesToLookAt.at(selectedIndex);
					newSelected.slice.FromTH3D(newSelected.fullHistogram, SliceZ);
				}
			}
			ImGui::PopID();
		}
		ImGui::EndListBox();
	}
}

void LabEnergyWindow::ShowSettings()
{
	if (ImGui::Button("Load lab energies"))
	{
		std::vector<std::filesystem::path> files = FileHandler::GetInstance().SelectFiles();
		for (const std::filesystem::path& file : files)
		{
			LoadToLookAt(file);
		}
	}

	ImGui::SetNextItemWidth(200.0f);
	if (ImGui::SliderFloat("slice z", &SliceZ, 0, 0.7))
	{
		if (selectedIndex >= 0)
		{
			LabEnergy& le = labEnergiesToLookAt.at(selectedIndex);
			le.slice.FromTH3D(le.fullHistogram, SliceZ);
		}
	}

	ImGui::Separator();
	ImGui::BeginDisabled(activeDist.simplifyParams.sliceLabEnergies);
	ImGui::Checkbox("uniform energies", activeDist.simplifyParams.uniformLabEnergies);
	ImGui::EndDisabled();

	ImGui::BeginDisabled(activeDist.simplifyParams.uniformLabEnergies);
	ImGui::Checkbox("fill energies with slice", activeDist.simplifyParams.sliceLabEnergies);
	ImGui::EndDisabled();
	ImGui::SameLine();
	ImGui::BeginDisabled(!activeDist.simplifyParams.sliceLabEnergies);
	ImGui::SetNextItemWidth(70.0f);
	ImGui::InputDouble("##z slice", activeDist.simplifyParams.sliceToFill, 0,0, "%.3f");
	ImGui::EndDisabled();
}

void LabEnergyWindow::ShowLabEnergyPlots()
{
	if(ImPlot::BeginSubplots("##labenergy subplots", 2, 3, ImVec2(-1, -1), ImPlotSubplotFlags_ShareItems))
	{
		if (ImPlot::BeginPlot("Projection X"))
		{
			for (const LabEnergy& labEnergy : labEnergiesToLookAt)
			{
				ImPlot::PlotLine(labEnergy.label.c_str(), labEnergy.xAxis.data(), labEnergy.projectionValuesX.data(), labEnergy.xAxis.size());
			}
			ImPlot::EndPlot();
		}

		if (ImPlot::BeginPlot("Projection Y"))
		{
			for (const LabEnergy& labEnergy : labEnergiesToLookAt)
			{
				ImPlot::PlotLine(labEnergy.label.c_str(), labEnergy.yAxis.data(), labEnergy.projectionValuesY.data(), labEnergy.yAxis.size());
			}
			ImPlot::EndPlot();
		}

		if (ImPlot::BeginPlot("Projection Z"))
		{
			for (const LabEnergy& labEnergy : labEnergiesToLookAt)
			{
				ImPlot::PlotLine(labEnergy.label.c_str(), labEnergy.zAxis.data(), labEnergy.projectionValuesZ.data(), labEnergy.zAxis.size());
			}
			ImPlot::EndPlot();
		}

		if (ImPlot::BeginPlot("Inside/Outside"))
		{
			for (const LabEnergy& labEnergy : labEnergiesToLookAt)
			{
				ImPlot::PlotLine(labEnergy.label.c_str(), labEnergy.zAxis.data(), labEnergy.labEnergyInside.data(), labEnergy.zAxis.size());
				ImPlot::PlotLine(labEnergy.label.c_str(), labEnergy.zAxis.data(), labEnergy.labEnergyOutside.data(), labEnergy.zAxis.size(), ImPlotLineFlags_Segments);
			}
			ImPlot::EndPlot();
		}

		ImPlot::PushColormap(9);
		if (selectedIndex >= 0)
		{
			const LabEnergy& sliceLE = labEnergiesToLookAt.at(selectedIndex);
		
			if (ImPlot::BeginPlot("XY Slice"))
			{
				ImPlot::PlotHeatmap(sliceLE.label.c_str(), sliceLE.slice.values.data(), sliceLE.slice.nRows,
					sliceLE.slice.nCols, sliceLE.slice.minValue, sliceLE.slice.maxValue, nullptr,
					sliceLE.slice.bottomLeft, sliceLE.slice.topRight);
				ImPlot::EndPlot();
			}
			ImGui::SameLine();
			ImPlot::ColormapScale("##HeatScale", sliceLE.slice.minValue, sliceLE.slice.maxValue, ImVec2(100, -1));
			
		}
		ImPlot::PopColormap();

		ImPlot::EndSubplots();
	}
}

void LabEnergyWindow::GenerateUniformLabEnergy()
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

void LabEnergyWindow::FillEnergiesWithXY_Slice()
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

void LabEnergy::FillData(const ElectronBeamWindow* eBeam)
{
	xAxis.reserve(fullHistogram->GetNbinsX());
	yAxis.reserve(fullHistogram->GetNbinsY());
	zAxis.reserve(fullHistogram->GetNbinsZ());

	projectionValuesX.reserve(fullHistogram->GetNbinsX());
	projectionValuesY.reserve(fullHistogram->GetNbinsY());
	projectionValuesZ.reserve(fullHistogram->GetNbinsZ());

	for (int i = 1; i <= fullHistogram->GetNbinsX(); i++)
	{
		xAxis.push_back(fullHistogram->GetXaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= fullHistogram->GetNbinsY(); i++)
	{
		yAxis.push_back(fullHistogram->GetYaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= fullHistogram->GetNbinsZ(); i++)
	{
		zAxis.push_back(fullHistogram->GetZaxis()->GetBinCenter(i));
	}

	TH1D* projectionX = fullHistogram->ProjectionX();
	TH1D* projectionY = fullHistogram->ProjectionY();
	TH1D* projectionZ = fullHistogram->ProjectionZ();

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

	int binInCenterX = fullHistogram->GetNbinsX() / 2;
	for (int i = 1; i <= fullHistogram->GetNbinsZ(); i++)
	{
		double zValue = fullHistogram->GetZaxis()->GetBinCenter(i);
		int binInCenterY = fullHistogram->GetYaxis()->FindBin(eBeam->Trajectory(zValue));

		double energyValueIn = fullHistogram->GetBinContent(binInCenterX, binInCenterY, i);
		double energyValueOut = fullHistogram->GetBinContent(1, 1, i);

		labEnergyInside.push_back(energyValueIn);
		labEnergyOutside.push_back(energyValueOut);
	}
}

LabEnergy::LabEnergy(LabEnergy&& other)
{
	fullHistogram = other.fullHistogram;
	other.fullHistogram = nullptr;

	xAxis = std::move(other.xAxis);
	yAxis = std::move(other.yAxis);
	zAxis = std::move(other.zAxis);

	projectionValuesX = std::move(other.projectionValuesX);
	projectionValuesY = std::move(other.projectionValuesY);
	projectionValuesZ = std::move(other.projectionValuesZ);

	slice = std::move(other.slice);

	labEnergyInside = std::move(other.labEnergyInside);
	labEnergyOutside = std::move(other.labEnergyOutside);

	label = std::move(other.label);
	//std::cout << "Lab energy move Constructor" << std::endl;
}

LabEnergy& LabEnergy::operator=(LabEnergy&& other)
{
	if (this == &other) return *this;

	delete fullHistogram;
	fullHistogram = other.fullHistogram;
	other.fullHistogram = nullptr;

	xAxis = std::move(other.xAxis);
	yAxis = std::move(other.yAxis);
	zAxis = std::move(other.zAxis);

	projectionValuesX = std::move(other.projectionValuesX);
	projectionValuesY = std::move(other.projectionValuesY);
	projectionValuesZ = std::move(other.projectionValuesZ);

	slice = std::move(other.slice);

	labEnergyInside = std::move(other.labEnergyInside);
	labEnergyOutside = std::move(other.labEnergyOutside);

	label = std::move(other.label);
	//std::cout << "Lab energy move assignemnt" << std::endl;
	return *this;
}
