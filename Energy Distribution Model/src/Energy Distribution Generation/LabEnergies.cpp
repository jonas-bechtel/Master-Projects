#include "pch.h"

#include "LabEnergies.h"
#include "FileUtils.h"
#include "HistUtils.h"

namespace LabEnergy
{
	LabEnergyParameters parameters;

	// 3D Hist with main data
	static TH3D* hist;

	// optional parameters
	static bool uniformLabEnergies = false;
	static bool noSpaceCharge = false;
	static bool interpolateEnergy = true;

	// plotting things
	static std::vector<HistData3D> plotEnergies;
	static int selectedIndex = -1;

	// z value for the xy slice of the lab energies
	static float SliceZ = 0.0f;
	static bool showMarkers = false;

	void Init()
	{
	}

	TH3D* Get()
	{
		return hist;
	}

	double GetValue(double x, double y, double z)
	{
		z = TMath::Abs(z);

		if(noSpaceCharge)
		{
			double minValue = std::numeric_limits<double>::max();  // start with large number

			int nx = hist->GetNbinsX();
			int ny = hist->GetNbinsY();
			int iz = hist->GetZaxis()->FindBin(z);

			for (int ix = 1; ix <= nx; ++ix) 
			{   
				for (int iy = 1; iy <= ny; ++iy)
				{
					double val = hist->GetBinContent(ix, iy, iz);
					if (val < minValue)
					{
						minValue = val;
					}
				}
			}
			return minValue;
			//return HistUtils::GetValueAtPosition(hist, { 0, 0, z }, interpolateEnergy);
		}

		return HistUtils::GetValueAtPosition(hist, {x, y, z}, interpolateEnergy);
	}

	LabEnergyParameters GetParameters()
	{
		return parameters;
	}

	void SetParameters(const LabEnergyParameters& params)
	{
		parameters = params;
	}

	double GetCenterLabEnergy()
	{
		return parameters.centerLabEnergy;
	}

	void SetDriftTubeVoltage(double voltage)
	{
		parameters.driftTubeVoltage = voltage;
	}

	void SetCenterEnergy(double energy)
	{
		parameters.centerLabEnergy = energy;
	}

	std::string GetTags()
	{
		std::string tags = "";
		if (uniformLabEnergies) tags += "uniform energy, ";
		if (!interpolateEnergy) tags += "no energy interpolation, ";
		if (noSpaceCharge) tags += "no space charge, ";

		return tags;
	}

	void SetupDistribution(std::filesystem::path energyfile)
	{
		delete hist;
		if (uniformLabEnergies && parameters.centerLabEnergy)
		{
			hist = GenerateLabEnergies();
		}
		else
		{
			hist = LoadLabEnergyFile(energyfile);
		}
	}

	TH3D* LoadLabEnergyFile(std::filesystem::path file)
	{
		if (!file.empty())
		{
			TH3D* result = FileUtils::LoadMatrixFile(file);
			result->SetTitle("lab energies");
			result->SetName("lab energies");

			parameters.energyFile.set(file);

			return result;
		}
		return nullptr;
	}

	TH3D* GenerateLabEnergies()
	{
		TH3D* result = new TH3D("uniform energies", "uniform energies", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, 0, 0.7);
		for (int x = 1; x <= result->GetNbinsX(); x++)
		{
			for (int y = 1; y <= result->GetNbinsY(); y++)
			{
				for (int z = 1; z <= result->GetNbinsZ(); z++)
				{
					result->SetBinContent(x, y, z, parameters.centerLabEnergy);
				}
			}
		}
		return result;
	}

	void SelectedItemChanged()
	{
		HistData3D& newlySelected = plotEnergies.at(selectedIndex);
		newlySelected.UpdateSlice(SliceZ);
	}

	void AddBeamToList(HistData3D& beamData)
	{
		plotEnergies.push_back(std::move(beamData));
		if (plotEnergies.size() == 1)
		{
			selectedIndex = 0;
			SelectedItemChanged();
		}
	}

	void RemoveBeamFromList(int index)
	{
		plotEnergies.erase(plotEnergies.begin() + index);
		selectedIndex = std::min(selectedIndex, (int)plotEnergies.size() - 1);

		if (selectedIndex >= 0)
		{
			SelectedItemChanged();
		}
	}

	void ShowWindow()
	{
		if (ImGui::Begin("Lab energy Window"))
		{
			if (ImGui::BeginChild("left side", ImVec2(300, -1), ImGuiChildFlags_ResizeX))
			{
				ShowList();

				if (ImGui::Button("Load lab energies"))
				{
					std::vector<std::filesystem::path> files = FileUtils::SelectFiles(FileUtils::GetDataFolder());
					for (const std::filesystem::path& file : files)
					{
						TH3D* hist = LoadLabEnergyFile(file);

						HistData3D newEnergy(hist);
						std::string index = file.filename().string().substr(0, 4);
						std::string label = file.parent_path().parent_path().filename().string() + ": index " + index;
						newEnergy.SetLabel(label);

						AddBeamToList(newEnergy);
					}
				}
				ImGui::SameLine();
				if (ImGui::Button("clear list"))
				{
					for (int i = plotEnergies.size() - 1; i >= 0; i--)
					{
						RemoveBeamFromList(i);
					}
				}

				ImGui::SetNextItemWidth(200.0f);
				if (ImGui::SliderFloat("slice z", &SliceZ, 0.0f, 0.7f))
				{
					if (selectedIndex >= 0)
					{
						HistData3D& energy = plotEnergies.at(selectedIndex);
						energy.UpdateSlice(SliceZ);
					}
				}

				ImGui::Checkbox("show markers", &showMarkers);
				ImGui::Separator();

				ShowParameterControls();

			}
			ImGui::EndChild();

			ImGui::SameLine();
			ShowPlots();

		}
		ImGui::End();
	}

	void ShowList()
	{
		if (ImGui::BeginListBox("##lab energies", ImVec2(-1, 400.0f)))
		{
			for (int i = 0; i < plotEnergies.size(); i++)
			{
				ImGui::PushID(i);
				HistData3D& le = plotEnergies.at(i);

				if (ImGui::Selectable(le.GetLabel().c_str(), i == selectedIndex, ImGuiSelectableFlags_AllowItemOverlap))
				{
					selectedIndex = i;
					SelectedItemChanged();
				}

				ImGui::SameLine();
				if (ImGui::SmallButton("x"))
				{
					RemoveBeamFromList(i);
				}
				ImGui::PopID();
			}
			ImGui::EndListBox();
		}
	}

	void ShowParameterControls()
	{
		ImGui::BeginGroup();

		ImGui::Checkbox("interpolate", &interpolateEnergy);
		ImGui::Checkbox("no space charge", &noSpaceCharge);
		ImGui::Checkbox("uniform energies", &uniformLabEnergies);
		ImGui::EndGroup();
	}

	void ShowPlots()
	{
		if (ImPlot::BeginSubplots("##labenergy subplots", 2, 3, ImVec2(-1, -1), ImPlotSubplotFlags_ShareItems))
		{
			if (showMarkers) ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Square);

			if (ImPlot::BeginPlot("Projection X"))
			{
				ImPlot::SetupAxes("x", "normalised value");
				for (const HistData3D& energy : plotEnergies)
				{
					energy.PlotProjectionX();
				}
				ImPlot::EndPlot();
			}

			if (ImPlot::BeginPlot("Projection Y"))
			{
				ImPlot::SetupAxes("y", "normalised value");
				for (const HistData3D& energy : plotEnergies)
				{
					energy.PlotProjectionY();
				}
				ImPlot::EndPlot();
			}

			if (ImPlot::BeginPlot("Projection Z"))
			{
				ImPlot::SetupAxes("z", "normalised value");
				for (const HistData3D& energy : plotEnergies)
				{
					energy.PlotProjectionZ();
				}
				ImPlot::EndPlot();
			}

			if (ImPlot::BeginPlot("Inside/Outside"))
			{
				ImPlot::SetupAxes("z", "value");
				for (const HistData3D& energy : plotEnergies)
				{
					energy.PlotInsideOutsideValue();
				}
				ImPlot::EndPlot();
			}
			if (showMarkers) ImPlot::PopStyleVar();

			if (selectedIndex >= 0)
			{
				const HistData3D& sliceLE = plotEnergies.at(selectedIndex);
				sliceLE.PlotSlice();
			}

			ImPlot::EndSubplots();
		}
	}
}
