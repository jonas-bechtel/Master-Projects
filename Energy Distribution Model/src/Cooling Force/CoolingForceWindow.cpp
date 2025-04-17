#include "pch.h"
#include "CoolingForceWindow.h"
#include "CoolingForceModel.h"
#include "CoolingForceCurve.h"
#include "FileUtils.h"

namespace CoolingForce
{
	// main storge data
	static std::vector<Curve> curveList;
	static int currentCurveIndex = -1;

	// optional parameter
	static int method = 3;
	static bool interpolate = true;

	// currently loaded description file
	static std::filesystem::path currentDescriptionFile;
	static int maxIndex = 0;

	// start/end index in description file to generate cooling force for
	static int startIndex = 1;
	static int endIndex = 1;
	static bool doAll = false;

	static bool showAllParamsWindow = false;
	static bool showForceDetailWindow = false;

	static float sliceZ = 0.0f;

	static Model::Parameter modelParameter;
	static bool showModelParameterWindow = false;


	void Init()
	{
		currentDescriptionFile = std::filesystem::path("data\\C60\\C60 0.012 peak\\100x100x100_Ie0.012_Ucath44.2_RelTol0_Ni0_mbrc2_energies.asc");
		maxIndex = FileUtils::GetMaxIndex(currentDescriptionFile);
	}

	void CreateNewCurve()
	{
		curveList.emplace_back();
		currentCurveIndex = curveList.size() - 1;
	}

	void SetupCurve(std::filesystem::path folder, std::filesystem::path subfolder)
	{
		if (curveList.empty() || curveList.at(currentCurveIndex).IsNumerical())
		{
			CreateNewCurve();
		}
		Curve& currentCurve = curveList.at(currentCurveIndex);

		if (currentCurve.Empty())
		{
			currentCurve.SetFolder(folder);
			if (!subfolder.empty()) currentCurve.SetSubfolder(subfolder);
		}
		else if (currentCurve.GetFolder() != folder || (currentCurve.GetSubfolder() != subfolder && !subfolder.empty()))
		{
			CreateNewCurve();
			curveList.at(currentCurveIndex).SetFolder(folder);
			if (!subfolder.empty()) curveList.at(currentCurveIndex).SetSubfolder(subfolder);
		}
	}

	float GetSliceValue()
	{
		return sliceZ;
	}

	void ShowWindow()
	{
		if (ImGui::Begin("Cooling Force Window"))
		{
			ShowTabs();

			ImGui::SameLine();

			ImGui::BeginGroup();
			ShowSettings();
			ShowPlots();
			ImGui::EndGroup();

			ShowAllParametersWindow();
			modelParameter.ShowWindow(showModelParameterWindow);
			ShowForceDetailWindow();
		}
		ImGui::End();
	}

	void ShowSettings()
	{
		ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX | ImGuiChildFlags_ResizeY;
		if (ImGui::BeginChild("##ccSettings", ImVec2(0.0f, 0.0f), flags))
		{
			ImGui::SeparatorText("input things");
			if (ImGui::Button("select description file"))
			{
				currentDescriptionFile = FileUtils::SelectFile();
				maxIndex = FileUtils::GetMaxIndex(currentDescriptionFile);
			}
			ImGui::Text("file: %s", (currentDescriptionFile.parent_path().filename() / currentDescriptionFile.filename()).string().c_str());
			
			ImGui::Checkbox("All Parameters", &showAllParamsWindow);
			ImGui::SeparatorText("output things");

			if (!curveList.empty())
			{
				Curve& currentCurve = curveList.at(currentCurveIndex);
				if (!currentCurve.IsNumerical())
				{
					char buf[64] = "";
					strncpy_s(buf, currentCurve.GetSubfolder().string().c_str(), sizeof(buf) - 1);
					ImGui::SetNextItemWidth(150.0f);
					if (ImGui::InputText("curve subfolder", buf, IM_ARRAYSIZE(buf)))
					{
						currentCurve.SetSubfolder(std::filesystem::path(buf));
					}
					ImGui::SameLine();
				}
				if (ImGui::Button("save"))
				{
					currentCurve.Save();
				}
			}

			ImGui::Separator();

			ImGui::BeginDisabled(currentDescriptionFile.empty());
			if (ImGui::Button("Generate Cooling curve from File"))
			{
				std::filesystem::path folder = currentDescriptionFile.parent_path().parent_path().filename() /
					currentDescriptionFile.parent_path().filename();
				
				SetupCurve(folder);

				int start = startIndex;
				int end = endIndex;
				if (doAll)
				{
					start = 1;
					end = maxIndex;
				}
				for (int i = start; i <= end; i++)
				{
					Value newValue;
					if (method == 0) newValue.CalculateOriginal(currentDescriptionFile, i);
					if (method == 1) newValue.CalculateHalfIntegrated(currentDescriptionFile, i, interpolate);
					if (method == 2) newValue.CalculateFullIntegrated(currentDescriptionFile, i);
					if (method == 3) newValue.CalculateFullIntegratedBetter(currentDescriptionFile, i);
					Curve& curve = curveList.at(currentCurveIndex);
					curve.AddForceValue(std::move(newValue));
				}
			}
			
			ImGui::EndDisabled();
			ImGui::SameLine();
			if (ImGui::Button("Generate Cooling curve from numerical Model"))
			{
				CreateNewCurve();
				Curve& curve = curveList.at(currentCurveIndex);
				curve.IntegrateNumerically(modelParameter);
			}
			ImGui::SameLine();
			ImGui::Checkbox("parameter", &showModelParameterWindow);

			if (ImGui::RadioButton("original", method == 0))
			{
				method = 0;
			}
			ImGui::SameLine();
			if (ImGui::RadioButton("half integrated", method == 1))
			{
				method = 1;
			}
			ImGui::SameLine();
			if (ImGui::RadioButton("full integrated", method == 2))
			{
				method = 2;
			}
			ImGui::SameLine();
			if (ImGui::RadioButton("full integrated better", method == 3))
			{
				method = 3;
			}
			ImGui::SameLine();
			ImGui::Checkbox("interpolate pre calc force", &interpolate);
			ImGui::SameLine();
			Value::ShowParallelPrecalculationCheckbox();

			ImGui::SetNextItemWidth(80.0f);
			ImGui::BeginDisabled(doAll);
			ImGui::InputInt("start", &startIndex);
			ImGui::SameLine();
			ImGui::SetNextItemWidth(80.0f);
			ImGui::InputInt("end", &endIndex);
			ImGui::EndDisabled();
			ImGui::SameLine();
			ImGui::Checkbox("all", &doAll);
			ImGui::BeginDisabled(!doAll);
			ImGui::SameLine();
			ImGui::Text("(max Index: %d)", maxIndex);
			ImGui::EndDisabled();

		}
		ImGui::EndChild();
	}

	void ShowTabs()
	{
		ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX | ImGuiChildFlags_ResizeY;
		if (ImGui::BeginChild("cc tabs", ImVec2(400.0f, -1), flags))
		{
			ImGui::Text("Cooling force curves");
			ImGui::PushStyleColor(ImGuiCol_TabActive, ImVec4(0.7f, 0.15f, 0.15f, 1.0f));
			ImGui::PushStyleColor(ImGuiCol_TabHovered, ImVec4(0.8f, 0.2f, 0.2f, 1.0f));
			if (ImGui::BeginTabBar("##cc tab bar", ImGuiTabBarFlags_AutoSelectNewTabs))
			{
				if (ImGui::TabItemButton("+", ImGuiTabItemFlags_Trailing | ImGuiTabItemFlags_NoTooltip))
				{
					CreateNewCurve();
				}
				for (int curveIndex = 0; curveIndex < curveList.size(); curveIndex++)
				{
					bool open = true;
					std::string label = "curve " + std::to_string(curveIndex);
					if (ImGui::BeginTabItem(label.c_str(), &open))
					{
						currentCurveIndex = curveIndex;
						curveList.at(curveIndex).ShowContent();
						ImGui::SameLine();
						ImGui::Checkbox("show force details", &showForceDetailWindow);
						ImGui::EndTabItem();
					}
					if (!open)
					{
						curveList.erase(curveList.begin() + curveIndex);
						currentCurveIndex = std::min(currentCurveIndex, (int)curveList.size() - 1);
					}
				}
				ImGui::EndTabBar();
			}
			ImGui::PopStyleColor(2);

			ImGui::Separator();
			if (ImGui::Button("load MCMC Cool Curve"))
			{
				std::filesystem::path folder = FileUtils::SelectFolder(FileUtils::GetCoolingForceCurveFolder());
				if (!folder.empty())
				{
					CreateNewCurve();
					curveList.back().Load(folder);
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("load numerical Cool Curves"))
			{
				std::vector<std::filesystem::path> filenames = FileUtils::SelectFiles(FileUtils::GetNumericalCoolingForceCurveFolder());
				if (!filenames.empty())
				{
					for (auto& filename : filenames)
					{
						CreateNewCurve();
						curveList.back().Load(filename);
					}
				}
			}
			
		}
		ImGui::EndChild();
	}

	void ShowPlots()
	{
		if (ImPlot::BeginPlot("cooling force curves", ImVec2(-1, -1)))
		{
			ImPlot::SetupAxis(ImAxis_X1, "detuning velocity [m/s]");
			ImPlot::SetupAxis(ImAxis_Y1, "cooling force [eV/m]");
			ImPlot::SetupLegend(ImPlotLocation_NorthEast);

			int i = 0;
			for (const Curve& curve : curveList)
			{
					ImVec4 color = ImPlot::GetColormapColor(i);
					ImPlot::PushStyleColor(ImPlotCol_Line, color);
					ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, color);

					ImGui::PushID(i++);
					//curve.PlotForceX();
					//curve.PlotForceY();
					curve.PlotForceZ();
					ImGui::PopID();

					ImPlot::PopStyleColor(2);
			}
			if (showModelParameterWindow)
			{
				modelParameter.ShowVelocityLines();
			}
			ImPlot::EndPlot();
		}
	}

	void ShowAllParametersWindow()
	{
		if (!showAllParamsWindow)
		{
			return;
		}
		if (ImGui::Begin("all cooling force parameters", &showAllParamsWindow, ImGuiWindowFlags_NoDocking))
		{
			ImGuiChildFlags flags = ImGuiChildFlags_ResizeX | ImGuiChildFlags_Border;
			if (ImGui::BeginChild("labe", ImVec2(100, -1), flags))
			{
				ImGui::SeparatorText("Lab Energy parameter");
				LabEnergy::ShowParameterControls();
			}
			ImGui::EndChild();
			ImGui::SameLine();
			if (ImGui::BeginChild("ebeam", ImVec2(100, -1), flags))
			{
				ImGui::SeparatorText("Electron Beam parameter");
				ElectronBeam::ShowParameterControls();
			}
			ImGui::EndChild();
			ImGui::SameLine();
			if (ImGui::BeginChild("ibeam", ImVec2(100, -1), flags))
			{
				ImGui::SeparatorText("Ion Beam parameter");
				IonBeam::ShowParameterControls();
				ImGui::Separator();
				IonBeam::ShowCoolingForceParameterControls();
			}
			ImGui::EndChild();

		}
		ImGui::End();
	}

	void ShowForceDetailWindow()
	{
		if (!showForceDetailWindow)
		{
			return;
		}
		if (ImGui::Begin("cooling force details", &showForceDetailWindow, ImGuiWindowFlags_NoDocking))
		{
			if (currentCurveIndex >= 0)
			{
				Curve& curve = curveList.at(currentCurveIndex);
				if (ImGui::SliderFloat("z slice", &sliceZ, -0.7, 0.7))
				{
					curve.UpdateSlice(sliceZ);
				}
				curve.PlotDetails();
			}
		}
		ImGui::End();
	}
}
