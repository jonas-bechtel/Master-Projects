#include "pch.h"
#include "CoolingForceWindow.h"

#include "CoolingForceCurve.h"
#include "FileUtils.h"

namespace CoolingForceWindow
{
	// main storge data
	static std::vector<CoolingForceCurve> curveList;
	static int currentCurveIndex = -1;

	// currently loaded description file
	static std::filesystem::path currentDescriptionFile = std::filesystem::path("data\\C60\\dataset1\\100x100x100_Ie0.95_Ucath44.2_RelTol0_mbrc1_energies.asc");
	static int maxIndex = 0;

	// start/end index in description file to generate cooling force for
	static int startIndex = 1;
	static int endIndex = 1;
	static bool doAll = false;

	static bool showAllParamsWindow = false;


	void Init()
	{
	}

	void CreateNewCurve()
	{
		curveList.emplace_back();
		currentCurveIndex = curveList.size() - 1;
	}

	void SetupCurve(std::filesystem::path folder, std::filesystem::path subfolder)
	{
		if (curveList.empty())
		{
			CreateNewCurve();
		}
		CoolingForceCurve& currentCurve = curveList.at(currentCurveIndex);
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
				CoolingForceCurve& currentCurve = curveList.at(currentCurveIndex);
				char buf[64] = "";
				strncpy_s(buf, currentCurve.GetSubfolder().string().c_str(), sizeof(buf) - 1);
				ImGui::SetNextItemWidth(150.0f);
				if (ImGui::InputText("set subfolder", buf, IM_ARRAYSIZE(buf)))
				{
					currentCurve.SetSubfolder(std::filesystem::path(buf));
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
					start = 0;
					end = maxIndex;
				}
				for (int i = start; i <= end; i++)
				{
					CoolingForceValue newValue;
					newValue.Calculate(currentDescriptionFile, i);
					CoolingForceCurve& curve = curveList.at(currentCurveIndex);
					curve.AddForceValue(std::move(newValue));
				}
			}
			ImGui::EndDisabled();

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
						curveList.at(curveIndex).ShowList();
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

		}
		ImGui::EndChild();
	}

	void ShowPlots()
	{
		if (ImPlot::BeginPlot("cooling force curves", ImVec2(-1, -1)))
		{
			ImPlot::SetupAxis(ImAxis_X1, "detuning velocity [m/s]");
			ImPlot::SetupAxis(ImAxis_Y1, "cooling force [eV/m]", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit);
			ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);

			int i = 0;
			for (const CoolingForceCurve& curve : curveList)
			{
					ImGui::PushID(i++);
					curve.Plot();
					ImGui::PopID();
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
				IonBeam::ShowCoolingForceParameterControls();
				IonBeam::ShowParameterControls();
			}
			ImGui::EndChild();

		}
		ImGui::End();
	}
}
