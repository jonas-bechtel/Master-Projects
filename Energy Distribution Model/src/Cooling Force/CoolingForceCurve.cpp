#include "pch.h"
#include "CoolingForceCurve.h"

CoolingForceCurve::CoolingForceCurve()
{
}

void CoolingForceCurve::AddForceValue(CoolingForceValue&& value)
{
	forceX.push_back(value.forceXValue);
	forceY.push_back(value.forceYValue);
	forceZ.push_back(value.forceZValue);

	// will call move Constructor
	values.emplace_back(std::move(value));
	detuningVelocites.push_back(value.eBeamParameter.detuningVelocity);
}

void CoolingForceCurve::RemoveForceValue(int index)
{
	forceX.erase(forceX.begin() + index);
	forceY.erase(forceY.begin() + index);
	forceZ.erase(forceZ.begin() + index);

	values.erase(values.begin() + index);
}

void CoolingForceCurve::SetFolder(std::filesystem::path path)
{
	folder = path;
}

void CoolingForceCurve::SetSubfolder(std::filesystem::path path)
{
	subFolder = path;
}

std::filesystem::path CoolingForceCurve::GetFolder() const
{
	return folder;
}

std::filesystem::path CoolingForceCurve::GetSubfolder() const
{
	return subFolder;
}

std::string CoolingForceCurve::GetLabel() const
{
	return (folder / subFolder).string();
}

bool CoolingForceCurve::Empty() const
{
	return values.empty();
}

void CoolingForceCurve::ShowList() 
{
	ImGui::PushID(this);

	float sizeY = ImGui::GetContentRegionAvail().y - 100.0f;
	if (ImGui::BeginListBox("listbox", ImVec2(-1, sizeY)))
	{
		for (int i = 0; i < values.size(); i++)
		{
			ImGui::PushID(i);
			const CoolingForceValue& value = values.at(i);
			if (value.ShowListItem(selectedIndex == i))
			{
				selectedIndex = i;
			}

			ImGui::SameLine();
			if (ImGui::SmallButton("x"))
			{
				RemoveForceValue(i);
			}
			ImGui::PopID();
		}

		ImGui::EndListBox();
	}
	ImGui::PopID();

	ImGui::Text("cooling force curve: %s", GetLabel().c_str());

	if (ImGui::Button("save cooling curve"))
	{
		//Save();
	}
}

void CoolingForceCurve::Plot() const
{
	ImPlot::PlotScatter(GetLabel().c_str(), detuningVelocites.data(), forceZ.data(), detuningVelocites.size());
}
