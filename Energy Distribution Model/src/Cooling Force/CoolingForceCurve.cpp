#include "pch.h"
#include "CoolingForceCurve.h"
#include "Constants.h"

CoolingForceCurve::CoolingForceCurve()
{
}

void CoolingForceCurve::IntegrateNumerically(NumericalIntegrationParameter& params)
{
	double deltaTrans = sqrt(2 * params.kT_trans * TMath::Qe() / PhysicalConstants::electronMass);
	double deltaLong = sqrt(params.kT_long * TMath::Qe() / PhysicalConstants::electronMass);

	double start = params.relativeVelocityRange[0];
	double step = (params.relativeVelocityRange[1] - params.relativeVelocityRange[0]) / (params.numberPoints - 1);

	detuningVelocites.reserve(params.numberPoints);
	forceZ.reserve(params.numberPoints);

	TF2 func("func", CoolingForceModel::NumericalIntegrand, 0.0, 5.0 * deltaTrans, -5.0 * deltaLong, 5.0 * deltaLong, 5, 2);
	for (int i = 0; i < params.numberPoints; i++)
	{
		params.relativeVelocity = start + i * step;
		func.SetParameters((double*)&params);

		// Now integrate the function over the specified range
		double result = func.Integral(0.0, 5.0 * deltaTrans, -5.0 * deltaLong, 5.0 * deltaLong);

		detuningVelocites.push_back(params.relativeVelocity);
		forceZ.push_back(result);
	}
	numerical = true;
}

void CoolingForceCurve::AddForceValue(CoolingForceValue&& value)
{
	forceX.push_back(value.forceXValue);
	forceY.push_back(value.forceYValue);
	forceZ.push_back(value.forceZValue);

	// will call move Constructor
	detuningVelocites.push_back(value.eBeamParameter.detuningVelocity);
	values.emplace_back(std::move(value));
}

void CoolingForceCurve::RemoveForceValue(int index)
{
	forceX.erase(forceX.begin() + index);
	forceY.erase(forceY.begin() + index);
	forceZ.erase(forceZ.begin() + index);

	detuningVelocites.erase(detuningVelocites.begin() + index);
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
	if (ImGui::BeginListBox("##cc listbox", ImVec2(-1, sizeY)))
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

void CoolingForceCurve::PlotForceZ() const
{
	if(numerical)
		ImPlot::PlotLine(GetLabel().c_str(), detuningVelocites.data(), forceZ.data(), detuningVelocites.size());
	else
		ImPlot::PlotScatter(GetLabel().c_str(), detuningVelocites.data(), forceZ.data(), detuningVelocites.size());
}
