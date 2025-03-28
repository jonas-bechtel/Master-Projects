#include "pch.h"
#include "CoolingForceCurve.h"
#include "Constants.h"
#include "FileUtils.h"

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
	forceZscaled.reserve(params.numberPoints);

	TF2 func("func", CoolingForceModel::NumericalIntegrand, 0.0, 5.0 * deltaTrans, -5.0 * deltaLong, 5.0 * deltaLong, 5, 2);
	for (int i = 0; i < params.numberPoints; i++)
	{
		params.relativeVelocity = start + i * step;
		func.SetParameters((double*)&params);

		// Now integrate the function over the specified range
		double result = func.Integral(0.0, 5.0 * deltaTrans, -5.0 * deltaLong, 5.0 * deltaLong, 1e-8);

		detuningVelocites.push_back(params.relativeVelocity);
		forceZ.push_back(result);
		forceZscaled.push_back(result);
	}
	numerical = true;
	numericalParams = params;
}

void CoolingForceCurve::AddForceValue(CoolingForceValue&& value)
{
	forceX.push_back(value.forceXValue);
	forceY.push_back(value.forceYValue);
	forceZ.push_back(value.forceZValue);
	forceZscaled.push_back(value.forceZValue);

	// will call move Constructor
	detuningVelocites.push_back(value.eBeamParameter.detuningVelocity);
	values.emplace_back(std::move(value));
}

void CoolingForceCurve::RemoveForceValue(int index)
{
	forceX.erase(forceX.begin() + index);
	forceY.erase(forceY.begin() + index);
	forceZ.erase(forceZ.begin() + index);
	forceZscaled.erase(forceZscaled.begin() + index);

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

	float sizeY = ImGui::GetContentRegionAvail().y - 200.0f;
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

	ImGui::SetNextItemWidth(200.0f);
	bool changed = ImGui::SliderFloat("scale", &scale, 0, 10);
	if (ImGui::IsItemHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Right))
	{
		scale = 1.0f;
		changed = true;
	}
	if (changed)
	{
		for (int i = 0; i < forceZ.size(); i++)
		{
			forceZscaled[i] = forceZ[i] * scale;
		}
	}
	
	ImGui::Text("cooling force curve: %s", GetLabel().c_str());
}

void CoolingForceCurve::PlotForceX() const
{
	ImPlot::PlotScatter(GetLabel().c_str(), detuningVelocites.data(), forceX.data(), forceX.size());
}

void CoolingForceCurve::PlotForceY() const
{
	ImPlot::PlotScatter(GetLabel().c_str(), detuningVelocites.data(), forceY.data(), forceY.size());
}

void CoolingForceCurve::PlotForceZ() const
{
	if(numerical)
		ImPlot::PlotLine(GetLabel().c_str(), detuningVelocites.data(), forceZscaled.data(), forceZscaled.size());
	else
		ImPlot::PlotScatter(GetLabel().c_str(), detuningVelocites.data(), forceZscaled.data(), forceZscaled.size());
}

void CoolingForceCurve::Save() const
{
	if (numerical)
	{
		std::filesystem::path outfolder = FileUtils::GetNumericalCoolingForceCurveFolder();

		if (!std::filesystem::exists(outfolder))
		{
			std::filesystem::create_directories(outfolder);
		}

		std::filesystem::path file = outfolder / (numericalParams.String() + ".asc");
		std::ofstream outfile(file);

		if (!outfile.is_open())
		{
			std::cerr << "Error opening file" << std::endl;
			return;
		}
		outfile << "# " << numericalParams.String() << "\n";
		outfile << "# detuning velocity [m/s]\tcooling force [eV/m]\n";
		for (int i = 0; i < detuningVelocites.size(); i++)
		{
			outfile << detuningVelocites.at(i) << "\t" << forceZ.at(i) <<  "\n";
		}

		outfile.close();
		return;
	}

	std::filesystem::path outfolder = FileUtils::GetCoolingForceCurveFolder() / folder / subFolder;

	if (!std::filesystem::exists(outfolder))
	{
		std::filesystem::create_directories(outfolder);
	}

	for (const CoolingForceValue& value : values)
	{
		value.Save(outfolder);
	}
}

void CoolingForceCurve::Load(const std::filesystem::path& input)
{
	if (!std::filesystem::exists(input))
	{
		std::cerr << "Invalid input path!" << std::endl;
		return;
	}

	// loading mcmc generated cooling curve
	if (std::filesystem::is_directory(input))
	{
		for (const auto& entry : std::filesystem::directory_iterator(input))
		{
			if (entry.is_regular_file() && entry.path().extension() == ".asc")
			{
				std::filesystem::path file = entry.path();
				CoolingForceValue value;
				value.Load(file);

				AddForceValue(std::move(value));
			}
		}
		folder = input.parent_path().parent_path().filename() / input.parent_path().filename();
		subFolder = input.filename();
	}
	// loading numerically integrated cooling curves
	else if (std::filesystem::is_regular_file(input))
	{
		// load the .asc file with the histogram data
		std::ifstream infile(input);

		// Check if the file was successfully opened
		if (!infile.is_open())
		{
			std::cerr << "Error: Could not open the file " << input << std::endl;
			return;
		}

		// skip the header
		FileUtils::GetHeaderFromFile(infile);

		std::string line;
		while (std::getline(infile, line))
		{
			std::vector<std::string> tokens = FileUtils::SplitLine(line, "\t");

			detuningVelocites.push_back(std::stod(tokens[0]));
			forceZ.push_back(std::stod(tokens[1]));
			forceZscaled.push_back(std::stod(tokens[1]));
		}
		infile.close();

		numericalParams.FromString(input.filename().replace_extension().string());
		numerical = true;
	}
}
