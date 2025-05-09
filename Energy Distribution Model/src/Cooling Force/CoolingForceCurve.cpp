#include "pch.h"
#include "CoolingForceCurve.h"
#include "CoolingForceWindow.h"
#include "Constants.h"
#include "FileUtils.h"

namespace CoolingForce
{
	Curve::Curve()
	{
	}

	void Curve::IntegrateNumerically(const Model::Parameter& params)
	{
		double step = (params.relativeVelocityRange[1] - params.relativeVelocityRange[0]) / (params.numberPoints - 1);

		detuningVelocites.reserve(params.numberPoints);
		detuningVelocitesMoved.reserve(params.numberPoints);
		forceZ.reserve(params.numberPoints);
		forceZmovedScaled.reserve(params.numberPoints);
		
		numerical = true;
		modelParams = params;

		for (int i = 0; i < params.numberPoints; i++)
		{
			modelParams.relativeVelocity.SetZ(params.relativeVelocityRange[0] + i * step);
			detuningVelocites.push_back(modelParams.relativeVelocity.z());
			detuningVelocitesMoved.push_back(modelParams.relativeVelocity.z());

			double force = Model::ForceZ(modelParams);
			
			forceZ.push_back(force);
			forceZmovedScaled.push_back(force);
		}

		
	}

	void Curve::AddForceValue(Value&& value)
	{
		forceX.push_back(value.forceXValue);
		forceY.push_back(value.forceYValue);
		forceZ.push_back(value.forceZValue);
		forceZmovedScaled.push_back(value.forceZValue);

		// will call move Constructor
		detuningVelocites.push_back(value.eBeamParameter.detuningVelocity);
		detuningVelocitesMoved.push_back(value.eBeamParameter.detuningVelocity);
		values.emplace_back(std::move(value));

		if (values.size() == 1)
		{
			selectedIndex = 0;
			SelectedItemChanged();
		}
	}

	void Curve::RemoveForceValue(int index)
	{
		forceX.erase(forceX.begin() + index);
		forceY.erase(forceY.begin() + index);
		forceZ.erase(forceZ.begin() + index);
		forceZmovedScaled.erase(forceZmovedScaled.begin() + index);

		detuningVelocites.erase(detuningVelocites.begin() + index);
		detuningVelocitesMoved.erase(detuningVelocitesMoved.begin() + index);
		values.erase(values.begin() + index);

		selectedIndex = std::min(selectedIndex, (int)values.size() - 1);
		if (selectedIndex >= 0)
		{
			SelectedItemChanged();
		}
	}

	void Curve::SetFolder(std::filesystem::path path)
	{
		folder = path;
	}

	void Curve::SetSubfolder(std::filesystem::path path)
	{
		subFolder = path;
	}

	std::filesystem::path Curve::GetFolder() const
	{
		return folder;
	}

	std::filesystem::path Curve::GetSubfolder() const
	{
		return subFolder;
	}

	std::string Curve::GetLabel() const
	{
		if (numerical)
			return modelParams.String();
		else if (measured)
			return subFolder.string();
		else
			return (folder / subFolder).string();
	}

	bool Curve::Empty() const
	{
		return values.empty();
	}

	bool Curve::IsNumerical() const
	{
		return numerical;
	}

	bool Curve::IsMeasured() const
	{
		return measured;
	}

	void Curve::ShowContent()
	{
		ImGui::Text("label: %s", GetLabel().c_str());
		if (numerical)
		{
			ImGui::Separator();
			modelParams.ShowValues();
			ImGui::Separator();
		}
		else if (measured)
		{
			//measuredParams.ShowValues();
		}
		else
		{
			ImGui::PushID(this);

			float sizeY = ImGui::GetContentRegionAvail().y - 200.0f;
			if (ImGui::BeginListBox("##cc listbox", ImVec2(-1, sizeY)))
			{
				for (int i = 0; i < values.size(); i++)
				{
					ImGui::PushID(i);
					const Value& value = values.at(i);
					if (value.ShowListItem(selectedIndex == i))
					{
						selectedIndex = i;
						SelectedItemChanged();
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
		}

		ImGui::PushItemWidth(200.0f);
		bool changed = ImGui::SliderFloat("scale", &scale, 0.0f, 10.0f);
		if (ImGui::IsItemHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Right))
		{
			scale = 1.0f;
			changed = true;
		}
		changed |= ImGui::SliderFloat("move vertically [eV/m]", &moveForce, -0.05f, 0.05f);
		if (ImGui::IsItemHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Right))
		{
			moveForce = 0.0f;
			changed = true;
		}
		changed |= ImGui::SliderFloat("move horizontally [m/s]", &moveVelocity, -1e4f, 1e4f);
		if (ImGui::IsItemHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Right))
		{
			moveVelocity = 0.0f;
			changed = true;
		}
		if (changed)
		{
			UpdateMovedScaledValues();
		}
		ImGui::PopItemWidth();
	}

	void Curve::SelectedItemChanged()
	{
		UpdateSlice(GetSliceValue());
	}

	void Curve::UpdateMovedScaledValues()
	{
		for (int i = 0; i < forceZ.size(); i++)
		{
			forceZmovedScaled[i] = forceZ[i] * scale + moveForce;
			if (measured)
				forceZErrorScaled[i] = forceZError[i] * scale;
			detuningVelocitesMoved[i] = detuningVelocites[i] + moveVelocity;
		}
	}

	void Curve::PlotForceX() const
	{
		ImPlot::PlotScatter(GetLabel().c_str(), detuningVelocites.data(), forceX.data(), forceX.size());
	}

	void Curve::PlotForceY() const
	{
		ImPlot::PlotScatter(GetLabel().c_str(), detuningVelocites.data(), forceY.data(), forceY.size());
	}

	void Curve::PlotForceZ() const
	{
		if (numerical)
		{
			ImPlot::PlotLine(GetLabel().c_str(), detuningVelocitesMoved.data(), forceZmovedScaled.data(), forceZmovedScaled.size());
		}
		else if (measured)
		{
			ImPlot::PlotScatter(GetLabel().c_str(), detuningVelocitesMoved.data(), forceZmovedScaled.data(), forceZmovedScaled.size());
			ImPlot::PlotErrorBars(GetLabel().c_str(), detuningVelocitesMoved.data(), forceZmovedScaled.data(), forceZErrorScaled.data(), forceZmovedScaled.size());
		}
		else
		{
			ImPlot::PlotScatter(GetLabel().c_str(), detuningVelocitesMoved.data(), forceZmovedScaled.data(), forceZmovedScaled.size());
		}
	}

	void Curve::PlotDetails() const
	{
		if (selectedIndex < 0) return;
		const Value& value = values.at(selectedIndex);

		//ImPlot::BeginPlot("##details");
		value.PlotPreForceSlize();
		//ImPlot::EndPlot();

	}

	void Curve::UpdateSlice(float zValue)
	{
		if (selectedIndex < 0) return;

		values.at(selectedIndex).UpdateSlice(zValue);
	}

	void Curve::Save() const
	{
		if (numerical)
		{
			std::filesystem::path outfolder = FileUtils::GetNumericalCoolingForceCurveFolder();

			if (!std::filesystem::exists(outfolder))
			{
				std::filesystem::create_directories(outfolder);
			}

			std::filesystem::path file = outfolder / (modelParams.String() + ".asc");
			std::ofstream outfile(file);

			if (!outfile.is_open())
			{
				std::cerr << "Error opening file" << std::endl;
				return;
			}
			outfile << "# " << modelParams.String() << "\n";
			outfile << "# detuning velocity [m/s]\tcooling force [eV/m]\n";
			for (int i = 0; i < detuningVelocites.size(); i++)
			{
				outfile << detuningVelocites.at(i) << "\t" << forceZ.at(i) << "\n";
			}

			outfile.close();
			return;
		}

		std::filesystem::path outfolder = FileUtils::GetCoolingForceCurveFolder() / folder / subFolder;

		if (!std::filesystem::exists(outfolder))
		{
			std::filesystem::create_directories(outfolder);
		}

		for (const Value& value : values)
		{
			value.Save(outfolder);
		}
	}

	void Curve::Load(const std::filesystem::path& input)
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
					Value value;
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
				detuningVelocitesMoved.push_back(std::stod(tokens[0]));
				forceZ.push_back(std::stod(tokens[1]));
				forceZmovedScaled.push_back(std::stod(tokens[1]));
			}
			infile.close();

			modelParams.FromString(input.filename().replace_extension().string());
			numerical = true;
		}
	}

	void Curve::LoadMeasured(const std::filesystem::path& input)
	{
		if (!std::filesystem::exists(input) || !std::filesystem::is_regular_file(input))
		{
			std::cerr << "Invalid input path!" << std::endl;
			return;
		}

		std::ifstream infile(input);

		// Check if the file was successfully opened
		if (!infile.is_open())
		{
			std::cerr << "Error: Could not open the file " << input << std::endl;
			return;
		}

		std::string header = FileUtils::GetHeaderFromFile(infile);
		//measuredParams.FromString(header);

		std::string line;
		while (std::getline(infile, line))
		{
			std::vector<std::string> tokens = FileUtils::SplitLine(line, "\t");

			detuningVelocites.push_back(std::stod(tokens[1]));
			detuningVelocitesMoved.push_back(std::stod(tokens[1]));
			forceZ.push_back(std::stod(tokens[2]));
			forceZmovedScaled.push_back(std::stod(tokens[2]));
			forceZError.push_back(std::stod(tokens[3]));
			forceZErrorScaled.push_back(std::stod(tokens[3]));
		}
		infile.close();
		
		subFolder = input.filename().replace_extension().string();

		measured = true;
	}

	/*void MeasuredCurveParameter::FromString(std::string& input)
	{
		std::istringstream stream(input);
		std::string token;

		while (std::getline(stream, token, '_'))
		{
			if (token.find("Q=") == 0)
			{
				ionCharge = std::stoi(token.substr(2));
			}
			else if (token.find("ne=") == 0)
			{
				electronDensity = std::stod(token.substr(3));
			}
			else if (token.find("Tlong=") == 0)
			{
				kT_long = std::stod(token.substr(6));
			}
			else if (token.find("Ttrans=") == 0)
			{
				kT_trans = std::stod(token.substr(7));
			}
			else if (token.find("B=") == 0)
			{
				magneticField = std::stod(token.substr(2));
			}
			else if (token.find("Veff=") == 0)
			{
				effectiveVelocity = std::stod(token.substr(5));
			}
		}
	}*/
}
