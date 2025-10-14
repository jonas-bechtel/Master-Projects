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
		detuningVelocitesMovedScaled.reserve(params.numberPoints);
		force.reserve(params.numberPoints);
		forceMovedScaled.reserve(params.numberPoints);
		
		simpleModel = true;
		modelParams = params;

		for (int i = 0; i < params.numberPoints; i++)
		{
			double forceValue = 0.0;
			if (Value::calculateTransverseForce)
			{
				modelParams.relativeVelocity.SetX(params.relativeVelocityRange[0] + i * step);
				forceValue = Model::Force(modelParams, 1);

				detuningVelocites.push_back(modelParams.relativeVelocity.x());
				detuningVelocitesMovedScaled.push_back(modelParams.relativeVelocity.x());
			}
			else
			{
				modelParams.relativeVelocity.SetZ(params.relativeVelocityRange[0] + i * step);
				forceValue = Model::Force(modelParams);

				detuningVelocites.push_back(modelParams.relativeVelocity.z());
				detuningVelocitesMovedScaled.push_back(modelParams.relativeVelocity.z());
			}
			
			force.push_back(forceValue);
			forceMovedScaled.push_back(forceValue);
		}
	}

	void Curve::AddForceValue(Value&& value)
	{
		force.push_back(value.forceValue);
		forceMovedScaled.push_back(value.forceValue);

		// will call move Constructor
		detuningVelocites.push_back(value.eBeamParameter.detuningVelocity);
		detuningVelocitesMovedScaled.push_back(value.eBeamParameter.detuningVelocity);
		values.emplace_back(std::move(value));

		if (values.size() == 1)
		{
			selectedIndex = 0;
			SelectedItemChanged();
		}
	}

	void Curve::RemoveForceValue(int index)
	{
		force.erase(force.begin() + index);
		forceMovedScaled.erase(forceMovedScaled.begin() + index);

		detuningVelocites.erase(detuningVelocites.begin() + index);
		detuningVelocitesMovedScaled.erase(detuningVelocitesMovedScaled.begin() + index);
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
		if (simpleModel)
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

	bool Curve::IsSimpleModel() const
	{
		return simpleModel;
	}

	bool Curve::IsMeasured() const
	{
		return measured;
	}

	void Curve::ShowContent()
	{
		ImGui::Text("label: %s", GetLabel().c_str());
		if (simpleModel)
		{
			ImGui::Separator();
			modelParams.ShowValues();

			//if (ImGui::BeginDragDropSource(ImGuiDragDropFlags_None))
			//{
			//	ImGui::SetDragDropPayload("Cooling Force Curve", this, sizeof(Curve));
			//	ImGui::Text("drag to set params");
			//	ImGui::EndDragDropSource();
			//}
			//ImGui::Separator();
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
		for (int i = 0; i < force.size(); i++)
		{
			forceMovedScaled[i] = force[i] * scale + moveForce;
			if (measured)
				forceErrorMovedScaled[i] = forceError[i] * scale;
			detuningVelocitesMovedScaled[i] = detuningVelocites[i] + moveVelocity;
		}
	}

	void Curve::PlotForce() const
	{
		if (simpleModel)
		{
			ImPlot::PlotLine(GetLabel().c_str(), detuningVelocitesMovedScaled.data(), forceMovedScaled.data(), forceMovedScaled.size());
		}
		else if (measured)
		{
			ImPlot::PlotScatter(GetLabel().c_str(), detuningVelocitesMovedScaled.data(), forceMovedScaled.data(), forceMovedScaled.size());
			ImPlot::PlotErrorBars(GetLabel().c_str(), detuningVelocitesMovedScaled.data(), forceMovedScaled.data(), forceErrorMovedScaled.data(), forceMovedScaled.size());
		}
		else
		{
			ImPlot::PlotScatter(GetLabel().c_str(), detuningVelocitesMovedScaled.data(), forceMovedScaled.data(), forceMovedScaled.size());
		}
	}

	void Curve::PlotDetails() const
	{
		if (selectedIndex < 0) return;
		const Value& value = values.at(selectedIndex);

		//ImPlot::BeginPlot("##details");
		value.PlotForceSlice();
		//ImPlot::EndPlot();

	}

	void Curve::UpdateSlice(float zValue)
	{
		if (selectedIndex < 0) return;

		values.at(selectedIndex).UpdateSlice(zValue);
	}

	void Curve::Save() const
	{
		if (simpleModel)
		{
			std::filesystem::path outfolder = FileUtils::GetSimpleModelCoolingForceCurveFolder();

			if (!std::filesystem::exists(outfolder))
			{
				std::filesystem::create_directories(outfolder);
			}

			std::filesystem::path file = outfolder / (modelParams.String() + ".curve");
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
				outfile << detuningVelocites.at(i) << "\t" << force.at(i) << "\n";
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

		// also save the important data 
		std::string filename = folder.parent_path().filename().string() + "_" + folder.filename().string() + "_" + subFolder.string() + ".curve";
		std::filesystem::path file = outfolder / filename;
		std::ofstream outfile(file);

		if (!outfile.is_open())
		{
			std::cerr << "Error opening file" << std::endl;
			return;
		}
		outfile << "# detuning velocity [m/s]\tcooling force [eV/m]\n";
		for (int i = 0; i < detuningVelocites.size(); i++)
		{
			outfile << detuningVelocites.at(i) << "\t" << force.at(i) << "\n";
		}

		outfile.close();
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
				if (entry.is_regular_file() && (entry.path().extension() == ".asc" || entry.path().extension() == ".root"))
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
				detuningVelocitesMovedScaled.push_back(std::stod(tokens[0]));
				force.push_back(std::stod(tokens[1]));
				forceMovedScaled.push_back(std::stod(tokens[1]));
			}
			infile.close();

			modelParams.FromString(input.filename().replace_extension().string());
			simpleModel = true;
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

		std::string line;
		while (std::getline(infile, line))
		{
			std::vector<std::string> tokens = FileUtils::SplitLine(line, "\t");

			detuningVelocites.push_back(std::stod(tokens[1]));
			detuningVelocitesMovedScaled.push_back(std::stod(tokens[1]));
			force.push_back(std::stod(tokens[2]));
			forceMovedScaled.push_back(std::stod(tokens[2]));
			forceError.push_back(std::stod(tokens[3]));
			forceErrorMovedScaled.push_back(std::stod(tokens[3]));
		}
		infile.close();
		
		subFolder = input.filename().replace_extension().string();

		measured = true;
	}
}
