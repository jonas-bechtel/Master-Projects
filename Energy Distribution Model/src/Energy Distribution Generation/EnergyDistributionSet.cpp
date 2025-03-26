#include "pch.h"
#include "EnergyDistributionSet.h"
#include "FileUtils.h"

void SetInformation::AddDistributionValues(const EnergyDistribution& dist)
{
	indeces.push_back(dist.index);
	centerLabEnergy.push_back(dist.labEnergiesParameter.centerLabEnergy.get());
	detuningEnergy.push_back(dist.eBeamParameter.detuningEnergy.get());

	fitDetuningEnergy.push_back(dist.outputParameter.fitDetuningEnergy.get());
	fitLongitudinalTemperature.push_back(dist.outputParameter.fitLongitudinalTemperature.get());
	fitTransverseTemperature.push_back(dist.outputParameter.fitTransverseTemperature.get());
	fitFWHM.push_back(dist.outputParameter.fitFWHM.get());
	fitScalingFactor.push_back(dist.outputParameter.fitScalingFactor.get());
	FWHM.push_back(dist.outputParameter.FWHM.get());
}

void SetInformation::RemoveDistributionValues(int index)
{
	indeces.erase(indeces.begin() + index);
	centerLabEnergy.erase(centerLabEnergy.begin() + index);
	detuningEnergy.erase(detuningEnergy.begin() + index);

	fitDetuningEnergy.erase(fitDetuningEnergy.begin() + index);
	fitLongitudinalTemperature.erase(fitLongitudinalTemperature.begin() + index);
	fitTransverseTemperature.erase(fitTransverseTemperature.begin() + index);
	fitFWHM.erase(fitFWHM.begin() + index);
	fitScalingFactor.erase(fitScalingFactor.begin() + index);
	FWHM.erase(FWHM.begin() + index);
}

void SetInformation::PlotFitEd(std::string setLabel)
{
	if (!plot) return;
	ImPlot::PlotLine(setLabel.c_str(), detuningEnergy.data(), fitDetuningEnergy.data(), detuningEnergy.size());
}

void SetInformation::PlotFitlongkT(std::string setLabel)
{
	if (!plot) return;
	ImPlot::PlotLine(setLabel.c_str(), detuningEnergy.data(), fitLongitudinalTemperature.data(), detuningEnergy.size());
}

void SetInformation::PlotFitTranskT(std::string setLabel)
{
	if (!plot) return;
	ImPlot::PlotLine(setLabel.c_str(), detuningEnergy.data(), fitTransverseTemperature.data(), detuningEnergy.size());
}

void SetInformation::PlotFitScalingFactor(std::string setLabel)
{
	if (!plot) return;
	ImPlot::PlotLine(setLabel.c_str(), detuningEnergy.data(), fitScalingFactor.data(), detuningEnergy.size());
}

void SetInformation::PlotFitFWHM(std::string setLabel)
{
	if (!plot) return;
	ImPlot::PlotLine(setLabel.c_str(), detuningEnergy.data(), fitFWHM.data(), detuningEnergy.size());
}

void SetInformation::PlotFWHM(std::string setLabel)
{
	if (!plot) return;
	ImPlot::PlotLine(setLabel.c_str(), detuningEnergy.data(), FWHM.data(), detuningEnergy.size());
}

void SetInformation::Save(std::filesystem::path folder)
{
	// set the output filepath
	std::filesystem::path file = FileUtils::GetEnergyDistSetFolder() / folder / "Info.txt";

	// Create the directories if they don't exist
	if (!std::filesystem::exists(file.parent_path()))
	{
		std::filesystem::create_directories(file.parent_path());
	}

	std::ofstream outfile(file);

	if (!outfile.is_open())
	{
		std::cerr << "Error opening file" << std::endl;
		return;
	}

	outfile << "# index\tE_lab\tE_d\tfit E_d\tfit kT_long\tfit kT_trans\tfit scaling factor\tfit FWHM\tFWHM\n";

	for (int i = 0; i < detuningEnergy.size(); i++)
	{
		outfile << indeces[i] << "\t" << centerLabEnergy[i] << "\t" << detuningEnergy[i] << "\t" << fitDetuningEnergy[i] << "\t" <<
			fitLongitudinalTemperature[i] << "\t" << fitTransverseTemperature[i] << "\t" << fitScalingFactor[i] << "\t" <<
			fitFWHM[i] << "\t" << FWHM[i] << "\n";
	}

	outfile.close();
}

void EnergyDistributionSet::AddDistribution(EnergyDistribution&& distribution)
{
	info.AddDistributionValues(distribution);

	// will call move Constructor
	distributions.emplace_back(std::move(distribution));
	EnergyDistribution& justMoved = distributions.back();

	EdToDistMap[justMoved.GetDetuningEnergy()] = &justMoved;
}

void EnergyDistributionSet::RemoveDistribution(int index)
{
	//EnergyDistribution& distToRemove = distributions.at(index);
	info.RemoveDistributionValues(index);


	// edist needs to be removed from map
	for (auto it = EdToDistMap.begin(); it != EdToDistMap.end(); it++)
	{
		if (it->second == &distributions.at(index))
		{
			EdToDistMap.erase(it);
			break;
		}
	}
	distributions.erase(distributions.begin() + index);
}

EnergyDistribution* EnergyDistributionSet::FindByEd(double detuningEnergy)
{
	if (EdToDistMap.find(detuningEnergy) == EdToDistMap.end())
	{
		std::cout << "no energy distribution with E_d = " << detuningEnergy << " was found\n";
		return nullptr;
	}
	return EdToDistMap.at(detuningEnergy);
}

void EnergyDistributionSet::SetFolder(std::filesystem::path path)
{
	folder = path;
}

void EnergyDistributionSet::SetSubfolder(std::filesystem::path folder)
{
	subFolder = folder;
}

void EnergyDistributionSet::SetAllPlotted(bool plotted)
{
	for (EnergyDistribution& eDist : distributions)
	{
		eDist.SetPlot(plotted);
	}
}

void EnergyDistributionSet::SetAllShowNormalised(bool showNormalised)
{
	for (EnergyDistribution& eDist : distributions)
	{
		eDist.SetNormalised(showNormalised);
	}
}

std::filesystem::path EnergyDistributionSet::GetFolder()
{
	return folder;
}

std::filesystem::path EnergyDistributionSet::GetSubfolder()
{
	return subFolder;
}

const std::vector<EnergyDistribution>& EnergyDistributionSet::GetDistributions() const
{
	return distributions;
}

SetInformation& EnergyDistributionSet::GetInfo()
{
	return info;
}

void EnergyDistributionSet::CalculatePsisFromBinning(TH1D* crossSection)
{
	for (EnergyDistribution& eDist : distributions)
	{
		eDist.CalculatePsisFromBinning(crossSection);
	}
}

void EnergyDistributionSet::ShowList()
{
	//ImGui::PushID(this);
	std::cout << ImGui::GetItemID() << std::endl;
	float sizeY = ImGui::GetContentRegionAvail().y - 100.0f;
	if (ImGui::BeginListBox("edist listbox", ImVec2(-1, sizeY)))
	{
		for (int i = 0; i < distributions.size(); i++)
		{
			ImGui::PushID(i);
			EnergyDistribution& eDist = distributions[i];
			eDist.ShowListItem();
			
			ImGui::SameLine();
			if (ImGui::SmallButton("x"))
			{
				RemoveDistribution(i);
			}
			ImGui::PopID();
		}
		std::cout << ImGui::GetItemID() << std::endl;
		ImGui::EndListBox();
	}
	
	//ImGui::PopID();
	std::cout << ImGui::GetItemID() << std::endl;
	if (ImGui::BeginPopupContextItem())
	{
		if (ImGui::Button("plot all"))
		{
			SetAllPlotted(true);
		}
		if (ImGui::Button("plot none"))
		{
			SetAllPlotted(false);
		}
		if (ImGui::Button("normalise all"))
		{
			SetAllShowNormalised(true);
		}
		if (ImGui::Button("unnormalise all"))
		{
			SetAllShowNormalised(false);
		}
		ImGui::EndPopup();
	}
	
	ImGui::Text("energy distrubtion set: %s", Label().c_str());

	if (ImGui::Button("save set"))
	{
		SaveHists();
		SaveSamples();
		info.Save(folder / subFolder);
	}
}

void EnergyDistributionSet::SaveSamples() const
{
	// set the output filepath
	std::filesystem::path outfolder = FileUtils::GetEnergyDistSetFolder() / folder / subFolder;                        
	std::cout << outfolder << std::endl;
	// Create the directories if they don't exist
	if (!std::filesystem::exists(outfolder))
	{
		std::filesystem::create_directories(outfolder);
	}

	for (const EnergyDistribution& edist : distributions)
	{
		edist.SaveSamples(outfolder);
	}
}

void EnergyDistributionSet::SaveHists() const
{
	// set the output filepath
	std::filesystem::path outfolder = FileUtils::GetEnergyDistSetFolder() / folder / subFolder;

	// Create the directories if they don't exist
	if (!std::filesystem::exists(outfolder))
	{
		std::filesystem::create_directories(outfolder);
	}

	for (const EnergyDistribution& edist : distributions)
	{
		edist.SaveHist(outfolder);
	}
}

void EnergyDistributionSet::Load(std::filesystem::path& infolder, bool loadSamples)
{
	if (!std::filesystem::exists(infolder) || !std::filesystem::is_directory(infolder))
	{
		std::cerr << "Invalid directory path!" << std::endl;
		return;
	}

	// Iterate through the directory
	for (const auto& entry : std::filesystem::directory_iterator(infolder))
	{
		if (entry.is_regular_file() && entry.path().extension() == ".asc")
		{
			std::filesystem::path file = entry.path();
			EnergyDistribution newDist;
			newDist.Load(file, loadSamples);

			AddDistribution(std::move(newDist));
		}
	}
	folder = infolder.parent_path().parent_path().filename() / infolder.parent_path().filename();
	subFolder = infolder.filename();
}

EnergyDistributionSet::EnergyDistributionSet()
{
	distributions.reserve(100);
}
//
//EnergyDistributionSet::EnergyDistributionSet(EnergyDistributionSet&& other) noexcept
//{
//	distributions = std::move(other.distributions);
//	EdToDistMap = std::move(other.EdToDistMap);
//	info = std::move(other.info);
//	folder = std::move(other.folder);
//	subFolder = std::move(other.subFolder);
//}
//
//EnergyDistributionSet& EnergyDistributionSet::operator=(EnergyDistributionSet&& other) noexcept
//{
//	distributions = std::move(other.distributions);
//	EdToDistMap = std::move(other.EdToDistMap);
//	info = std::move(other.info);
//	folder = std::move(other.folder);
//	subFolder = std::move(other.subFolder);
//
//	return *this;
//}

std::string EnergyDistributionSet::Label()
{
	return (folder / subFolder).string();
}
