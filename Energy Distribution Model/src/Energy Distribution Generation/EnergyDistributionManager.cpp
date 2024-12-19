#include "pch.h"

#include "EnergyDistributionManager.h"
#include "Constants.h"
#include "FileHandler.h"

EnergyDistributionManager::EnergyDistributionManager()
	: EnergyDistributionModule("Energy Distribution Manager")
{
	manager = this;

	m_mainCanvas->cd();
	m_mainCanvas->Clear();
	m_mainCanvas->Divide(3,2);
	m_mainCanvas->SetWindowSize(1500, 800);

	//CreateNewSet();
}

std::vector<EnergyDistributionSet>& EnergyDistributionManager::GetEnergyDistributionSets()
{
	return energyDistributionSets;
}

void EnergyDistributionManager::ShowUI()
{
	ShowSettings();
	ImGui::SameLine();
	ShowTabsWithSets();
	
	//ImGui::Separator();
	ShowEnergyDistributionPlot();
}

void EnergyDistributionManager::ShowSettings()
{
	ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX | ImGuiChildFlags_ResizeY;
	if (ImGui::BeginChild("##Settings", ImVec2(0.0f, 0.0f), flags))
	{
		if (ImGui::Button("select description file"))
		{
			currentDescriptionFile = FileHandler::GetInstance().SelectFile();
			maxIndex = FileHandler::GetInstance().GetMaxIndex(currentDescriptionFile);
		}
		ImGui::SameLine();
		ImGui::Text(currentDescriptionFile.filename().string().c_str());

		ImGui::BeginDisabled(currentDescriptionFile.empty());
		if (ImGui::Button("Generate Distributions from File"))
		{
			GenerateEnergyDistributionsFromFile(currentDescriptionFile);

			PlotEnergyDistributions();
			PLotZweightByEnergy();
			PlotLongkTDistribution();
			PlotLongVelAddition();
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

		ImGui::SetNextItemWidth(200.0f);
		ImGui::BeginDisabled(!activeDist.simplifyParams.cutOutZValues);
		ImGui::InputFloat2("", activeDist.simplifyParams.cutOutRange);
		ImGui::EndDisabled();
		ImGui::SameLine();
		ImGui::Checkbox("cut out z range", activeDist.simplifyParams.cutOutZValues);

		ImGui::SeparatorText("input things");
		ImGui::SeparatorText("output things");
		if(!energyDistributionSets.empty())
		{
			EnergyDistributionSet& currentSet = energyDistributionSets.at(currentSetIndex);
			char buf[64] = "";
			strncpy_s(buf, currentSet.subFolder.string().c_str(), sizeof(buf) - 1);
			ImGui::SetNextItemWidth(150.0f);
			if (ImGui::InputText("set subfolder", buf, IM_ARRAYSIZE(buf)))
			{
				currentSet.subFolder = std::filesystem::path(buf);
			}
		}

		ImGui::SeparatorText("Binning options");
		ImGui::SetNextItemWidth(150.0f);
		ImGui::InputFloat2("energy range", binSettings.energyRange, "%.1e");
		ImGui::SameLine();
		ImGui::Checkbox("more bins at peaks", &binSettings.increasePeakResolution);

		if (ImGui::Checkbox("factor binning", &binSettings.factorBinning))
		{
			binSettings.constantBinSize = !binSettings.factorBinning;
		}
		ImGui::SameLine();
		ImGui::BeginDisabled(!binSettings.factorBinning);
		ImGui::SetNextItemWidth(50.0f);
		ImGui::InputInt("bins per decade", &binSettings.binsPerDecade, 0);
		ImGui::SameLine();
		ImGui::BeginDisabled(!binSettings.increasePeakResolution);
		ImGui::SetNextItemWidth(50.0f);
		ImGui::InputInt("bins at peak", &binSettings.binsAtPeak, 0);
		ImGui::EndDisabled();
		ImGui::EndDisabled();


		if (ImGui::Checkbox("constant bin size", &binSettings.constantBinSize))
		{
			binSettings.factorBinning = !binSettings.constantBinSize;
		}
		ImGui::BeginDisabled(!binSettings.constantBinSize);
		ImGui::SameLine();
		ImGui::SetNextItemWidth(50.0f);
		ImGui::InputDouble("step size", &binSettings.normalStepSize, 0, 0, "%.3f");
		ImGui::SameLine();
		ImGui::BeginDisabled(!binSettings.increasePeakResolution);
		ImGui::SetNextItemWidth(50.0f);
		ImGui::InputDouble("peak step size", &binSettings.peakStepSize, 0, 0, "%.3f");
		ImGui::EndDisabled();
		ImGui::EndDisabled();

		ImGui::SeparatorText("Peak fitting options");
		ImGui::Checkbox("fix kT trans", &fixKT_trans);
		ImGui::Checkbox("fix detuning energy", &fixDetuningEnergy);
		ImGui::Checkbox("show limited fit range", &showLimitedFit);

		ImGui::EndChild();
	}
}

void EnergyDistributionManager::ShowTabsWithSets()
{
	ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeY;
	if (ImGui::BeginChild("listbox", ImVec2(0, 0), flags))
	{
		ImGui::Text("energy distributions");
		ImGui::PushStyleColor(ImGuiCol_TabActive, ImVec4(0.7f, 0.15f, 0.15f, 1.0f));
		ImGui::PushStyleColor(ImGuiCol_TabHovered, ImVec4(0.8f, 0.2f, 0.2f, 1.0f));
		if (ImGui::BeginTabBar("##tab bar", ImGuiTabBarFlags_AutoSelectNewTabs))
		{
			if (ImGui::TabItemButton("+", ImGuiTabItemFlags_Trailing | ImGuiTabItemFlags_NoTooltip))
			{
				CreateNewSet();
			}
			for (int setIndex = 0; setIndex < energyDistributionSets.size(); setIndex++)
			{
				bool open = true;
				std::string label = "set " + std::to_string(setIndex);
				if (ImGui::BeginTabItem(label.c_str(), &open))
				{
					currentSetIndex = setIndex;
					ShowEnergyDistributionSet(setIndex);
					ImGui::EndTabItem();
				}
				if (!open)
				{
					RemoveSet(setIndex);
				}
			}
			ImGui::EndTabBar();
		}
		ImGui::PopStyleColor(2);

		ImGui::SeparatorText("general options");
		if (ImGui::Button("clear plot"))
		{
			for (EnergyDistributionSet& set : energyDistributionSets)
			{
				set.SetAllPlotted(false);
			}
		}
		ImGui::SameLine();
		if (ImGui::Button("load hists"))
		{
			std::vector<std::filesystem::path> filenames = FileHandler::GetInstance().SelectFiles("output\\");
			if (!filenames.empty())
			{
				for (auto& filename : filenames)
				{
					EnergyDistribution energyDist = FileHandler::GetInstance().LoadEnergyDistribution(filename, loadSamples);
					std::filesystem::path subfolder = filename.parent_path().filename();
					std::filesystem::path folder = filename.parent_path().parent_path().filename();
					std::cout << folder << " " << subfolder << std::endl;

					PrepareCurrentSet(folder, subfolder);
					AddDistributionToSet(energyDist, currentSetIndex);
				}
			}
		}
		ImGui::SameLine();
		ImGui::Checkbox("load samples", &loadSamples);

		ImGui::EndChild();
	}
}

void EnergyDistributionManager::ShowEnergyDistributionSet(int setIndex)
{
	EnergyDistributionSet& set = energyDistributionSets.at(setIndex);
	std::vector<EnergyDistribution>& energyDistributionList = set.distributions;

	std::string listboxLabel = std::to_string(setIndex);
	if (ImGui::BeginListBox(listboxLabel.c_str(), ImVec2(-1, 250)))
	{
		for (int i = 0; i < energyDistributionList.size(); i++)
		{
			ImGui::PushID(i);
			EnergyDistribution& eDist = energyDistributionList[i];

			std::string label = eDist.label;
			if (!eDist.tags.empty())
			{
				label += "\n";
				label += eDist.tags;
			}

			// Render each item as selectable
			if (ImGui::Selectable(label.c_str(), eDist.plotted, ImGuiSelectableFlags_AllowItemOverlap))
			{
				eDist.plotted = !eDist.plotted;
			}

			if (ImGui::BeginItemTooltip())
			{
				ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
				ImGui::TextUnformatted(eDist.String().c_str());
				ImGui::PopTextWrapPos();
				ImGui::EndTooltip();
			}

			ImGui::SameLine();
			ImGui::Checkbox("normalised", &eDist.showNormalisedByWidth);
				
			ImGui::SameLine();
			if (ImGui::SmallButton("x"))
			{
				RemoveDistributionFromSet(i, setIndex);
			}

			ImGui::PopID();
		}
		
		ImGui::EndListBox();
	}
	if (ImGui::BeginPopupContextItem())
	{
		//ImGui::SameLine();
		if (ImGui::Button("plot all"))
		{
			set.SetAllPlotted(true);
		}
		//ImGui::SameLine();
		if (ImGui::Button("plot none"))
		{
			set.SetAllPlotted(false);
		}
		//ImGui::SameLine();
		if (ImGui::Button("normalise all"))
		{
			set.SetAllShowNormalised(true);
		}
		//ImGui::SameLine();
		if (ImGui::Button("unnormalise all"))
		{
			set.SetAllShowNormalised(false);
		}
		ImGui::EndPopup();
	}

	ImGui::Text("energy distrubtion set: %s/%s", set.folder.string().c_str(), set.subFolder.string().c_str());	

	if (ImGui::Button("save set"))
	{
		FileHandler::GetInstance().SaveEnergyDistributionSetAsHist(set);
		FileHandler::GetInstance().SaveEnergyDistributionSetAsSamples(set);
	}
}

void EnergyDistributionManager::ShowEnergyDistributionPlot()
{
	ImGui::Checkbox("log X", &logX);
	ImGui::SameLine();
	ImGui::Checkbox("log Y", &logY);
	ImGui::SameLine();
	ImGui::Checkbox("show markers", &showMarkers);
	ImGui::SameLine();
	ImGui::Checkbox("show fit", &showFits);

	if (ImPlot::BeginPlot("collision Energy distribution"))
	{
		ImPlot::SetupAxis(ImAxis_X1, "collision energy [eV]");
		ImPlot::SetupAxis(ImAxis_Y1, "f(E)", ImPlotAxisFlags_AutoFit);
		if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
		if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
		ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);

		int i = 0;
		for (const EnergyDistributionSet& set : energyDistributionSets)
		{
			for (const EnergyDistribution& eDist : set.distributions)
			{
				if (eDist.plotted)
				{
					ImGui::PushID(i);
					if (showMarkers) ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);

					if (eDist.showNormalisedByWidth)
					{
						// Get the automatic color for this pair
						ImVec4 color = ImPlot::GetColormapColor(i % ImPlot::GetColormapSize());

						// Plot the first line with the automatic color
						ImPlot::PushStyleColor(ImPlotCol_Line, color);
						ImPlot::PlotLine(eDist.label.c_str(), eDist.binCenters.data(), eDist.binValuesNormalised.data(), eDist.binValuesNormalised.size());
						ImPlot::PopStyleColor();
						
						if (showFits)
						{
							color.x *= 2;
							color.y *= 2;
							color.z *= 2;

							// Plot the second line with a lighter color and dashed
							ImPlot::PushStyleColor(ImPlotCol_Line, color);
							ImPlot::PlotLine("##", eDist.fitX.data(), eDist.fitY.data(), eDist.fitX.size(), ImPlotLineFlags_Segments);
							ImPlot::PopStyleColor();
						}
					}
					else
					{
						ImPlot::PlotLine(eDist.label.c_str(), eDist.binCenters.data(), eDist.binValues.data(), eDist.binCenters.size());
					}
				}
				i++;
			}
		}
		ImPlot::EndPlot();
	}
}

void EnergyDistributionManager::CreateNewSet()
{
	energyDistributionSets.emplace_back();
	currentSetIndex = energyDistributionSets.size() - 1;
}

void EnergyDistributionManager::RemoveSet(int setIndex)
{
	ClearDistributionsInSet(setIndex);
	energyDistributionSets.erase(energyDistributionSets.begin() + setIndex);
	currentSetIndex = std::min(currentSetIndex, (int)energyDistributionSets.size() - 1);
}

void EnergyDistributionManager::AddDistributionToSet(EnergyDistribution& distribution, int setIndex)
{
	if (energyDistributionSets.empty())
	{
		CreateNewSet();
		setIndex = 0;
	}
	EnergyDistributionSet& set = energyDistributionSets.at(setIndex);
	std::vector<EnergyDistribution>& list = set.distributions;

	// will call move Constructor
	list.emplace_back(std::move(distribution));
	EnergyDistribution& justMoved = list.back();
	//std::cout << "E_d: " << justMoved.eBeamParameter.detuningEnergy << std::endl;
	set.EdToDistMap[justMoved.eBeamParameter.detuningEnergy] = &justMoved;
}

void EnergyDistributionManager::RemoveDistributionFromSet(int index, int setIndex)
{
	EnergyDistributionSet& set = energyDistributionSets.at(setIndex);
	std::vector<EnergyDistribution>& list = set.distributions;

	// edist needs to be removed from map
	for (auto it = set.EdToDistMap.begin(); it != set.EdToDistMap.end(); it++)
	{
		if (it->second == &list.at(index))
		{
			set.EdToDistMap.erase(it);
			break;
		}
	}
	list.erase(list.begin() + index);
}

void EnergyDistributionManager::PrepareCurrentSet(std::filesystem::path folder, std::filesystem::path subfolder)
{
	if (energyDistributionSets.empty())
	{
		CreateNewSet();
	}
	EnergyDistributionSet& currentSet = energyDistributionSets.at(currentSetIndex);
	if (currentSet.distributions.empty())
	{
		currentSet.folder = folder;
		if(!subfolder.empty()) currentSet.subFolder = subfolder;
	}
	else if (currentSet.folder != folder || (currentSet.subFolder != subfolder && !subfolder.empty()))
	{
		CreateNewSet();
		energyDistributionSets.at(currentSetIndex).folder = folder;
		if (!subfolder.empty()) energyDistributionSets.at(currentSetIndex).subFolder = subfolder;
	}
}

void EnergyDistributionManager::GenerateEnergyDistribution()
{
	// sample positions from electron density multiplied with ion density given from outside
	std::vector<Point3D> positionSamples = mcmc->GetSamples();
	activeDist.collisionEnergies.reserve(positionSamples.size());

	if (positionSamples.empty())
	{
		std::cout << "no sampled positions were given\n";
		return; 
	}
	
	if (!labEnergies->GetDistribution())
	{
		std::cout << "no lab energies were given\n";
		return; 
	}

	SetupSecondaryPlots();

	for (const Point3D& point : positionSamples)
	{
		double x = point.x;
		double y = point.y;
		double z = point.z;
		if (activeDist.simplifyParams.cutOutZValues)
		{
			if (z < activeDist.simplifyParams.cutOutRange.get().x ||
				z > activeDist.simplifyParams.cutOutRange.get().y)
				continue;
		}
		
		// calculate velocity magnitude from lab energy given as matrix (TH3D) from outside
		double labEnergy = labEnergies->Get(x, y, z);
		double electronVelocityMagnitude = TMath::Sqrt(2 * labEnergy * TMath::Qe() / PhysicalConstants::electronMass);

		zPositions->Fill(z);
		zWeightByEnergy->Fill(z, labEnergy);

		// determine direction of velocity based on beam trajectory function
		TVector3 longitudinalDirection = eBeam->GetDirection(point.z);
		TVector3 transverseDirection = longitudinalDirection.Orthogonal();

		// add random values to velocity in transverse and longitudinal directions:
		// - calculate longitudinal kT, transverse kT is fixed
		double long_kT = eBeam->GetLongitudinal_kT(labEnergy);
		double trans_kT = eBeam->GetTransverse_kT();
		long_ktDistribution->Fill(long_kT);

		// - use kT to calculate sigmas of gaussians
		double longSigma = TMath::Sqrt(long_kT * TMath::Qe() / PhysicalConstants::electronMass);
		double transSigma = TMath::Sqrt(trans_kT * TMath::Qe() / PhysicalConstants::electronMass);

		// - sample from gaussians with these sigmas and add that to the electron velocity
		longitudinalNormalDistribution = std::normal_distribution<double>(0, longSigma);
		transverseNormalDistribution = std::normal_distribution<double>(0, transSigma);
		double longitudinalAddition = longitudinalNormalDistribution(generator);
		double transverseAddition = transverseNormalDistribution(generator);
		double transverseAdditionAngle = angleDistribution(generator);
		long_VelAddition->Fill(longitudinalAddition);

		transverseDirection.Rotate(transverseAdditionAngle, longitudinalDirection);

		TVector3 finalElectronVelocity = electronVelocityMagnitude * longitudinalDirection
										 + longitudinalAddition * longitudinalDirection
										 + sqrt(2) * transverseAddition * transverseDirection;

		// calculate collision velocity vector and magnitude using a fixed ion beam velocity
		double ionVelocityMagnitude = TMath::Sqrt(2 * activeDist.eBeamParameter.coolingEnergy * TMath::Qe() / PhysicalConstants::electronMass); // calc from cooling energy;
		TVector3 ionVelocity(0, 0, ionVelocityMagnitude);

		double collisionVelocity = (finalElectronVelocity - ionVelocity).Mag();

		// calculate collision energy [eV] and put it in a histogram
		double collisionEnergy = 0.5 * PhysicalConstants::electronMass * collisionVelocity * collisionVelocity / TMath::Qe();
		activeDist.Fill(collisionEnergy);
		activeDist.collisionEnergies.push_back(collisionEnergy);
	}

	activeDist.FillVectorsFromHist();
	activeDist.RemoveEdgeZeros();
	activeDist.CalculateFWHM();
	activeDist.FitAnalyticalToPeak(fixKT_trans, fixDetuningEnergy, showLimitedFit);

	//std::cout << "Ed1: " << activeDist.eBeamParameter.detuningEnergy << "\n";
}

void EnergyDistributionManager::GenerateEnergyDistributionsFromFile(std::filesystem::path file)
{
	// get all necessary modules
	FileHandler fileHandler = FileHandler::GetInstance();
	std::filesystem::path folder = file.parent_path().filename();
	std::cout << folder << std::endl;

	int end = endIndex;
	int start = startIndex;
	if (doAll)
	{
		start = 1;
		end = maxIndex;
	}

	for (int index = start; index <= end; index++)
	{
		// get 3 parameters: U drift tube, electron current, center E lab if index is in file
		std::array<float, 3> additionalParameter = fileHandler.GetParamtersFromDescriptionFileAtIndex(file, index);

		// if they are not found the index is not in the file
		if (!additionalParameter[0]) continue;

		// set read electron current and center lab energy
		activeDist.labEnergiesParameter.driftTubeVoltage = additionalParameter[0];
		activeDist.eBeamParameter.electronCurrent = additionalParameter[1];
		activeDist.labEnergiesParameter.centerLabEnergy = additionalParameter[2];
		activeDist.eBeamParameter.longitudinal_kT_estimate = eBeam->GetLongitudinal_kT(additionalParameter[2]);
		eBeam->CalculateDetuningEnergy();
		
		// full procedure to generate one energy distribution 
		// 1. setup necessary distributions
		std::filesystem::path densityfile = fileHandler.FindFileWithIndex(file.parent_path() / "e-densities", index);
		if (densityfile.empty()) continue;
		eBeam->SetupDistribution(densityfile);

		std::filesystem::path energyfile = FileHandler::GetInstance().FindFileWithIndex(file.parent_path() / "lab-energies", index);
		if (energyfile.empty()) continue;
		labEnergies->SetupDistribution(energyfile);

		ionBeam->SetupDistribution();
		mcmc->SetupDistribution();

		// 2. sample from this distribution
		mcmc->GenerateSamples();
		
		// final setup of current distribution
		activeDist.SetupLabellingThings();
		activeDist.SetupBinning(binSettings);

		// 3. generate energy distribution
		GenerateEnergyDistribution();

		// store/move and save current distribution that has been worked on
		PrepareCurrentSet(folder);
		AddDistributionToSet(activeDist, currentSetIndex);
	}
}

void EnergyDistributionManager::ShowSetListWindow()
{
	if (ImGui::Begin("energy distribution sets"))
	{
		if (ImGui::BeginListBox("##setlist", ImVec2(-1, 150)))
		{
			for (int i = 0; i < energyDistributionSets.size(); i++)
			{
				EnergyDistributionSet& set = energyDistributionSets.at(i);
				std::string label = (set.folder / set.subFolder).string();

				ImGui::PushID(i);
				bool selected = i == currentSetIndex;
				if (ImGui::Selectable(label.c_str(), selected, ImGuiSelectableFlags_AllowItemOverlap))
				{
					currentSetIndex = i;
				}
				ImGui::SameLine();
				if (ImGui::SmallButton("x"))
				{
					RemoveSet(i);
				}
				ImGui::PopID();

			}
			ImGui::EndListBox();
		}
		if (ImGui::Button("load energy distribution set"))
		{
			std::filesystem::path folder = FileHandler::GetInstance().SelectFolder("output\\Energy Distribution Sets\\");
			if (!folder.empty())
			{
				EnergyDistributionSet set = FileHandler::GetInstance().LoadEnergyDistributionSet(folder);
				if (set.distributions.empty())
				{
					std::cout << "there was no energy distribution in that folder" << std::endl;
				}
				else
				{
					energyDistributionSets.emplace_back(std::move(set));
					currentSetIndex = energyDistributionSets.size() - 1;
				}
			}
		}

		ImGui::End();
	}
}

void EnergyDistributionManager::SetupSecondaryPlots()
{
	double kTLongGuess = eBeam->GetLongitudinal_kT(activeDist.labEnergiesParameter.centerLabEnergy);
	double sigmaGuess = TMath::Sqrt(kTLongGuess * TMath::Qe() / PhysicalConstants::electronMass);
	std::cout << "long kT guess: " << kTLongGuess << "\n";

	delete zPositions;
	delete zWeightByEnergy;
	delete long_ktDistribution;
	delete long_VelAddition;
	zPositions = new TH1D("z-positions", "z-positions", 100, 0, 0.65);
	zWeightByEnergy = new TH1D("z weight by energy", "z weight by energy", 100, 0, 0.65);
	long_ktDistribution = new TH1D("long kT", "long kT", 1000, 0, kTLongGuess * 3);
	long_VelAddition = new TH1D("long vel add", "long vel add", 500, -4 * sigmaGuess, 4 * sigmaGuess);

}

void EnergyDistributionManager::PlotEnergyDistributions()
{
	m_mainCanvas->cd(3);
	int colors[5] = { kRed, kBlue, kGreen, kOrange, kMagenta };

	gPad->SetLogy();
	gPad->SetLogx();

	// Create a legend
	TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
	int i = 0;
	for (EnergyDistributionSet& set : energyDistributionSets)
	{
		for (EnergyDistribution& eDist : set.distributions)
		{
			//if (!energyDistributionList[i]) return;

			eDist.SetLineColor(colors[i % 5]);
			legend->AddEntry(&eDist, eDist.label.c_str(), "l");

			if (i == 0)
			{
				eDist.Draw("HIST");
			}
			else
			{
				eDist.Draw("HIST SAME");
			}
			i++;
		}
	}
	
	legend->Draw();
}

void EnergyDistributionManager::PLotZweightByEnergy()
{
	if (!zWeightByEnergy) return;

	m_secondCanvas->cd(4);
	zPositions->Draw();

	m_secondCanvas->cd(5);
	zWeightByEnergy->Divide(zPositions);
	zWeightByEnergy->Draw("Hist");
}

void EnergyDistributionManager::ClearDistributionsInSet(int setIndex)
{
	std::vector<EnergyDistribution>& list = energyDistributionSets.at(setIndex).distributions;

	for (int i = list.size() - 1; i >= 0; i--)
	{
		RemoveDistributionFromSet(i, setIndex);
	}
}

void EnergyDistributionManager::PlotLongkTDistribution()
{
	if (!long_ktDistribution) return;

	m_secondCanvas->cd(1);

	gPad->SetLogy();
	gPad->SetLogx();

	long_ktDistribution->Draw();
}

void EnergyDistributionManager::PlotLongVelAddition()
{
	if (!long_VelAddition) return;

	m_secondCanvas->cd(2);

	long_VelAddition->Draw();
}

