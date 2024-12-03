#include "pch.h"

#include "EnergyDistributionManager.h"
#include "Constants.h"
#include "FileHandler.h"

EnergyDistributionManager::EnergyDistributionManager()
	: EnergyDistributionModule("Energy Distribution Manager")
{
	m_mainCanvas->cd();
	m_mainCanvas->Clear();
	m_mainCanvas->Divide(3,2);
	m_mainCanvas->SetWindowSize(1500, 800);

	energyDistributions.reserve(5);
}

std::vector<EnergyDistribution>& EnergyDistributionManager::GetEnergyDistributions()
{
	return energyDistributions;
}

void EnergyDistributionManager::ShowUI()
{
	ShowSettings();
	ImGui::SameLine();
	ShowEnergyDistributionList();
	
	//ImGui::Separator();
	ShowEnergyDistributionPlot();
}

void EnergyDistributionManager::ShowSettings()
{
	ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX | ImGuiChildFlags_ResizeY;
	if (ImGui::BeginChild("##Settings", ImVec2(0.0f, 0.0f), flags))
	{
		if (ImGui::Button("Generate Single Energy Distribution"))
		{
			GenerateEnergyDistribution();

			PlotEnergyDistributions();
			PLotZweightByEnergy();
		}

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
		ImGui::SameLine();
		if (ImGui::Checkbox("save as hist", &saveAsHist))
		{
			if (!saveAsHist) saveSamplesToFile = false;
		}
		ImGui::SameLine();
		if (ImGui::Checkbox("save samples", &saveSamplesToFile))
		{
			if (saveSamplesToFile) saveAsHist = true;
		}

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

		//ImGui::SeparatorText("analytical distribution");
		//ImGui::PushID("analytical energy distribution");
		//if (ImGui::Button("generate analytical distribution"))
		//{
		//	GenerateAnalyticalDistribution();
		//}
		//ImGui::SetNextItemWidth(200.0f);
		//ImGui::InputFloat2("energy range", analyticalEnergyRange, "%.1e");
		//ImGui::SameLine();
		//ImGui::SetNextItemWidth(100.0f);
		//ImGui::InputInt("number bins", &analyticalNumberBins);
		//ImGui::SetNextItemWidth(100.0f);
		//ImGui::InputDouble("detuning energy [eV]", analyticalParameter.detuningEnergy);
		//ImGui::SetNextItemWidth(100.0f);
		//ImGui::InputDouble("transverse kT [eV]", analyticalParameter.transverseTemperature);
		//ImGui::SetNextItemWidth(100.0f);
		//ImGui::InputDouble("longitudinal kT [eV]", analyticalParameter.longitudinalTemperature);
		//ImGui::PopID();

		ImGui::EndChild();
	}
}

void EnergyDistributionManager::ShowEnergyDistributionList()
{
	ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeY;
	if (ImGui::BeginChild("##listbox", ImVec2(0, 0), flags))
	{
		ImGui::Text("loaded distributions");
		if (ImGui::BeginListBox("##", ImVec2(-1, 250)))
		{
			for (int i = 0; i < energyDistributions.size(); i++)
			{
				ImGui::PushID(i);
				EnergyDistribution& eDist = energyDistributions[i];

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
					RemoveDistributionFromList(i);
				}

				ImGui::PopID();
			}
			ImGui::EndListBox();
		}
		
		if (ImGui::Button("load hists"))
		{
			std::vector<std::filesystem::path> filenames = FileHandler::GetInstance().SelectFiles("output\\");
			if (!filenames.empty())
			{
				for (auto& filename : filenames)
				{
					EnergyDistribution energyDist = FileHandler::GetInstance().LoadEnergyDistribution(filename, loadSamples);
					AddDistributionToList(std::move(energyDist));
				}
			}
		}
		ImGui::SameLine();
		ImGui::Checkbox("load samples", &loadSamples);

		if (ImGui::Button("clear list"))
		{
			ClearDistributionList();
		}
		ImGui::SameLine();
		if (ImGui::Button("clear plot"))
		{
			for (EnergyDistribution& eDist : energyDistributions)
			{
				eDist.plotted = false;
			}
		}
		ImGui::SameLine();
		if (ImGui::Button("plot all"))
		{
			for (EnergyDistribution& eDist : energyDistributions)
			{
				eDist.plotted = true;
			}
		}
		ImGui::SameLine();
		if (ImGui::Button("normalise all"))
		{
			for (EnergyDistribution& eDist : energyDistributions)
			{
				eDist.showNormalisedByWidth = true;
			}
		}
		ImGui::SameLine();
		if (ImGui::Button("unnormalise all"))
		{
			for (EnergyDistribution& eDist : energyDistributions)
			{
				eDist.showNormalisedByWidth = false;
			}
		}
		ImGui::EndChild();
	}
}

void EnergyDistributionManager::ShowEnergyDistributionPlot()
{
	ImGui::Checkbox("log X", &logX);
	ImGui::SameLine();
	ImGui::Checkbox("log Y", &logY);
	ImGui::SameLine();
	ImGui::Checkbox("show markers", &showMarkers);

	if (ImPlot::BeginPlot("collision Energy distribution"))
	{
		ImPlot::SetupAxis(ImAxis_X1, "collision energy [eV]");
		ImPlot::SetupAxis(ImAxis_Y1, "f(E)", ImPlotAxisFlags_AutoFit);
		if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
		if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
		ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);
		int i = 0;
		for (EnergyDistribution& eDist : energyDistributions)
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
					color.x *= 2;
					color.y *= 2;
					color.z *= 2;

					// Plot the second line with a lighter color and dashed
					ImPlot::PushStyleColor(ImPlotCol_Line, color);
					ImPlot::PlotLine("##", eDist.fitX.data(), eDist.fitY.data(), eDist.fitX.size(), ImPlotLineFlags_Segments);
					ImPlot::PopStyleColor();
				}
				else
				{
					ImPlot::PlotLine(eDist.label.c_str(), eDist.binCenters.data(), eDist.binValues.data(), eDist.binCenters.size());
				}
			}
			i++;
		}
		ImPlot::EndPlot();
	}
}


void EnergyDistributionManager::AddDistributionToList(EnergyDistribution&& distribution)
{
	// will call move Constructor
	energyDistributions.emplace_back(std::move(distribution));
	EnergyDistribution& justMoved = energyDistributions.back();
	EnergyDistribution::s_allDistributions[justMoved.eBeamParameter.detuningEnergy] = &justMoved;
}

void EnergyDistributionManager::RemoveDistributionFromList(int index)
{
	//delete energyDistributions[index];
	energyDistributions.erase(energyDistributions.begin() + index);
}

void EnergyDistributionManager::GenerateEnergyDistribution()
{
	// final setup of current distribution
	activeDist.SetupLabellingThings();
	activeDist.SetupBinning(binSettings);

	MCMC* mcmc = (MCMC*)EnergyDistributionModule::Get("MCMC");
	ElectronBeam* eBeam = (ElectronBeam*)EnergyDistributionModule::Get("Electron Beam");
	LabEnergies* labEnergies = (LabEnergies*)EnergyDistributionModule::Get("Lab Energies");

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
	activeDist.FitAnalyticalToPeak();

	std::cout << "Ed1: " << activeDist.eBeamParameter.detuningEnergy << "\n";

	if (saveAsHist)
	{
		FileHandler::GetInstance().SaveEnergyDistributionHistToFile(activeDist);
	}
	
	if (saveSamplesToFile)
	{
		FileHandler::GetInstance().SaveEnergyDistributionSamplesToFile(activeDist);
	}

	// store/move and save current distribution that has been worked on
	AddDistributionToList(std::move(activeDist));

	// create new distribution object that will be worked on now
	//activeDist = new EnergyDistribution();
}

void EnergyDistributionManager::GenerateEnergyDistributionsFromFile(std::filesystem::path file)
{
	// get all necessary modules
	FileHandler fileHandler = FileHandler::GetInstance();
	ElectronBeam* eBeam = (ElectronBeam*)EnergyDistributionModule::Get("Electron Beam");
	IonBeam* ionBeam = (IonBeam*)EnergyDistributionModule::Get("Ion Beam");
	MCMC* mcmc = (MCMC*)EnergyDistributionModule::Get("MCMC");
	LabEnergies* labEnergies = (LabEnergies*)EnergyDistributionModule::Get("Lab Energies");

	int end = endIndex;
	int start = startIndex;
	if (doAll)
	{
		start = 1;
		end = maxIndex;
	}

	for (int index = start; index <= end; index++)
	{
		// get 3 parameters: U drift tube, electron current, center E lab
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
		
		// 3. generate energy distribution
		GenerateEnergyDistribution();		
	}
}

void EnergyDistributionManager::SetupSecondaryPlots()
{
	ElectronBeam* eBeam = (ElectronBeam*)EnergyDistributionModule::Get("Electron Beam");

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

	for (int i = 0; i < energyDistributions.size(); i++)
	{
		//if (!energyDistributions[i]) return;
		
		energyDistributions[i].SetLineColor(colors[i % 5]);
		legend->AddEntry(&energyDistributions[i], energyDistributions[i].label.c_str(), "l");

		if (i == 0)
		{
			energyDistributions[i].Draw("HIST");
		}
		else
		{
			energyDistributions[i].Draw("HIST SAME");
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

void EnergyDistributionManager::ClearDistributionList()
{
	for (int i = energyDistributions.size() - 1; i >= 0; i--)
	{
		RemoveDistributionFromList(i);
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

