#include "pch.h"

#include "EnergyDistributionWindow.h"

#include "MCMC.h"
#include "ElectronBeam.h"
#include "IonBeam.h"
#include "LabEnergies.h"
#include "EnergyDistribution.h"
#include "AnalyticalDistribution.h"

#include "Constants.h"
#include "FileUtils.h"

#include "CoolingForceModel.h"

namespace EnergyDistributionWindow
{
	// main data storage
	static std::vector<EnergyDistributionSet> setList;
	static int currentSetIndex;

	// currently loaded description file
	static std::filesystem::path currentDescriptionFile = std::filesystem::path("data\\dataset1\\100x100x100_Ie0.95_Ucath44.2_RelTol0_mbrc1_energies.asc");
	static int maxIndex = 0;

	// start/end index in description file to generate distribution for
	static int startIndex = 1;
	static int endIndex = 1;
	static bool doAll = false;

	// save all sampled values in a file to load it again
	static bool loadSamples = true;

	// set information window things
	static bool showSetInformation = false;
	static bool infoPlotsLogX = false;
	static bool infoPlotsLogY = false;

	// binning related parameters
	static bool showBinningSettings = false;
	static BinningSettings binSettings;

	// fit options
	static bool showPeakFitSettings = false;
	static PeakFitSettings peakFitSettings;

	// plot parameters
	static bool logX = true;
	static bool logY = true;
	static bool showMarkers = false;
	static bool showFits = true;
	static bool showAnalytical = false;

	static bool showAllParamsWindow = false;

	void Init()
	{
		AnalyticalDistribution::Update();
	}

	EnergyDistributionSet& GetCurrentSet()
	{
		return setList.at(currentSetIndex);
	}

	std::vector<EnergyDistributionSet>& GetSetList()
	{
		return setList;
	}

	int GetCurrentSetIndex()
	{
		return currentSetIndex;
	}

	void CreateNewSet()
	{
		setList.emplace_back();
		currentSetIndex = setList.size() - 1;
	}

	void SetupSet(std::filesystem::path folder, std::filesystem::path subfolder)
	{
		if (setList.empty())
		{
			CreateNewSet();
		}
		EnergyDistributionSet& currentSet = setList.at(currentSetIndex);
		if (currentSet.GetDistributions().empty())
		{
			currentSet.SetFolder(folder);
			if (!subfolder.empty()) currentSet.SetSubfolder(subfolder);
		}
		else if (currentSet.GetFolder() != folder || (currentSet.GetSubfolder() != subfolder && !subfolder.empty()))
		{
			CreateNewSet();
			setList.at(currentSetIndex).SetFolder(folder);
			if (!subfolder.empty()) setList.at(currentSetIndex).SetSubfolder(subfolder);
		}
	}

	void ShowWindow()
	{
		if (ImGui::Begin("Energy Distribution Window"))
		{
			ShowTabs();

			ImGui::SameLine();

			ImGui::BeginGroup();
			ShowSettings();
			ShowPlot();
			ImGui::EndGroup();

			//ShowSetInformationWindow();
			AnalyticalDistribution::ShowWindow(showAnalytical);
			peakFitSettings.ShowWindow(showPeakFitSettings);
			binSettings.ShowWindow(showBinningSettings);
			ShowAllParamtersWindow();
			ShowSetInformationWindow();

			ImGui::End();
		}
	}

	void ShowSettings()
	{
		ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX | ImGuiChildFlags_ResizeY;
		if (ImGui::BeginChild("##Settings", ImVec2(0.0f, 0.0f), flags))
		{
			ImGui::SeparatorText("input things");
			if (ImGui::Button("select description file"))
			{
				currentDescriptionFile = FileUtils::SelectFile();
				maxIndex = FileUtils::GetMaxIndex(currentDescriptionFile);
			}
			ImGui::Text("file: %s", (currentDescriptionFile.parent_path().filename() / currentDescriptionFile.filename()).string().c_str());

			ImGui::Checkbox("Binning settings", &showBinningSettings);
			ImGui::SameLine();
			ImGui::Checkbox("All Parameters", &showAllParamsWindow);


			ImGui::SeparatorText("output things");
			ImGui::Checkbox("Peak fit settings", &showPeakFitSettings);

			if (!setList.empty())
			{
				EnergyDistributionSet& currentSet = setList.at(currentSetIndex);
				char buf[64] = "";
				strncpy_s(buf, currentSet.GetSubfolder().string().c_str(), sizeof(buf) - 1);
				ImGui::SetNextItemWidth(150.0f);
				if (ImGui::InputText("set subfolder", buf, IM_ARRAYSIZE(buf)))
				{
					currentSet.SetSubfolder(std::filesystem::path(buf));
				}
			}

			ImGui::Separator();

			ImGui::BeginDisabled(currentDescriptionFile.empty());
			if (ImGui::Button("Generate Distributions from File"))
			{
				std::filesystem::path folder = currentDescriptionFile.parent_path().filename();
				SetupSet(folder);

				for (int i = startIndex; i <= endIndex; i++)
				{
					EnergyDistribution newDist;
					newDist.Generate(currentDescriptionFile, i, binSettings, peakFitSettings);
					EnergyDistributionSet& set = GetCurrentSet();
					set.AddDistribution(std::move(newDist));
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

			ImGui::EndChild();
		}
	}

	void ShowTabs()
	{
		ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX | ImGuiChildFlags_ResizeY;
		if (ImGui::BeginChild("listbox", ImVec2(400.0f, -1), flags))
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
				for (int setIndex = 0; setIndex < setList.size(); setIndex++)
				{
					bool open = true;
					std::string label = "set " + std::to_string(setIndex);
					if (ImGui::BeginTabItem(label.c_str(), &open))
					{
						currentSetIndex = setIndex;
						setList.at(setIndex).ShowList();
						ImGui::EndTabItem();
					}
					if (!open)
					{
						setList.erase(setList.begin() + setIndex);
						currentSetIndex = std::min(currentSetIndex, (int)setList.size() - 1);
					}
				}
				ImGui::EndTabBar();
			}
			ImGui::PopStyleColor(2);

			ImGui::SeparatorText("general options");
			if (ImGui::Button("clear plot"))
			{
				for (EnergyDistributionSet& set : setList)
				{
					set.SetAllPlotted(false);
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("load hists"))
			{
				std::vector<std::filesystem::path> filenames = FileUtils::SelectFiles(FileUtils::GetEnergyDistSetFolder());
				if (!filenames.empty())
				{
					for (auto& filename : filenames)
					{
						EnergyDistribution energyDist;
						energyDist.Load(filename, loadSamples);
						std::filesystem::path subfolder = filename.parent_path().filename();
						std::filesystem::path folder = filename.parent_path().parent_path().filename();

						SetupSet(folder, subfolder);
						GetCurrentSet().AddDistribution(std::move(energyDist));
					}
				}
			}
			ImGui::SameLine();
			ImGui::Checkbox("load samples", &loadSamples);
			ImGui::SameLine();
			ImGui::Checkbox("show set information", &showSetInformation);
			ImGui::EndChild();
		}
	}

	void ShowPlot()
	{
		ImGui::Checkbox("log X", &logX);
		ImGui::SameLine();
		ImGui::Checkbox("log Y", &logY);
		ImGui::SameLine();
		ImGui::Checkbox("show markers", &showMarkers);
		ImGui::SameLine();
		ImGui::Checkbox("show fit", &showFits);
		ImGui::SameLine();
		ImGui::Checkbox("show analytical", &showAnalytical);

		if (ImPlot::BeginPlot("collision Energy distribution", ImVec2(-1, -1)))
		{
			ImPlot::SetupAxis(ImAxis_X1, "collision energy [eV]");
			ImPlot::SetupAxis(ImAxis_Y1, "f(E)", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit);
			if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
			if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
			ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);

			int i = 0;
			for (const EnergyDistributionSet& set : setList)
			{
				for (const EnergyDistribution& eDist : set.GetDistributions())
				{
					ImGui::PushID(i++);
					eDist.Plot(showMarkers, showFits);
					ImGui::PopID();
				}
			}

			if (showAnalytical)
			{
				AnalyticalDistribution::Plot();
			}
			ImPlot::EndPlot();
		}
	}

	void ShowAllParamtersWindow()
	{
		if (showAllParamsWindow && ImGui::Begin("all parameters", &showAllParamsWindow, ImGuiWindowFlags_NoDocking))
		{
			ImGuiChildFlags flags = ImGuiChildFlags_ResizeX | ImGuiChildFlags_Border;
			if (ImGui::BeginChild("mcmc", ImVec2(100, -1), flags))
			{
				ImGui::SeparatorText("MCMC parameter");
				MCMC::ShowParameterControls();
				ImGui::EndChild();
			}
			ImGui::SameLine();
			if (ImGui::BeginChild("labe", ImVec2(100, -1), flags))
			{
				ImGui::SeparatorText("Lab Energy parameter");
				LabEnergy::ShowParameterControls();
				ImGui::EndChild();
			}
			ImGui::SameLine();
			if (ImGui::BeginChild("ebeam", ImVec2(100, -1), flags))
			{
				ImGui::SeparatorText("Electron Beam parameter");
				ElectronBeam::ShowParameterControls();
				ImGui::EndChild();
			}
			ImGui::SameLine();
			if (ImGui::BeginChild("ibeam", ImVec2(100, -1), flags))
			{
				ImGui::SeparatorText("Ion Beam parameter");
				IonBeam::ShowParameterControls();
				ImGui::EndChild();
			}
			
			
			ImGui::End();
		}
	}

	void ShowSetInformationWindow()
	{
		if (showSetInformation && ImGui::Begin("Set Information", &showSetInformation, ImGuiWindowFlags_NoDocking))
		{
			ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX;
			if (ImGui::BeginChild("Selection", ImVec2(100, -1), flags))
			{
				int i = 0;
				for (EnergyDistributionSet& set : setList)
				{
					ImGui::PushID(i++);
					ImGui::Checkbox(set.Label().c_str(), &set.GetInfo().plot);
					ImGui::PopID();
				}

				ImGui::Checkbox("log X", &infoPlotsLogX);
				ImGui::Checkbox("log Y", &infoPlotsLogY);
				ImGui::EndChild();
			}
			ImGui::SameLine();

			int i = 0;
			ImPlotSubplotFlags subplotflags = ImPlotSubplotFlags_NoTitle | ImPlotSubplotFlags_ShareItems;
			if (ImPlot::BeginSubplots("set info", 2, 3, ImVec2(-1, -1), subplotflags))
			{
				if (ImPlot::BeginPlot("fit E_d"))
				{
					ImPlot::SetupAxes("detuning energy [eV]", "fit E_d");
					if (infoPlotsLogX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
					if (infoPlotsLogY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

					int i = 0;
					for (EnergyDistributionSet& set : setList)
					{
						set.GetInfo().PlotFitEd(set.Label() + "##" + std::to_string(i++));
					}
					
					ImPlot::EndPlot();
				}
				if (ImPlot::BeginPlot("fit kT_long"))
				{
					ImPlot::SetupAxes("detuning energy [eV]", "fit kT_long");
					if (infoPlotsLogX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
					if (infoPlotsLogY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

					int i = 0;
					for (EnergyDistributionSet& set : setList)
					{
						set.GetInfo().PlotFitlongkT(set.Label() + "##" + std::to_string(i++));
					}

					ImPlot::EndPlot();
				}
				if (ImPlot::BeginPlot("fit kT_trans"))
				{
					ImPlot::SetupAxes("detuning energy [eV]", "fit kT_trans");
					if (infoPlotsLogX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
					if (infoPlotsLogY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

					int i = 0;
					for (EnergyDistributionSet& set : setList)
					{
						set.GetInfo().PlotFitTranskT(set.Label() + "##" + std::to_string(i++));
					}

					ImPlot::EndPlot();
				}
				if (ImPlot::BeginPlot("fit scaling factor"))
				{
					ImPlot::SetupAxes("detuning energy [eV]", "fit scaling factor");
					if (infoPlotsLogX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
					if (infoPlotsLogY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

					int i = 0;
					for (EnergyDistributionSet& set : setList)
					{
						set.GetInfo().PlotFitScalingFactor(set.Label() + "##" + std::to_string(i++));
					}

					ImPlot::EndPlot();
				}
				if (ImPlot::BeginPlot("fit FWHM"))
				{
					ImPlot::SetupAxes("detuning energy [eV]", "fit FWHM");
					if (infoPlotsLogX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
					if (infoPlotsLogY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

					int i = 0;
					for (EnergyDistributionSet& set : setList)
					{
						set.GetInfo().PlotFitFWHM(set.Label() + "##" + std::to_string(i++));
					}

					ImPlot::EndPlot();
				}
				if (ImPlot::BeginPlot("FWHM"))
				{
					ImPlot::SetupAxes("detuning energy [eV]", "FWHM");
					if (infoPlotsLogX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
					if (infoPlotsLogY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

					int i = 0;
					for (EnergyDistributionSet& set : setList)
					{
						set.GetInfo().PlotFWHM(set.Label() + "##" + std::to_string(i++));
					}

					ImPlot::EndPlot();
				}
			
				ImPlot::EndSubplots();
			}
			ImGui::End();
		}
	}

	void ShowSetList()
	{
		ImGuiChildFlags flags = ImGuiChildFlags_Border | ImGuiChildFlags_ResizeX | ImGuiChildFlags_ResizeY;
		if (ImGui::BeginChild("energy distribution sets", ImVec2(100, 100), flags))
		{
			if (ImGui::BeginListBox("##setlist", ImVec2(-1, 150)))
			{
				for (int i = 0; i < setList.size(); i++)
				{
					EnergyDistributionSet& set = setList.at(i);
					std::string label = set.Label();

					ImGui::PushID(i);
					bool selected = i == currentSetIndex;
					if (ImGui::Selectable(label.c_str(), selected, ImGuiSelectableFlags_AllowItemOverlap))
					{
						currentSetIndex = i;
					}
					ImGui::SameLine();
					if (ImGui::SmallButton("x"))
					{
						setList.erase(setList.begin() + i);
						currentSetIndex = std::min(currentSetIndex, (int)setList.size() - 1);
					}
					ImGui::PopID();

				}
				ImGui::EndListBox();
			}
			if (ImGui::Button("load energy distribution set"))
			{
				std::filesystem::path folder = FileUtils::SelectFolder(FileUtils::GetEnergyDistSetFolder());
				if (!folder.empty())
				{
					EnergyDistributionSet set;
					set.Load(folder, loadSamples);
					if (set.GetDistributions().empty())
					{
						std::cout << "there was no energy distribution in that folder" << std::endl;
					}
					else
					{
						setList.emplace_back(std::move(set));
						currentSetIndex = setList.size() - 1;
					}
				}
			}

			ImGui::EndChild();
		}
	}
}

//
//void EnergyDistributionManager::ShowCoolingForceWindow()
//{
//	if (showCoolingForce && ImGui::Begin("Cooling force", &showCoolingForce, ImGuiWindowFlags_NoDocking))
//	{
//		if (currentSetIndex >= 0 && energyDistributionSets.at(currentSetIndex).distributions.size() > 0)
//		{
//			CoolingForceData& cfData = energyDistributionSets.at(currentSetIndex).distributions.at(0).cfData;
//
//			if (ImGui::SliderFloat("z value", &zValue, -0.7, 0.7))
//			{
//				cfData.forceXSlice.FromTH3D(cfData.forceX, zValue);
//				cfData.forceYSlice.FromTH3D(cfData.forceY, zValue);
//				cfData.forceZSlice.FromTH3D(cfData.forceZ, zValue);
//			}
//
//			if (ImPlot::BeginPlot("cool curve", ImVec2(600, 0)))
//			{
//				int i = 0;
//				for (const EnergyDistributionSet& set : energyDistributionSets)
//				{
//					ImPlot::PlotScatter(("set " + std::to_string(i++)).c_str(), set.info.detuningVelocity.data(), set.info.longitudinalCoolingForce.data(), set.info.detuningVelocity.size());
//				}
//				ImPlot::EndPlot();
//			}
//
//			if (ImPlot::BeginSubplots("##cf subplots", 1, 3, ImVec2(-1, 400)))
//			{
//				if (ImPlot::BeginPlot("Projection X"))
//				{
//					ImPlot::PlotLine("##", cfData.xAxis.data(), cfData.forceZProjectionX.data(), cfData.xAxis.size());
//					ImPlot::EndPlot();
//				}
//
//				if (ImPlot::BeginPlot("Projection Y"))
//				{
//					ImPlot::PlotLine("##", cfData.yAxis.data(), cfData.forceZProjectionY.data(), cfData.yAxis.size());
//					ImPlot::EndPlot();
//				}
//
//				if (ImPlot::BeginPlot("Projection Z"))
//				{
//					ImPlot::PlotLine("##", cfData.zAxis.data(), cfData.forceZProjectionZ.data(), cfData.zAxis.size());
//					ImPlot::EndPlot();
//				}
//
//				ImPlot::EndSubplots();
//			}
//			cfData.forceXSlice.Plot("x component");
//			ImGui::SameLine();
//			cfData.forceYSlice.Plot("y component");
//			ImGui::SameLine();
//			cfData.forceZSlice.Plot("z component");
//
//		}
//		ImGui::End();
//	}
//}
//
//void EnergyDistributionManager::GenerateCoolingCurveFromFile(std::filesystem::path file)
//{
//	// get all necessary modules
//	FileUtils FileUtils = FileUtils::GetInstance();
//	std::filesystem::path folder = file.parent_path().filename();
//	std::cout << folder << std::endl;
//
//	int end = endIndex;
//	int start = startIndex;
//	if (doAll)
//	{
//		start = 1;
//		end = maxIndex;
//	}
//
//	for (int index = start; index <= end; index++)
//	{
//		// get 3 parameters: U drift tube, electron current, center E lab if index is in file
//		std::array<float, 3> additionalParameter = FileUtils.GetParamtersFromDescriptionFileAtIndex(file, index);
//
//		// if they are not found the index is not in the file
//		if (!additionalParameter[0]) continue;
//
//		// set read electron current and center lab energy
//		activeDist.labEnergiesParameter.driftTubeVoltage = additionalParameter[0];
//		activeDist.eBeamParameter.electronCurrent = additionalParameter[1];
//		activeDist.labEnergiesParameter.centerLabEnergy = additionalParameter[2];
//		activeDist.eBeamParameter.longitudinal_kT_estimate = eBeam->GetLongitudinal_kT(additionalParameter[2]);
//		eBeam->CalculateDetuningEnergy();
//		eBeam->CalculateDetuningVelocity();
//
//		// full procedure to generate one energy distribution 
//		// 1. setup necessary distributions
//		std::filesystem::path densityfile = FileUtils.FindFileWithIndex(file.parent_path() / "e-densities", index);
//		if (densityfile.empty()) continue;
//		eBeam->SetupDistribution(densityfile);
//
//		std::filesystem::path energyfile = FileUtils::GetInstance().FindFileWithIndex(file.parent_path() / "lab-energies", index);
//		if (energyfile.empty()) continue;
//		labEnergies->SetupDistribution(energyfile);
//
//		//mcmc->SetupDistribution();
//
//		// 2. sample from this distribution
//		//mcmc->GenerateSamples();
//
//		// final setup of current distribution
//		activeDist.SetupLabellingThings();
//		//activeDist.SetupBinning(binSettings);
//
//		// 3. generate energy distribution
//		//GenerateEnergyDistribution();
//
//		// calculate cooling force
//		CalculateCoolingForce();
//
//		// store/move and save current distribution that has been worked on
//		PrepareCurrentSet(folder);
//		AddDistributionToSet(std::move(activeDist), currentSetIndex);
//	}
//}
//
//void EnergyDistributionManager::CalculateCoolingForce()
//{
//	std::vector<Point3D> positionSamples = ionBeam->GeneratePositions((int)1e6);
//
//	if (positionSamples.empty())
//	{
//		std::cout << "sampling positions failed\n";
//		return;
//	}
//
//	if (!labEnergies->GetDistribution())
//	{
//		std::cout << "no lab energies were given\n";
//		return;
//	}
//
//	activeDist.cfData.SetupHistograms(ionBeam->GetDistribution());
//
//	for (const Point3D& point : positionSamples)
//	{
//		double x = point.x;
//		double y = point.y;
//		double z = point.z;
//
//		if (activeDist.simplifyParams.cutOutZValues)
//		{
//			if (z < activeDist.simplifyParams.cutOutRange.get().x ||
//				z > activeDist.simplifyParams.cutOutRange.get().y)
//				continue;
//		}
//
//		// calculate velocity magnitude from lab energy given as matrix (TH3D) from outside
//		double labEnergy = labEnergies->Get(x, y, z);
//		double electronVelocityMagnitude = TMath::Sqrt(2 * labEnergy * TMath::Qe() / PhysicalConstants::electronMass);
//
//		// determine direction of velocity based on beam trajectory function
//		TVector3 longitudinalDirection = eBeam->GetDirection(point.z);
//
//		// add random values to velocity in transverse and longitudinal directions:
//		// - calculate longitudinal kT, transverse kT is fixed
//		double long_kT = eBeam->GetLongitudinal_kT(labEnergy);
//		double trans_kT = eBeam->GetTransverse_kT();
//
//		// - use kT to calculate sigmas of gaussians
//		double longSigma = TMath::Sqrt(long_kT * TMath::Qe() / PhysicalConstants::electronMass);
//		double transSigma = TMath::Sqrt(trans_kT * TMath::Qe() / PhysicalConstants::electronMass);
//
//		// - sample from gaussians with these sigmas and add that to the electron velocity
//		longitudinalNormalDistribution = std::normal_distribution<double>(0, longSigma);
//		transverseNormalDistribution = std::normal_distribution<double>(0, transSigma);
//		double longitudinalAddition = longitudinalNormalDistribution(generator);
//
//		TVector3 finalElectronVelocity = electronVelocityMagnitude * longitudinalDirection
//			+ longitudinalAddition * longitudinalDirection;
//
//		double transverseAdditionX = transverseNormalDistribution(generator);
//		double transverseAdditionY = transverseNormalDistribution(generator);
//
//		// we need a vector that is never in line with the longitudinalDirection
//		TVector3 helpVector = TVector3(1, 0, 0);
//		TVector3 transverseDirection1 = longitudinalDirection.Cross(helpVector);
//		TVector3 transverseDirection2 = longitudinalDirection.Cross(transverseDirection1);
//
//		finalElectronVelocity += transverseAdditionX * transverseDirection1;
//		finalElectronVelocity += transverseAdditionY * transverseDirection2;
//
//		double ionVelocityMagnitude = TMath::Sqrt(2 * activeDist.eBeamParameter.coolingEnergy * TMath::Qe() / PhysicalConstants::electronMass); // calc from cooling energy;
//		TVector3 ionVelocityDirection = ionBeam->GetDirection();
//		TVector3 ionVelocity = ionVelocityDirection * ionVelocityMagnitude;
//
//		TVector3 collisionVelocity = ionVelocity - finalElectronVelocity;
//		double collisionVelocityMagnitude = collisionVelocity.Mag();
//
//		double electronDensity = eBeam->GetDensity(point);
//		//std::cout << "density: " << electronDensity << std::endl;
//		int ionCharge = 1;
//		TVector3 coolingforce = cf::CoolingForce(collisionVelocity, kT_trans, electronDensity, ionCharge);
//
//		activeDist.cfData.positionSamples->Fill(x, y, z);
//
//		int bin = activeDist.cfData.forceX->FindBin(x, y, z);
//		activeDist.cfData.forceX->AddBinContent(bin, coolingforce.x());
//		activeDist.cfData.forceY->AddBinContent(bin, coolingforce.y());
//		activeDist.cfData.forceZ->AddBinContent(bin, coolingforce.z());
//	}
//
//	activeDist.cfData.forceX->Divide(activeDist.cfData.positionSamples);
//	activeDist.cfData.forceY->Divide(activeDist.cfData.positionSamples);
//	activeDist.cfData.forceZ->Divide(activeDist.cfData.positionSamples);
//
//	activeDist.cfData.FillData();
//
//	//m_mainCanvas->cd(1);
//	//TH1D* forceProj = activeDist.coolingForceHist->ProjectionZ();
//	//forceProj->Draw();
//	//
//	//m_mainCanvas->cd(4);
//	//TH3D*  forceSmall = (TH3D*)activeDist.coolingForceHist->Rebin3D(s_rebinningFactors[0],
//	//	s_rebinningFactors[1],
//	//	s_rebinningFactors[2], "dist small");
//	//
//	//forceSmall->GetXaxis()->SetTitle("x-axis");
//	//forceSmall->GetYaxis()->SetTitle("y-axis");
//	//forceSmall->GetZaxis()->SetTitle("z-axis");
//	//forceSmall->Draw("BOX2 COLZ");
//	////activeDist.coolingForceHist->Draw("BOX2 COLZ");
//	//
//	//m_mainCanvas->cd(2);
//	//TH1D* posProj = positionHist->ProjectionZ();
//	//posProj->Draw();
//	//
//	//m_mainCanvas->cd(5);
//	//TH3D* positionsnSmall = (TH3D*)positionHist->Rebin3D(s_rebinningFactors[0],
//	//	s_rebinningFactors[1],
//	//	s_rebinningFactors[2], "dist small");
//	//
//	//positionsnSmall->GetXaxis()->SetTitle("x-axis");
//	//positionsnSmall->GetYaxis()->SetTitle("y-axis");
//	//positionsnSmall->GetZaxis()->SetTitle("z-axis");
//	//positionsnSmall->Draw("BOX2 COLZ");
//	////positionHist->Draw("BOX2 COLZ");
//}



