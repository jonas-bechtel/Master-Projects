#include "pch.h"

#include "EnergyDistributionManager.h"
#include "PhysicalConstants.h"
#include "FileHandler.h"

EnergyDistributionManager::EnergyDistributionManager()
	: Module("Energy Distribution Manager")
{
	m_mainCanvas->cd();
	m_mainCanvas->Clear();
	m_mainCanvas->Divide(3,2);
	m_mainCanvas->SetWindowSize(1500, 800);

	energyDistributions.reserve(5);
}

float* EnergyDistributionManager::GetEnergyRange()
{
	return energyRange;
}

int EnergyDistributionManager::GetBinsPerDecade()
{
	return binsPerDecade;
}

EnergyDistributionParameters& EnergyDistributionManager::GetParameter()
{
	return parameter;
}

std::vector<EnergyDistribution*>& EnergyDistributionManager::GetEnergyDistributions()
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

		ImGui::BeginDisabled(parameter.limitBinSize);
		ImGui::SetNextItemWidth(150.0f);
		ImGui::InputFloat2("energy range", energyRange, "%.1e");
		ImGui::SameLine();
		ImGui::SetNextItemWidth(100.0f);
		ImGui::InputInt("bins per decade", &binsPerDecade);
		ImGui::EndDisabled();

		ImGui::Separator();
		ImGui::SetNextItemWidth(200.0f);
		ImGui::Checkbox("limit lower bin size", parameter.limitBinSize);
		ImGui::SameLine();
		ImGui::SetNextItemWidth(100.0f);
		ImGui::InputDouble("min bin size", parameter.minBinSize, 0, 0, "%.1e");

		ImGui::SetNextItemWidth(200.0f);
		ImGui::BeginDisabled(!parameter.cutOutZValues);
		ImGui::InputFloat2("", parameter.cutOutRange);
		ImGui::EndDisabled();
		ImGui::SameLine();
		ImGui::Checkbox("cut out z range", parameter.cutOutZValues);

		ImGui::Separator();
		ImGui::PushID("analytical energy distribution");
		if (ImGui::Button("generate analytical distribution"))
		{
			GenerateAnalyticalDistribution();
		}
		ImGui::SetNextItemWidth(200.0f);
		ImGui::InputFloat2("energy range", analyticalEnergyRange, "%.1e");
		ImGui::SameLine();
		ImGui::SetNextItemWidth(100.0f);
		ImGui::InputInt("number bins", &analyticalNumberBins);
		ImGui::SetNextItemWidth(100.0f);
		ImGui::InputDouble("detuning energy [eV]", analyticalParameter.detuningEnergy);
		ImGui::SetNextItemWidth(100.0f);
		ImGui::InputDouble("transverse kT [eV]", analyticalParameter.transverseTemperature);
		ImGui::SetNextItemWidth(100.0f);
		ImGui::InputDouble("longitudinal kT [eV]", analyticalParameter.longitudinalTemperature);
		ImGui::PopID();

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
				EnergyDistribution* eDist = energyDistributions[i];

				std::string label = eDist->label;
				if (!eDist->tags.empty())
				{
					label += "\n";
					label += eDist->tags;
				}

				// Render each item as selectable
				if (ImGui::Selectable(label.c_str(), eDist->plotted, ImGuiSelectableFlags_AllowItemOverlap))
				{
					eDist->plotted = !eDist->plotted;
				}

				if (ImGui::BeginItemTooltip())
				{
					ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
					ImGui::TextUnformatted(eDist->String().c_str());
					ImGui::PopTextWrapPos();
					ImGui::EndTooltip();
				}

				if (!eDist->isAnalytical)
				{
					ImGui::SameLine();
					ImGui::Checkbox("normalised", &eDist->showNormalisedByWidth);
				}
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
					EnergyDistribution* energyDist = FileHandler::GetInstance().LoadEnergyDistribution(filename, loadSamples);
					AddDistributionToList(energyDist);
				}
			}
		}
		ImGui::SameLine();
		ImGui::Checkbox("load samples", &loadSamples);
		//ImGui::SameLine();
		//if (ImGui::Button("load samples"))
		//{
		//	std::vector<std::filesystem::path> filenames = FileHandler::GetInstance().SelectFiles("output\\");
		//	if (!filenames.empty())
		//	{
		//		for (const auto& filename : filenames)
		//		{
		//			EnergyDistribution* energyDist = FileHandler::GetInstance().LoadEnergyDistributionSamples(filename);
		//			AddDistributionToList(energyDist);
		//		}
		//	}
		//}

		if (ImGui::Button("clear list"))
		{
			ClearDistributionList();
		}
		ImGui::SameLine();
		if (ImGui::Button("clear plot"))
		{
			for (EnergyDistribution* eDist : energyDistributions)
			{
				eDist->plotted = false;
			}
		}
		ImGui::SameLine();
		if (ImGui::Button("normalise all"))
		{
			for (EnergyDistribution* eDist : energyDistributions)
			{
				if (eDist->isAnalytical) continue;
				eDist->showNormalisedByWidth = true;
			}
		}
		ImGui::SameLine();
		if (ImGui::Button("unnormalise all"))
		{
			for (EnergyDistribution* eDist : energyDistributions)
			{
				if (eDist->isAnalytical) continue;
				eDist->showNormalisedByWidth = false;
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

	if (ImPlot::BeginPlot("collision Energy distribution"))
	{
		ImPlot::SetupAxis(ImAxis_X1, "collision energy [eV]");
		ImPlot::SetupAxis(ImAxis_Y1, "f(E)", ImPlotAxisFlags_AutoFit);
		if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
		if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
		ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);
		int i = 0;
		for (EnergyDistribution* eDist : energyDistributions)
		{
			if (eDist->plotted)
			{
				ImGui::PushID(i);

				if (eDist->showNormalisedByWidth || eDist->isAnalytical)
				{
					// Get the automatic color for this pair
					ImVec4 color = ImPlot::GetColormapColor(i % ImPlot::GetColormapSize());

					// Plot the first line with the automatic color
					ImPlot::PushStyleColor(ImPlotCol_Line, color);
					ImPlot::PlotLine(eDist->label.c_str(), eDist->binCenters.data(), eDist->binValuesNormalised.data(), eDist->binValuesNormalised.size());
					ImPlot::PopStyleColor();
					color.x *= 2;
					color.y *= 2;
					color.z *= 2;

					// Plot the second line with a lighter color and dashed
					ImPlot::PushStyleColor(ImPlotCol_Line, color);
					ImPlot::PlotLine("##", eDist->fitX.data(), eDist->fitY.data(), eDist->fitX.size(), ImPlotLineFlags_Segments);
					ImPlot::PopStyleColor();
				}
				else
				{
					ImPlot::PlotLine(eDist->label.c_str(), eDist->binCenters.data(), eDist->binValues.data(), eDist->binCenters.size());
				}
				i++;
			}
		}
		ImPlot::EndPlot();
	}
}

void EnergyDistributionManager::AddDistributionToList(EnergyDistribution* distribution)
{
	energyDistributions.push_back(distribution);
	EnergyDistribution::s_allDistributions[distribution->eBeamParameter.detuningEnergy] = distribution;
}

void EnergyDistributionManager::RemoveDistributionFromList(int index)
{
	delete energyDistributions[index];
	energyDistributions.erase(energyDistributions.begin() + index);
}

void EnergyDistributionManager::GenerateEnergyDistribution()
{
	currentDistribution = new EnergyDistribution();
	currentDistribution->SetupFromCurrentEnvironment();
	EnergyDistributionParameters& parameter = currentDistribution->eDistParameter;

	MCMC* mcmc = (MCMC*)Module::Get("MCMC");
	ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");
	LabEnergies* labEnergies = (LabEnergies*)Module::Get("Lab Energies");

	// sample positions from electron density multiplied with ion density given from outside
	std::vector<Point3D> positionSamples = mcmc->GetSamples();
	currentDistribution->collisionEnergies.reserve(positionSamples.size());

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

	double kTLongGuess = eBeam->GetLongitudinal_kT(currentDistribution->labEnergiesParameter.centerLabEnergy);
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

	for (const Point3D& point : positionSamples)
	{
		double x = point.x;
		double y = point.y;
		double z = point.z;
		if (parameter.cutOutZValues)
		{
			if (z < parameter.cutOutRange.get().x ||
				z > parameter.cutOutRange.get().y)
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
		double ionVelocityMagnitude = TMath::Sqrt(2 * eBeam->GetParameter().coolingEnergy * TMath::Qe() / PhysicalConstants::electronMass); // calc from cooling energy;
		TVector3 ionVelocity(0, 0, ionVelocityMagnitude);

		double collosionVelocity = (finalElectronVelocity - ionVelocity).Mag();

		// calculate collision energy [eV] and put it in a histogram
		double collisionEnergy = 0.5 * PhysicalConstants::electronMass * collosionVelocity * collosionVelocity / TMath::Qe();
		currentDistribution->Fill(collisionEnergy);
		currentDistribution->collisionEnergies.push_back(collisionEnergy);
	}

	// write non normalised values to the vector
	for (int i = 1; i <= currentDistribution->GetNbinsX(); i++)
	{
		currentDistribution->binCenters.push_back(currentDistribution->GetBinCenter(i));
		currentDistribution->binValues.push_back(currentDistribution->GetBinContent(i));
	}

	// normalise the distribution to the width and then to one
	for (int i = 1; i <= currentDistribution->GetNbinsX(); i++)
	{
		currentDistribution->SetBinContent(i,
			currentDistribution->GetBinContent(i) / currentDistribution->GetBinWidth(i));
	}
	if (currentDistribution->Integral("width"))
		currentDistribution->Scale(1.0 / currentDistribution->Integral("width"));

	// write normalised values to the vector
	for (int i = 1; i <= currentDistribution->GetNbinsX(); i++)
	{
		currentDistribution->binValuesNormalised.push_back(currentDistribution->GetBinContent(i));
	}
	
	currentDistribution->RemoveEdgeZeros();
	FitAnalyticalToGeneratedDistribution(currentDistribution);

	// store and save current distribution that has been worked on
	AddDistributionToList(currentDistribution);
	
	std::cout << "Ed1: " << currentDistribution->eBeamParameter.detuningEnergy << "\n";

	if (saveAsHist)
	{
		FileHandler::GetInstance().SaveEnergyDistributionHistToFile(currentDistribution);
	}
	
	if (saveSamplesToFile)
	{
		FileHandler::GetInstance().SaveEnergyDistributionSamplesToFile(currentDistribution);
	}

	// create new distribution object that will be worked on now
	//currentDistribution = new EnergyDistribution();
}

void EnergyDistributionManager::GenerateEnergyDistributionsFromFile(std::filesystem::path file)
{
	// get all necessary modules
	FileHandler fileHandler = FileHandler::GetInstance();
	ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");
	IonBeam* ionBeam = (IonBeam*)Module::Get("Ion Beam");
	MCMC* mcmc = (MCMC*)Module::Get("MCMC");
	LabEnergies* labEnergies = (LabEnergies*)Module::Get("Lab Energies");

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
		//currentDistribution->eDistParameter.driftTubeVoltage = additionalParameter[0];
		labEnergies->GetParameter().driftTubeVoltage = additionalParameter[0];
		eBeam->GetParameter().electronCurrent = additionalParameter[1];
		labEnergies->GetParameter().centerLabEnergy = additionalParameter[2];
		eBeam->GetParameter().longitudinal_kT = eBeam->GetLongitudinal_kT(additionalParameter[2]);
		
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
		if (!energyDistributions[i]) return;
		
		energyDistributions[i]->SetLineColor(colors[i % 5]);
		legend->AddEntry(energyDistributions[i], energyDistributions[i]->label.c_str(), "l");

		if (i == 0)
		{
			energyDistributions[i]->Draw("HIST");
		}
		else
		{
			energyDistributions[i]->Draw("HIST SAME");
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

void EnergyDistributionManager::FitAnalyticalToGeneratedDistribution(EnergyDistribution* distribution)
{
	double detuningEnergy = distribution->eBeamParameter.detuningEnergy;
	double kt_trans = distribution->eBeamParameter.transverse_kT;
	double kT_long = distribution->eBeamParameter.longitudinal_kT;

	double peakWidthGuess = sqrt(pow((kt_trans * log(2)), 2) + 16 * log(2) * kT_long * detuningEnergy);;
	double energyMin = distribution->eBeamParameter.detuningEnergy + peakWidthGuess;
	double energyMax = distribution->eBeamParameter.detuningEnergy - peakWidthGuess;
	const int nParameter = 4;
	// first parameter is a scaling factor
	double initialGuess[nParameter] = {1, detuningEnergy, kt_trans, kT_long };

	TF1* fitFunction = new TF1("fit function", this, &EnergyDistributionManager::AnalyticalEnergyDistributionFit, energyMin, energyMax, nParameter);
	fitFunction->SetParameters(initialGuess); 

	TFitResultPtr result = distribution->Fit(fitFunction, "RN");

	double* fitParameter = fitFunction->GetParameters();
	distribution->analyticalParameter.detuningEnergy = fitParameter[1];
	distribution->analyticalParameter.transverseTemperature = fitParameter[2];
	distribution->analyticalParameter.longitudinalTemperature = fitParameter[3];

	// Fill distribution Fit data with these parameters
	double energyStep = (energyMax - energyMin) / distribution->fitX.size();
	for (int i = 0; i < distribution->fitX.size(); i++)
	{
		double energy = energyMin + i * energyStep;
		distribution->fitX[i] = energy;
		distribution->fitY[i] = AnalyticalEnergyDistributionFit(&energy, fitParameter);
	}
	fitFunction->Delete();
}

void EnergyDistributionManager::GenerateAnalyticalDistribution()
{
	EnergyDistribution* analyticalDist = new EnergyDistribution();
	analyticalDist->analyticalParameter = analyticalParameter;
	analyticalDist->label = Form("E_d: %.3f analytical distribution", analyticalParameter.detuningEnergy.get());
	analyticalDist->isAnalytical = true;

	double eMin = analyticalEnergyRange[0];
	double eMax = analyticalEnergyRange[1];
	double step = (eMax - eMin) / analyticalNumberBins;

	analyticalDist->binCenters.reserve(analyticalNumberBins);
	analyticalDist->binValuesNormalised.reserve(analyticalNumberBins);

	for (int i = 0; i < analyticalNumberBins; i++)
	{
		double energy = eMin + i * step;
		double value = AnalyticalEnergyDistribution(energy, analyticalParameter.detuningEnergy,
			analyticalParameter.transverseTemperature, analyticalParameter.longitudinalTemperature);

		// skip too small values because they ruin the plot
		if (value < 1e-8) continue;

		analyticalDist->binCenters.push_back(energy);
		analyticalDist->binValuesNormalised.push_back(value);
	}
	AddDistributionToList(analyticalDist);
}

double EnergyDistributionManager::AnalyticalEnergyDistributionFit(double* x, double* params)
{
	double factor = params[0];
	return factor * AnalyticalEnergyDistribution(x, &params[1]);
}

double EnergyDistributionManager::ComplexErrorFunction(double* x, double* par)
{
	// Special implementation of the complex errorfunction:
	// Returns real part of complex errorfunction (Voigt-Faddeeva function) of the argument representing complex part of the parameter
	//
	// x[0] - complex part of the parameter
	Double_t yy = x[0];
	//=============================
	Double_t factor1 = 2.82842712474619029; // 2.*TMath::Sqrt(2.)
	Double_t factor2 = 2.50662827463100024; // TMath::Sqrt(TMath::TwoPi())
	Double_t f = factor2 * TMath::Voigt(0., 1., factor1 * yy, 5);
	return f;
}

double EnergyDistributionManager::DawsonIntegral(double* xarg, double* par)
{
	// Returns Dawson's integral for any real x.
	//        (From Numerical Recipes ch. 6.10)
	Int_t i, n0;
	#define NMAX 6        // #define NMAX 6 
	const Double_t H = 0.4;      // #define H 0.4 
	const Double_t A1 = 2.0 / 3.0; // #define A1 (2.0/3.0) 
	const Double_t A2 = 0.4;   // #define A2 0.4 
	const Double_t A3 = 2.0 / 7.0; // #define A3 (2.0/7.0) 
	Double_t d1, d2, e1, e2, sum, x2, xp, xx, ans;
	Double_t x = xarg[0];
	static Double_t c[NMAX + 1];
	static Int_t init = 0; //  Flag is 0 if we need to initialize, else 1.

	if (init == 0)
	{
		init = 1;
		for (i = 1; i <= NMAX; i++)
			c[i] = TMath::Exp(-H * H * (Double_t)((2 * i - 1) * (2 * i - 1)));
	}
	if (TMath::Abs(x) < 0.2)
	{   // Use series expansion. 
		x2 = x * x;
		ans = x * (1.0 - A1 * x2 * (1.0 - A2 * x2 * (1.0 - A3 * x2)));
	}
	else
	{
		//Use sampling theorem representation. 
		xx = TMath::Abs(x);
		//    n0=2*(Int_t)(0.5*xx/H+0.5); 
		n0 = 2 * TMath::Nint(0.5 * xx / H);
		xp = xx - (Double_t)n0 * H;
		e1 = TMath::Exp(2.0 * xp * H);
		e2 = e1 * e1;
		d1 = (Double_t)n0 + 1.0;
		d2 = d1 - 2.0;
		sum = 0.0;
		//    for (i=1;i<=NMAX;i++,d1+=2.0,d2-=2.0,e1*=e2) {
		for (i = 1; i <= NMAX; i++, d1 += 2.0, d2 -= 2.0, e1 *= e2)
		{
			sum += c[i] * (e1 / d1 + 1.0 / (d2 * e1));
		}

		ans = 0.5641895835 * TMath::Sign(TMath::Exp(-xp * xp), x) * sum;
		//Constant is 1/sqrt(pi) 
	}
	return ans;
}

double EnergyDistributionManager::ExpDiff(double* x, double* par)
{
	Double_t yy = x[0];
	const Double_t A = 1. / 3.;
	//=============================
	Double_t f;
	if (TMath::Abs(yy) > 1e-6) 
	{
		f = (TMath::Exp(yy) - TMath::Exp(-yy)) / yy;
	}
	else 
	{
		f = 2. + A * yy * yy;
	}
	return f;
}

double EnergyDistributionManager::AnalyticalEnergyDistribution(double* x, double* par)
{
	// Flaterned energy distribution, expects transverse mean velocity = 0!
	// Original version written by Andreas Wolf (edistr.C) returns the same result normalised by 1/sqrt(ee)!!!
	//
	// x[0] energy in CM (the variable of the distribution) [eV]
	// par[0] detuning energy (average CM energy) [eV]
	// par[1] transversal electron temperature [eV]
	// par[2] longitudinal electron temperature [eV]

	Double_t ee = x[0];			//energy in CM (variable of the distribution) [eV]
	Double_t edet = par[0]; 	//detuning energy (average CM energy) [eV]
	Double_t tperp = par[1];	//transversal temperature [eV]
	Double_t tpar = par[2];		//longitudinal temperature [eV]
	//=============================
	Double_t rpar[1];
	rpar[0] = 0.;
	Double_t rx[1];
	//=============================
	const Double_t sqrtpi = 1.77245385090551588;  // TMath::Sqrt(TMath::Pi())
	const Double_t eps = 1e-5;
	//cout<<"edis: "<<ee<<" "<<edet<<" "<<tperp<<" "<<tpar<<" "<<endl;
	ee = TMath::Max(ee, 1e-20);
	tperp = TMath::Max(tperp, 1e-20);
	tpar = TMath::Max(tpar, 1e-20);
	//=============================
	Double_t eta = 1.0 - tpar / tperp;
	Double_t ksi = TMath::Abs(eta);

	if (ksi > eps)
	{
		Double_t aa = TMath::Sqrt(ksi * ee / tpar);
		Double_t bb = TMath::Sqrt(edet / ksi / tpar);
		//=============================
		Double_t exp1 = TMath::Sqrt(ee) - TMath::Sqrt(edet);
		exp1 = exp1 * exp1 / tpar;
		Double_t exp2 = TMath::Sqrt(ee) + TMath::Sqrt(edet);
		exp2 = exp2 * exp2 / tpar;
		//=============================
		if (eta > 0)
		{
			// Tpar < Tperp
			Double_t f = TMath::Exp(-exp1);
			if (bb >= aa)
			{
				rx[0] = bb - aa; f = f * ComplexErrorFunction(rx, rpar);
			}
			else
			{
				rx[0] = aa - bb; Double_t cc = (ee - edet / ksi) / tperp;
				f = 2. * TMath::Exp(-cc) - f * ComplexErrorFunction(rx, rpar);
			}
			rx[0] = bb + aa; f = f - TMath::Exp(-exp2) * ComplexErrorFunction(rx, rpar);
			f = f / (2. * tperp * TMath::Sqrt(tpar) * aa);  //factor for f(E)/sqrt(E)

			return f * sqrt(ee);
		}
		else
		{
			// Tpar > Tperp
			Double_t f = TMath::Exp(-exp1);
			rx[0] = bb + aa; f = f * DawsonIntegral(rx, rpar);
			rx[0] = bb - aa; f = f - TMath::Exp(-exp2) * DawsonIntegral(rx, rpar);
			f = f / (sqrtpi * tperp * TMath::Sqrt(tpar) * aa);  //factor for f(E)/sqrt(E)

			return f * sqrt(ee);
		}
	}
	else
	{
		// Tpar = Tperp
		Double_t aa = 2.0 * TMath::Sqrt(ee * edet) / tperp;
		Double_t f = 1 / (sqrtpi * tperp * TMath::Sqrt(tperp));  //factor for f(E)/sqrt(E)

		if (aa > 1)
		{
			Double_t exp1 = (TMath::Sqrt(ee) - TMath::Sqrt(edet));
			exp1 = exp1 * exp1 / tperp;
			Double_t exp2 = (TMath::Sqrt(ee) + TMath::Sqrt(edet));
			exp2 = exp2 * exp2 / tperp;
			f = f * (TMath::Exp(-exp1) - TMath::Exp(-exp2)) / aa;
		}
		else
		{
			rx[0] = aa;
			f = f * TMath::Exp(-(ee + edet) / tperp) * ExpDiff(rx, rpar);
		}
		return f * sqrt(ee);
	}
}

double EnergyDistributionManager::AnalyticalEnergyDistribution(double Ecm, double Ed, double Ttr, double Tlong)
{
	// The same as previous, different parameter representation
	//
	// Ecm energy in CM (the variable of the distribution) [eV]
	// Ed detuning energy (average CM energy) [eV]
	// Ttr transversal electron temperature [eV]
	// Tlong longitudinal electron temperature [eV]

	double x[1];
	x[0] = Ecm;

	double par[3];
	par[0] = Ed;
	par[1] = Ttr;
	par[2] = Tlong;

	return AnalyticalEnergyDistribution(x, par);
}

void EnergyDistributionManager::PlotLongVelAddition()
{
	if (!long_VelAddition) return;

	m_secondCanvas->cd(2);

	long_VelAddition->Draw();
}

