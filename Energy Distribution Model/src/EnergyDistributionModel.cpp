#include "EnergyDistributionModel.h"
#include "PhysicalConstants.h"
#include "FileHandler.h"

#include <math.h>
#include <filesystem>
#include <sstream>

#include <TRootCanvas.h>
#include <TLegend.h>


EnergyDistributionModel::EnergyDistributionModel()
	: Module("Energy Distribution Model")
{
	m_mainCanvas->cd();
	m_mainCanvas->Clear();
	m_mainCanvas->Divide(3,2);
	m_mainCanvas->SetWindowSize(1500, 800);
}

void EnergyDistributionModel::ShowUI()
{
	if (ImGui::BeginChild("##bla2", ImVec2(0, 300.0f)))
	{
		ImGui::Columns(2);

		if (ImGui::Button("Generate Single Energy Distribution"))
		{
			//delete currentDistribution.distribution;
			SetupEnergyDistribution();
			GenerateEnergyDistribution();
			PlotCurrentEnergyDistribution();
			PLotZweightByEnergy();
			FileHandler::GetInstance().SaveEnergyDistributionToFile(currentDistribution);
		}

		if (ImGui::Button("select description file"))
		{
			currentDescriptionFile = FileHandler::GetInstance().OpenFileExplorer();
			maxIndex = FileHandler::GetInstance().GetMaxIndex(currentDescriptionFile);
		}
		ImGui::SameLine();
		ImGui::Text(currentDescriptionFile.filename().string().c_str());

		ImGui::BeginDisabled(currentDescriptionFile.empty());
		if (ImGui::Button("Generate Distributions from description File"))
		{
			GenerateEnergyDistributionsFromFile(currentDescriptionFile);

			PlotEnergyDistributions();
			PLotZweightByEnergy();
			PlotRateCoefficients();
			PlotLongkTDistribution();
			PlotLongVelAddition();
			PlotDistribution();
		}
		ImGui::EndDisabled();

		ImGui::SetNextItemWidth(200.0f);
		ImGui::InputFloat2("energy range", energyRange, "%.1e");

		//ImGui::SameLine();
		ImGui::Checkbox("normalise", &normalise);

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
		ImGui::BeginDisabled(!parameter.cutOutZValues);
		ImGui::InputFloat2("", parameter.cutOutRange);
		ImGui::EndDisabled();
		ImGui::SameLine();
		ImGui::Checkbox("cut out z range", &parameter.cutOutZValues);

		ImGui::NextColumn();
		ShowEnergyDistributionList();
		if (ImGui::Button("clear list"))
		{
			ClearDistributionList();
		}
		ImGui::SameLine();
		if (ImGui::Button("clear plot"))
		{
			for (EnergyDistribution& eDist : energyDistributions)
			{
				eDist.plotted = !eDist.plotted;
			}
		}
		ImGui::EndChild();
	}
	ImGui::Separator();
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

		for (EnergyDistribution& eDist : energyDistributions)
		{
			if (eDist.distribution && eDist.plotted)
			{
				std::string label = eDist.label + Form("##%d", (int)eDist.distribution);
				//std::cout << label << "\n";
				ImPlot::PlotLine(label.c_str(), eDist.binCenters.data(), eDist.binValues.data(), eDist->GetNbinsX());
			}
		}

		ImPlot::EndPlot();
	}
}

void EnergyDistributionModel::ShowPlots()
{
	
}

void EnergyDistributionModel::ShowEnergyDistributionList()
{
	ImGui::Text("loaded distributions");
	if (ImGui::BeginListBox("", ImVec2(-1, 250)))
	{
		for (int i = 0; i < energyDistributions.size(); i++)
		{
			EnergyDistribution& eDist = energyDistributions[i];

			std::string label = eDist.label;
			if (!eDist.tags.empty())
			{
				label += "\n";
				label += eDist.tags;
			}
			label += Form("##%d", (int)eDist.distribution);

			// Render each item as selectable
			if (ImGui::Selectable(label.c_str(), energyDistributions[i].plotted))
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
		}
	}
	ImGui::EndListBox();
}

void EnergyDistributionModel::SetupEnergyDistribution()
{
	MCMC* mcmc = (MCMC*)Module::Get("MCMC");
	ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");
	IonBeam* ionBeam = (IonBeam*)Module::Get("Ion Beam");
	LabEnergies* labEnergies = (LabEnergies*)Module::Get("Lab Energies");

	currentDistribution = EnergyDistribution();

	currentDistribution.mcmcParameter = mcmc->GetParameter();
	currentDistribution.eBeamParameter = eBeam->GetParameter();
	currentDistribution.ionBeamParameter = ionBeam->GetParameter();
	currentDistribution.labEnergiesParameter = labEnergies->GetParameter();
	currentDistribution.eDistParameter = parameter;

	if (!currentDistribution.eBeamParameter.densityFile.empty() && !currentDistribution.eDistParameter.energyFile.empty())
	{
		currentDistribution.folder = currentDistribution.eBeamParameter.densityFile.parent_path().parent_path();
		currentDistribution.index = std::stoi(currentDistribution.eBeamParameter.densityFile.filename().string().substr(0, 4));
	}
	
	float min = std::max(energyRange[0], (float)1e-8);
	float max = energyRange[1];
	int binsPerDecade = 2000;
	float numberBins = (int)((max - min) / 10 * binsPerDecade);
	double factor = TMath::Power((max / min), (1 / numberBins));

	currentDistribution.binCenters.reserve(numberBins);
	currentDistribution.binValues.reserve(numberBins);

	std::vector<double> binEdges;
	binEdges.reserve(numberBins + 1);
	binEdges.push_back(min);
	for (int i = 0; i < numberBins; i++)
	{
		binEdges.push_back(binEdges[i] * factor);
	}

	std::string histDescription = currentDistribution.folder.filename().string() + " " + std::to_string(currentDistribution.index);
	currentDistribution.distribution = new TH1D(("Energy Distribution " + histDescription).c_str(),
		("Energy Distribution " + histDescription).c_str(), numberBins, binEdges.data());
}

void EnergyDistributionModel::GenerateEnergyDistribution()
{
	MCMC* mcmc = (MCMC*)Module::Get("MCMC");
	ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");
	LabEnergies* labEnergies = (LabEnergies*)Module::Get("Lab Energies");

	// sample positions from electron density multiplied with ion density given from outside
	std::vector<Point3D> positionSamples = mcmc->GetSamples();
	currentDistribution.collisionEnergies.reserve(positionSamples.size());

	delete zPositions;
	delete zWeightByEnergy;
	delete long_ktDistribution;
	delete long_VelAddition;
	zPositions = new TH1D("z-positions", "z-positions", 100, 0, 0.65);
	zWeightByEnergy = new TH1D("z weight by energy", "z weight by energy", 100, 0, 0.65);
	long_ktDistribution = new TH1D("long kT", "long kT", 500, 5e-6, 2e-5);
	long_VelAddition = new TH1D("long vel add", "long vel add", 500, -1e4, 1e4);

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

	for (const Point3D& point : positionSamples)
	{
		double x = point.x;
		double y = point.y;
		double z = point.z;
		if (parameter.cutOutZValues)
		{
			if (z < parameter.cutOutRange[0] ||
				z > parameter.cutOutRange[1])
				continue;
		}
		
		// calculate velocity magnitude from lab energy given as matrix (TH3D) from outside
		double labEnergy = labEnergies->Get(x, y, z);
		double electronVelocityMagnitude = TMath::Sqrt(2 * labEnergy * TMath::Qe() / PhysicalConstants::electronMass);

		zPositions->Fill(z);
		zWeightByEnergy->Fill(z, labEnergy);

		// determine direction of velocity based on beam trajectory function
		TVector3 longitudinalDirection = eBeam->GetDirection(point.z);
		//longitudinalDirection.Print();
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
		//std::cout << "kT: " << long_kT << " sigma: " << longSigma << " vel: " << electronVelocityMagnitude << "\n";
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
		currentDistribution.collisionEnergies.push_back(collisionEnergy);

		double sigma = 1 / collisionEnergy;
		currentDistribution.rateCoefficient += sigma * collosionVelocity;

		//std::cout << "position: (" << x << " " << y << " " << z << ")\n";
		////std::cout << "modified: (" << x_modified << " " << y_modified << " " << z_modified << ")\n";
		//std::cout << "lab energy: " << labEnergy << " eV\n";
		//std::cout << "electron velocity: " << electronVelocityMagnitude << " m/s\n";
		//std::cout << "kT trans: " << trans_kT << " eV, kT long: " << long_kT << " eV\n";
		//std::cout << "sigma trans: " << transSigma << " m/s, long sigma: " << longSigma << " m/s\n";
		//std::cout << "collision energy: " << collisionEnergy << " eV\n";
		//std::cout << "\n";
	}
	currentDistribution.rateCoefficient /= positionSamples.size();

	// normalisation
	if (normalise)
	{
		for (int i = 1; i <= currentDistribution->GetNbinsX(); i++)
		{
			currentDistribution->SetBinContent(i, currentDistribution->GetBinContent(i) / currentDistribution->GetBinWidth(i));
			currentDistribution->SetBinError(i, currentDistribution->GetBinError(i) / currentDistribution->GetBinWidth(i));
		}
		if(currentDistribution->Integral())
			currentDistribution->Scale(1.0 / currentDistribution->Integral());
	}

	for (int i = 1; i <= currentDistribution->GetNbinsX(); i++)
	{
		currentDistribution.binCenters.push_back(currentDistribution->GetBinCenter(i));
		currentDistribution.binValues.push_back(currentDistribution->GetBinContent(i));
	}
	
}

void EnergyDistributionModel::GenerateEnergyDistributionsFromFile(std::filesystem::path file)
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
		eBeam->SetCurrent(additionalParameter[1]);
		labEnergies->SetCenterLabEnergy(additionalParameter[2]);
		
		// full procedure to generate one energy distribution 
		// 1. load files if necessary
		std::filesystem::path densityfile = fileHandler.FindFileWithIndex(file.parent_path() / "e-densities", index);
		if (densityfile.empty()) continue;
		eBeam->SetupDistribution(densityfile);

		std::filesystem::path energyfile = FileHandler::GetInstance().FindFileWithIndex(file.parent_path() / "lab-energies", index);
		if (energyfile.empty()) continue;
		labEnergies->SetupDistribution(energyfile);

		// 2. multiply ion and electron beam
		TH3D* result = ionBeam->MultiplyWithElectronDensities(eBeam->GetDistribution());

		// 3. sample from this distribution
		mcmc->SetTargetDistribution(result);
		mcmc->GenerateSamples();

		// 4. prepare current distribution
		parameter.driftTubeVoltage = additionalParameter[0];
		SetupEnergyDistribution();

		if (currentDistribution.eBeamParameter.hasGaussianShape) currentDistribution.tags += "e-gaus ";
		if (currentDistribution.eBeamParameter.hasNoBending) currentDistribution.tags += "no bend ";
		if (currentDistribution.eBeamParameter.hasFixedLongitudinalTemperature) currentDistribution.tags += "fixed kT|| ";
		if (currentDistribution.labEnergiesParameter.useUniformEnergies) currentDistribution.tags += "uniform energy ";
		if (currentDistribution.labEnergiesParameter.useOnlySliceXY) currentDistribution.tags += Form("energy sliced %.3f ", currentDistribution.labEnergiesParameter.sliceToFill);
		if (parameter.cutOutZValues) currentDistribution.tags += Form("z samples %.3f - %.3f", parameter.cutOutRange[0], parameter.cutOutRange[1]);

		currentDistribution.label = Form("%d: U drift = %.2fV", currentDistribution.index, parameter.driftTubeVoltage);

		// 5. generate energy distribution
		GenerateEnergyDistribution();

		// 6. save it to a file
		FileHandler::GetInstance().SaveEnergyDistributionToFile(currentDistribution);

		energyDistributions.push_back(currentDistribution);
	}
}

void EnergyDistributionModel::PlotEnergyDistributions()
{
	m_mainCanvas->cd(3);
	int colors[5] = { kRed, kBlue, kGreen, kOrange, kMagenta };

	gPad->SetLogy();
	gPad->SetLogx();
	
	// Create a legend
	TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

	for (int i = 0; i < energyDistributions.size(); i++)
	{
		if (!energyDistributions[i].distribution) return;
		
		energyDistributions[i]->SetLineColor(colors[i % 5]);
		legend->AddEntry(energyDistributions[i].distribution, energyDistributions[i].label.c_str(), "l");

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

void EnergyDistributionModel::PlotRateCoefficients()
{
	delete rateCoefficients;
	rateCoefficients = new TGraph(energyDistributions.size());
	rateCoefficients->SetTitle("rate Coefficients");

	for (int i = 0; i < energyDistributions.size(); i++)
	{
		double E_d = pow(sqrt(energyDistributions[i].labEnergiesParameter.centerLabEnergy) - sqrt(energyDistributions[i].eBeamParameter.coolingEnergy), 2);
		rateCoefficients->SetPoint(i, E_d, energyDistributions[i].rateCoefficient);
	}

	m_secondCanvas->cd(6);
	rateCoefficients->SetMarkerStyle(21);
	rateCoefficients->GetXaxis()->SetTitle("E_d [eV]");
	rateCoefficients->GetYaxis()->SetTitle("alpha [m^3/s]");
	rateCoefficients->Draw("ALP");
}

void EnergyDistributionModel::PlotCurrentEnergyDistribution()
{
	m_mainCanvas->cd(3);

	currentDistribution->Draw("HIST");

	gPad->SetLogy();
	gPad->SetLogx();
}

void EnergyDistributionModel::PLotZweightByEnergy()
{
	if (!zWeightByEnergy) return;

	m_secondCanvas->cd(4);
	zPositions->Draw();

	m_secondCanvas->cd(5);
	zWeightByEnergy->Divide(zPositions);
	zWeightByEnergy->Draw("Hist");
}

void EnergyDistributionModel::ClearDistributionList()
{
	for (EnergyDistribution& distribution : energyDistributions) 
	{
		delete distribution.distribution;

		//distribution->Delete();
	}
	energyDistributions.clear();
}

void EnergyDistributionModel::PlotLongkTDistribution()
{
	m_secondCanvas->cd(1);

	long_ktDistribution->Draw();
}

void EnergyDistributionModel::PlotLongVelAddition()
{
	m_secondCanvas->cd(2);

	long_VelAddition->Draw();
}

std::string EnergyDistribution::String()
{
	std::string string = Form("# folder: %s\n", folder.filename().string()) + 
						 eDistParameter.String() +
						 eBeamParameter.String() +
						 ionBeamParameter.String() +
						 mcmcParameter.String();

	return string;
}

std::string EnergyDistribution::Filename()
{
	std::ostringstream indexSS;
	std::ostringstream eCoolSS;
	eCoolSS << std::fixed << std::setprecision(3) << eBeamParameter.coolingEnergy;
	indexSS << std::setw(4) << std::setfill('0') << index;

	std::string string = indexSS.str() + "_" + folder.filename().string() +
		"_energyDist_IonBeam" + std::to_string(int(ionBeamParameter.radius * 1000)) + "mm"
		+ "_Ecool" + eCoolSS.str() + "eV" + ".asc";

	return string;
}

std::string EnergyDistributionParameters::String()
{
	std::string string = std::string(Form("# energy distribution parameter:\n"));

	if (cutOutZValues)
	{
		string += std::string(Form("# cut out z sample values between: %f - %f\n", cutOutRange[0], cutOutRange[1]));
	}
						 
	string += std::string(Form("# drift tube voltage: %e V\n", driftTubeVoltage));

	return string;
}
