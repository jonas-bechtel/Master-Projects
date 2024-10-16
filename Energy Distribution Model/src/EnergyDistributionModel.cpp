#include "EnergyDistributionModel.h"
#include "PhysicalConstants.h"
#include "FileHandler.h"

#include "imgui.h"

#include <math.h>
#include <filesystem>
#include <sstream>

#include <TRootCanvas.h>


EnergyDistributionModel::EnergyDistributionModel()
	: Module("Energy Distribution Model")
{
}

void EnergyDistributionModel::ShowUI()
{
	if (ImGui::Button("Load lab energies"))
	{
		std::filesystem::path file = FileHandler::GetInstance().OpenFileExplorer();
		LoadLabEnergyFile(file);
	}

	ImGui::InputFloat2("energy range", energyRange, "%.1e");

	if (ImGui::Button("Generate Single Energy Distribution"))
	{
		//delete currentDistribution.distribution;
		SetupEnergyDistribution();
		GenerateEnergyDistribution();
		PlotCurrentEnergyDistribution();
		FileHandler::GetInstance().SaveEnergyDistributionToFile(currentDistribution);
	}
	ImGui::SameLine();
	ImGui::Checkbox("normalise", &normalise);

	if (ImGui::Button("Generate Distributions from description File"))
	{
		std::filesystem::path descriptionFile = FileHandler::GetInstance().OpenFileExplorer();
		GenerateEnergyDistributionsFromFile(descriptionFile);
	}
}

void EnergyDistributionModel::LoadLabEnergyFile(std::filesystem::path file)
{
	if (!file.empty())
	{
		delete m_distribution;
		m_distribution = FileHandler::GetInstance().LoadMatrixFile(file);
		m_distribution->SetTitle("lab energies");
		m_distribution->SetName("lab energies");
		
		loadedEnergyFile = file;

		PlotDistribution();
		PlotLabEnergyProjections();
	}
}

void EnergyDistributionModel::SetupEnergyDistribution()
{
	MCMC* mcmc = (MCMC*)Module::Get("MCMC");
	ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");
	IonBeam* ionBeam = (IonBeam*)Module::Get("Ion Beam");

	currentDistribution = EnergyDistribution();

	currentDistribution.mcmcParameter = mcmc->GetParameter();
	currentDistribution.eBeamParameter = eBeam->GetParameter();
	currentDistribution.ionBeamParameter = ionBeam->GetParameter();

	currentDistribution.densityFile = eBeam->GetLoadedDensityFile();
	currentDistribution.energyFile = loadedEnergyFile;

	currentDistribution.folder = currentDistribution.densityFile.parent_path().parent_path();
	currentDistribution.index = std::stoi(currentDistribution.densityFile.filename().string().substr(0, 4));

	std::ostringstream indexSS;
	std::ostringstream eCoolSS;
	eCoolSS << std::fixed << std::setprecision(3) << currentDistribution.eBeamParameter.coolingEnergy;
	indexSS << std::setw(4) << std::setfill('0') << currentDistribution.index;

	// name will be used as the filename
	currentDistribution.name = indexSS.str() + "_" + currentDistribution.folder.filename().string() +
							"energyDist_IonBeam" + int(currentDistribution.ionBeamParameter.radius * 1000) + "mm"
							+ "_Ecool" + eCoolSS.str() + "eV";

	float min = std::max(energyRange[0], (float)1e-8);
	float max = energyRange[1];
	int binsPerDecade = 2000;
	float numberBins = (int)((max - min) / 10 * binsPerDecade);
	double factor = TMath::Power((max / min), (1 / numberBins));
	//std::cout << "N: " << numberBins << " factor: " << factor << "\n";

	std::vector<double> binEdges;
	binEdges.reserve(numberBins + 1);
	binEdges.push_back(min);
	for (int i = 0; i < numberBins; i++)
	{
		//std::cout << binEdges[i] * factor << "\n";
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

	// sample positions from electron density multiplied with ion density given from outside
	std::vector<Point3D> positionSamples = mcmc->GetSamples();

	if (positionSamples.empty())
	{
		std::cout << "no sampled positions were given\n";
		return; 
	}

	if (!m_distribution)
	{
		std::cout << "no lab energies were given\n";
		return; 
	}
	//TH1D* energy = new TH1D("bla", "bla", 100, 0, 0.6);

	for (const Point3D& point : positionSamples)
	{
		double x = point.x;
		double y = point.y;
		double z = point.z;

		int numberBinsX = m_distribution->GetXaxis()->GetNbins();
		int numberBinsY = m_distribution->GetYaxis()->GetNbins();
		int numberBinsZ = m_distribution->GetZaxis()->GetNbins();

		double x_modified = std::min(std::max(x, m_distribution->GetXaxis()->GetBinCenter(1)), m_distribution->GetXaxis()->GetBinCenter(numberBinsX) - 1e-4);
		double y_modified = std::min(std::max(y, m_distribution->GetYaxis()->GetBinCenter(1)), m_distribution->GetYaxis()->GetBinCenter(numberBinsY) - 1e-4);
		double z_modified = std::min(std::max(z, m_distribution->GetZaxis()->GetBinCenter(1)), m_distribution->GetZaxis()->GetBinCenter(numberBinsZ) - 1e-4);

		// calculate velocity magnitude from lab energy given as matrix (TH3D) from outside
		double labEnergy = m_distribution->Interpolate(x_modified, y_modified, z_modified);
		double electronVelocityMagnitude = TMath::Sqrt(2 * labEnergy * TMath::Qe() / PhysicalConstants::electronMass);

		//energy->Fill(z_modified, labEnergy);

		// determine direction of velocity based on beam trajectory function
		TVector3 longitudinalDirection = eBeam->GetDirection(point.z);
		//longitudinalDirection.Print();
		TVector3 transverseDirection = longitudinalDirection.Orthogonal();

		// add random values to velocity in transverse and longitudinal directions:
		// - calculate longitudinal kT, transverse kT is fixed
		double long_kT = eBeam->GetLongitudinal_kT(labEnergy);
		double trans_kT = eBeam->GetTransverse_kT();

		// - use kT to calculate sigmas of gaussians
		double longSigma = TMath::Sqrt(long_kT * TMath::Qe() / PhysicalConstants::electronMass);
		double transSigma = TMath::Sqrt(trans_kT * TMath::Qe() / PhysicalConstants::electronMass);

		// - sample from gaussians with these sigmas and add that to the electron velocity
		longitudinalNormalDistribution = std::normal_distribution<double>(0, longSigma);
		transverseNormalDistribution = std::normal_distribution<double>(0, transSigma);

		double longitudinalAddition = longitudinalNormalDistribution(generator);
		double transverseAddition = transverseNormalDistribution(generator);
		double transverseAdditionAngle = angleDistribution(generator);

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

		//std::cout << "position: (" << x << " " << y << " " << z << ")\n";
		////std::cout << "modified: (" << x_modified << " " << y_modified << " " << z_modified << ")\n";
		//std::cout << "lab energy: " << labEnergy << " eV\n";
		//std::cout << "electron velocity: " << electronVelocityMagnitude << " m/s\n";
		//std::cout << "kT trans: " << trans_kT << " eV, kT long: " << long_kT << " eV\n";
		//std::cout << "sigma trans: " << transSigma << " m/s, long sigma: " << longSigma << " m/s\n";
		//std::cout << "collision energy: " << collisionEnergy << " eV\n";
		//std::cout << "\n";
	}

	// normalisation
	if (normalise)
	{
		for (int i = 1; i <= currentDistribution->GetNbinsX(); i++)
		{
			currentDistribution->SetBinContent(i, currentDistribution->GetBinContent(i) / currentDistribution->GetBinWidth(i));
			currentDistribution->SetBinError(i, currentDistribution->GetBinError(i) / currentDistribution->GetBinWidth(i));
		}
		currentDistribution->Scale(1.0 / currentDistribution->Integral());
	}
	
	//energy->Divide(samplesInZ);
	//m_secondCanvas->cd(4);
	//energy->Draw("Hist");
}

void EnergyDistributionModel::GenerateEnergyDistributionsFromFile(std::filesystem::path file)
{
	ClearDistributionList();

	// get all necessary modules
	FileHandler fileHandler = FileHandler::GetInstance();
	ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");
	IonBeam* ionBeam = (IonBeam*)Module::Get("Ion Beam");
	MCMC* mcmc = (MCMC*)Module::Get("MCMC");

	for (int index = 1; true; index++)
	{
		if (index > 5) break;

		// get 3 parameters: U drift tube, electron current, center E lab
		float* additionalParameter = fileHandler.GetParamtersFromDescriptionFileAtIndex(file, index);

		// if they are not found the index is not in the file
		if (!additionalParameter) break;
		
		// full procedure to generate one energy distribution 
		// 1. load files
		std::filesystem::path densityfile = fileHandler.FindFileWithIndex(file.parent_path() / "e-densities", index);
		std::filesystem::path energyfile = fileHandler.FindFileWithIndex(file.parent_path() / "lab-energies", index);
		eBeam->LoadDensityFile(densityfile);
		LoadLabEnergyFile(energyfile);

		// 2. multiply ion and electron beam
		//ionBeam->MultiplyWithElectronDensities(eBeam->GetDistribution());

		// 3. sample from this distribution
		mcmc->GenerateSamples();

		SetupEnergyDistribution();
		currentDistribution.driftTubeVoltage = additionalParameter[0];
		currentDistribution.centerLabEnergy = additionalParameter[2];

		// 4. generate energy distribution
		GenerateEnergyDistribution();

		// 5. save it to a file
		FileHandler::GetInstance().SaveEnergyDistributionToFile(currentDistribution);

		energyDistributions.push_back(currentDistribution);
	}

	PlotEnergyDistributions();
}

void EnergyDistributionModel::PlotEnergyDistributions()
{
	m_mainCanvas->cd(1);
	int colors[5] = { kRed, kBlue, kGreen, kOrange, kMagenta };

	gPad->SetLogy();
	gPad->SetLogx();
	
	for (int i = 0; i < energyDistributions.size(); i++)
	{
		if (!energyDistributions[i].distribution) return;
		
		energyDistributions[i]->SetLineColor(colors[i % 5]);

		if (i == 0)
		{
			energyDistributions[i]->Draw("HIST");
			energyDistributions[i]->GetXaxis()->SetRangeUser(energyRange[0], energyRange[1]);
		}
		else
		{
			energyDistributions[i]->Draw("HIST SAME");
		}
		std::cout << "bla\n";
	}
}

void EnergyDistributionModel::PlotCurrentEnergyDistribution()
{
	m_mainCanvas->cd(1);

	currentDistribution->Draw("HIST");

	gPad->SetLogy();
	gPad->SetLogx();
}

void EnergyDistributionModel::PlotLabEnergyProjections()
{
	if (!m_distribution) return;

	delete labEnergyProjectionX; 
	delete labEnergyProjectionY;
	delete labEnergyProjectionZ;

	m_secondCanvas->cd(1);
	labEnergyProjectionX = m_distribution->ProjectionX();
	labEnergyProjectionX->Draw();

	m_secondCanvas->cd(2);
	labEnergyProjectionY = m_distribution->ProjectionY();
	labEnergyProjectionY->Draw();

	m_secondCanvas->cd(3);
	labEnergyProjectionZ = m_distribution->ProjectionZ();
	labEnergyProjectionZ->Draw();
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

std::string EnergyDistribution::String()
{
	std::string string = "# electron beam parameter:\n" +
						  eBeamParameter.String() +
						  "# ion beam parameter:\n" +
						  ionBeamParameter.String() +
						  "# mcmc sampling parameter:\n" +
						  mcmcParameter.String();

	string += "# additional parameter:\n";
	string += Form("# drift tube voltage: %e V\n", driftTubeVoltage);
	string += Form("# lab energy in center: %e eV\n", centerLabEnergy);
	string += Form("# density file: %s\n", densityFile.filename().string().c_str());
	string += Form("# energy file: %s\n", energyFile.filename().string().c_str());
	string += Form("# folder: %s\n", folder.filename().string());

	return string;
}
