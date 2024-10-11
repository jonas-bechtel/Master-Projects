#include "EnergyDistributionModel.h"
#include "PhysicalConstants.h"
#include "FileLoader.h"

#include "imgui.h"

#include <math.h>
#include <filesystem>
#include <TRootCanvas.h>


EnergyDistributionModel::EnergyDistributionModel()
	: Module("Energy Distribution Model")
{
}

void EnergyDistributionModel::ShowUI()
{
	mcmcSampler.ShowWindow();
	electronBeam.ShowWindow();
	ionBeam.ShowWindow();

	if (electronBeam.HasDensityChanged() || ionBeam.HasDensityChanged())
	{
		if (electronBeam.GetDistribution() || ionBeam.GetDistribution())
		{
			ionBeam.CreateIonBeam(electronBeam.GetDistribution());
			TH3D* ionElectronDensity = MultiplyElectronAndIonDensities(electronBeam.GetDistribution(), ionBeam.GetDistribution());
			mcmcSampler.SetTargetDistribution(ionElectronDensity);
		}
	}

	if (ImGui::Button("Load lab energies"))
	{
		std::filesystem::path file = FileLoader::getInstance().openFileExplorer("data\\lab-energies\\");
		LoadLabEnergyFile(file);
	}

	ImGui::InputFloat2("energy range", energyRange, "%.1e");

	if (ImGui::Button("Generate Energy Distribution"))
	{
		GenerateEnergyDistribution();
		PlotEnergyDistribution();
	}
}

TH3D* EnergyDistributionModel::MultiplyElectronAndIonDensities(TH3D* electronDensities, TH3D* ionDensities)
{
	TH3D* result = (TH3D*)electronDensities->Clone("e-ion density");
	result->SetTitle("electron-ion density");
	result->Reset();

	double nXBins = electronDensities->GetXaxis()->GetNbins();
	double nYBins = electronDensities->GetYaxis()->GetNbins();
	double nZBins = electronDensities->GetZaxis()->GetNbins();

	for (int i = 1; i <= nXBins; i++) 
	{
		for (int j = 1; j <= nYBins; j++)
		{
			for (int k = 1; k <= nZBins; k++)
			{
				double value = electronDensities->GetBinContent(i, j, k) * ionDensities->GetBinContent(i, j, k);
				result->SetBinContent(i, j, k, value);
			}
		}
	}
	return result;
}	

void EnergyDistributionModel::LoadLabEnergyFile(std::filesystem::path file)
{
	if (!file.empty())
	{
		delete m_distribution;
		m_distribution = FileLoader::getInstance().LoadMatrixFile(file);
		m_distribution->SetTitle("lab energies");
	
		PlotDistribution();
		PlotLabEnergyProjections();
	}
}

void EnergyDistributionModel::CreateEnergyDistributionHistogram()
{
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

	delete energyDistribution;
	energyDistribution = new TH1D("Energy Distribution", "Energy Distribution", numberBins, binEdges.data());

}

void EnergyDistributionModel::GenerateEnergyDistribution()
{
	CreateEnergyDistributionHistogram();

	// sample positions from electron density multiplied with ion density given from outside
	std::vector<Point3D> positionSamples = mcmcSampler.GetSamples();
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

	for (const Point3D& point : positionSamples)
	{
		double x = point.x;
		double y = point.y;
		double z = point.z;

		double x_modified = std::min(std::max(x, m_distribution->GetXaxis()->GetBinCenter(1)), m_distribution->GetXaxis()->GetBinCenter(10) - 1e-3);
		double y_modified = std::min(std::max(y, m_distribution->GetYaxis()->GetBinCenter(1)), m_distribution->GetYaxis()->GetBinCenter(10) - 1e-3);
		double z_modified = std::min(std::max(z, m_distribution->GetZaxis()->GetBinCenter(1)), m_distribution->GetZaxis()->GetBinCenter(10) - 1e-3);

		// calculate velocity magnitude from lab energy given as matrix (TH3D) from outside
		double labEnergy = m_distribution->Interpolate(x_modified, y_modified, z_modified);
		double electronVelocityMagnitude = TMath::Sqrt(2 * labEnergy * TMath::Qe() / PhysicalConstants::electronMass);

		// determine direction of velocity based on beam trajectory function
		TVector3 longitudinalDirection = electronBeam.GetDirection(point.z);
		//longitudinalDirection.Print();
		TVector3 transverseDirection = longitudinalDirection.Orthogonal();

		// add random values to velocity in transverse and longitudinal directions:
		// - calculate longitudinal kT, transverse kT is fixed
		double long_kT = electronBeam.GetLongitudinal_kT(labEnergy);
		double trans_kT = electronBeam.GetTransverse_kT();

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
			+ 1 / sqrt(2) * transverseAddition * transverseDirection;

		// calculate collision velocity vector and magnitude using a fixed ion beam velocity
		double ionVelocityMagnitude = TMath::Sqrt(2 * electronBeam.GetParamters().coolingEnergy * TMath::Qe() / PhysicalConstants::electronMass); // calc from cooling energy;
		TVector3 ionVelocity(0, 0, ionVelocityMagnitude);

		double collosionVelocity = (finalElectronVelocity - ionVelocity).Mag();

		// calculate collision energy [eV] and put it in a histogram
		double collisionEnergy = 0.5 * PhysicalConstants::electronMass * collosionVelocity * collosionVelocity / TMath::Qe();
		energyDistribution->Fill(collisionEnergy);

		//std::cout << "position: (" << x << " " << y << " " << z << ")\n";
		////std::cout << "modified: (" << x_modified << " " << y_modified << " " << z_modified << ")\n";
		//std::cout << "lab energy: " << labEnergy << " eV\n";
		//std::cout << "electron velocity: " << electronVelocityMagnitude << " m/s\n";
		//std::cout << "kT trans: " << trans_kT << " eV, kT long: " << long_kT << " eV\n";
		//std::cout << "sigma trans: " << transSigma << " m/s, long sigma: " << longSigma << " m/s\n";
		//std::cout << "collision energy: " << collisionEnergy << " eV\n";
		//std::cout << "\n";
	}
	
}

void EnergyDistributionModel::PlotEnergyDistribution()
{
	if (!energyDistribution) return;

	m_mainCanvas->cd(1);
	energyDistribution->Draw();
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

