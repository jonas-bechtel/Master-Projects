#include "EnergyDistribution.h"

#include <sstream>

std::unordered_map<double, EnergyDistribution*> EnergyDistribution::s_allDistributions;

EnergyDistribution::EnergyDistribution()
	: TH1D()
{
}

std::vector<double>& EnergyDistribution::GetCollisionEnergies()
{
	return collisionEnergies;
}

std::vector<double>& EnergyDistribution::GetBinCenters()
{
	return binCenters;
}

std::vector<double>& EnergyDistribution::GetBinValues()
{
	return binValues;
}

std::vector<double>& EnergyDistribution::GetNormalisedBinValues()
{
	return binValuesNormalised;
}

MCMC_Parameters& EnergyDistribution::GetMCMC_Parameter()
{
	return mcmcParameter;
}

ElectronBeamParameters& EnergyDistribution::GetElectronBeamParameter()
{
	return eBeamParameter;
}

IonBeamParameters& EnergyDistribution::GetIonBeamParameter()
{
	return ionBeamParameter;
}

LabEnergiesParameters& EnergyDistribution::GetLabEnergyParameter()
{
	return labEnergiesParameter;
}

EnergyDistributionParameters& EnergyDistribution::GetEnergyDistributionParameter()
{
	return eDistParameter;
}

void EnergyDistribution::SetLabel(std::string label)
{
	m_label = label;
}

double& EnergyDistribution::GetRateCoefficient()
{
	return rateCoefficient;
}

std::vector<double>& EnergyDistribution::GetPsis()
{
	return psi;
}

std::string& EnergyDistribution::GetLabel()
{
	return m_label;
}

std::string& EnergyDistribution::GetTags()
{
	return tags;
}

std::filesystem::path& EnergyDistribution::GetFolder()
{
	return folder;
}

std::filesystem::path& EnergyDistribution::GetSubFolder()
{
	return subFolder;
}

bool& EnergyDistribution::IsPlotted()
{
	return plotted;
}

bool& EnergyDistribution::IsPlottedNormalised()
{
	return showNormalisedByWidth;
}

void EnergyDistribution::SetupFromCurrentEnvironment()
{
	MCMC* mcmc = (MCMC*)Module::Get("MCMC");
	ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");
	IonBeam* ionBeam = (IonBeam*)Module::Get("Ion Beam");
	LabEnergies* labEnergies = (LabEnergies*)Module::Get("Lab Energies");
	EnergyDistributionManager* eDistManager = (EnergyDistributionManager*)Module::Get("Energy Distribution Manager");

	mcmcParameter = mcmc->GetParameter();
	eBeamParameter = eBeam->GetParameter();
	ionBeamParameter = ionBeam->GetParameter();
	labEnergiesParameter = labEnergies->GetParameter();
	eDistParameter = eDistManager->GetParameter();
	eDistParameter.detuningEnergy = pow(sqrt(labEnergiesParameter.centerLabEnergy) - sqrt(eBeamParameter.coolingEnergy), 2);

	if (!eBeamParameter.densityFile.empty() && !labEnergiesParameter.energyFile.empty())
	{
		folder = eBeamParameter.densityFile.parent_path().parent_path();
		subFolder = Form("E_cool %.3feV I_e %.2eA r_ion %.4fm", eBeamParameter.coolingEnergy,
			eBeamParameter.electronCurrent, ionBeamParameter.radius);
		index = std::stoi(eBeamParameter.densityFile.filename().string().substr(0, 4));
	}

	if (eBeamParameter.hasGaussianShape) tags += "e-gaus, ";
	if (eBeamParameter.hasNoBending) tags += "no bend, ";
	if (eBeamParameter.hasFixedLongitudinalTemperature) tags += "fixed kT||, ";
	if (labEnergiesParameter.useUniformEnergies) tags += "uniform energy, ";
	if (labEnergiesParameter.useOnlySliceXY) tags += Form("energy sliced %.3f, ", labEnergiesParameter.sliceToFill);
	if (eDistParameter.cutOutZValues) tags += Form("z samples %.3f - %.3f, ", eDistParameter.cutOutRange[0], eDistParameter.cutOutRange[1]);
	m_label = Form("%d: U drift = %.2fV, E_d = %.4f", index, eDistParameter.driftTubeVoltage,
		eDistParameter.detuningEnergy);

	//setup binning
	int numberBins;
	std::vector<double> binEdges;

	//if (parameter.energyDefinedBinning)
	//{
	//	double kT_trans = eBeamParameter.transverse_kT;
	//	double kT_long = eBeam->GetLongitudinal_kT(labEnergiesParameter.centerLabEnergy);
	//	double maxEnergy = 100; //just a gues, not fixed
	//	double factor = 1;
	//
	//	binEdges.push_back(0);
	//	binEdges.push_back(factor * kT_trans / 20);
	//
	//	for (int i = 1; binEdges[i] < kT_trans; i++)
	//	{
	//		binEdges.push_back(factor * 2 * binEdges[i]);
	//		//std::cout << binEdges[i + 1] << "\n";
	//	}
	//	while (binEdges[binEdges.size() - 1] < maxEnergy)
	//	{
	//		double previousEdge = binEdges[binEdges.size() - 1];
	//		double delta_E = sqrt(pow((kT_trans * log(2)), 2) + 16 * log(2) * kT_long * previousEdge);
	//		binEdges.push_back(previousEdge + factor * delta_E);
	//
	//		//std::cout << previousEdge << " " << delta_E << "\n";
	//		//std::cout << binEdges[binEdges.size()] << "\n";
	//		//std::cout << (binEdges[binEdges.size()] < maxEnergy) << "\n";
	//	}
	//
	//	numberBins = binEdges.size() - 1;
	//	std::cout << numberBins << "\n";
	//}
	//else
	if (!eDistParameter.limitBinSize)
	{
		float min = std::max(eDistManager->GetEnergyRange()[0], 1e-9f);
		float max = eDistManager->GetEnergyRange()[1];
		numberBins = log10(max / min) * eDistManager->GetBinsPerDecade();
		double factor = TMath::Power((max / min), (1.0 / numberBins));

		binCenters.reserve(numberBins);
		binValues.reserve(numberBins);

		binEdges.reserve(numberBins + 1);
		binEdges.push_back(min);
		for (int i = 0; i < numberBins; i++)
		{
			binEdges.push_back(binEdges[i] * factor);
		}
		std::cout << "number bins " << binEdges.size() - 1 << "\n";
	}
	else
	{
		double min = 1e-6; // std::max(energyRange[0], 1e-9f);
		float max = eDistManager->GetEnergyRange()[1];
		numberBins = log10(max / min) * eDistManager->GetBinsPerDecade();
		double factor = TMath::Power((max / min), (1.0 / numberBins));

		binCenters.reserve(numberBins);
		binValues.reserve(numberBins);

		binEdges.reserve(numberBins + 1);
		binEdges.push_back(0);
		for (int i = 0; binEdges[binEdges.size() - 1] < max; i++)
		{
			double nextBin = binEdges[i] * factor;
			double difference = std::max(nextBin - binEdges[i], eDistParameter.minBinSize);
			binEdges.push_back(binEdges[i] + difference);
			//std::cout << binEdges[binEdges.size() - 1] << " " << nextBin << " " << difference << " \n";
		}
		std::cout << "number bins " << binEdges.size() - 1 << "\n";
	}

	std::string histDescription = folder.filename().string() + " " + std::to_string(index);
	SetName(("Energy Distribution " + histDescription).c_str());
	SetTitle(("Energy Distribution " + histDescription).c_str());
	SetBins(binEdges.size() - 1, binEdges.data());
}

void EnergyDistribution::RemoveEdgeZeros()
{
	std::vector<double>& values = binValues;
	std::vector<double>& valuesNorm = binValuesNormalised;
	std::vector<double>& centers = binCenters;

	// Find the first non-zero element
	auto start = std::find_if(values.begin(), values.end(), [](double x) { return x != 0; });

	// Find the last non-zero element
	auto end = std::find_if(values.rbegin(), values.rend(), [](double x) { return x != 0; }).base();

	// Check if we have any non-zero elements at all
	if (start < end)
	{
		// Create new vectors with the range [start, end)
		int startIndex = std::distance(values.begin(), start);
		int endIndex = std::distance(values.begin(), end);
		values = std::vector<double>(start, end);

		valuesNorm = std::vector<double>(valuesNorm.begin() + startIndex,
			valuesNorm.begin() + endIndex);

		centers = std::vector<double>(centers.begin() + startIndex,
			centers.begin() + endIndex);
	}
	else
	{
		// Clear the vector if it's all zeros
		values.clear();
		valuesNorm.clear();
		centers.clear();
	}
}

std::string EnergyDistribution::String()
{
	std::string string = Form("# folder: %s\n", folder.filename().string()) +
		eBeamParameter.String() +
		labEnergiesParameter.String() +
		eDistParameter.String() +
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

	std::string string = indexSS.str() + std::string(Form(" E_d %.4feV.asc", eDistParameter.detuningEnergy));

	return string;
}

EnergyDistribution* EnergyDistribution::FindByEd(double detuningEnergy)
{
	if (s_allDistributions.find(detuningEnergy) == s_allDistributions.end())
	{
		std::cout << "no energy distribution with E_d = " << detuningEnergy << " was found\n";
		return nullptr;
	}
	return s_allDistributions.at(detuningEnergy);
}
