#include "pch.h"

#include "EnergyDistributionManager.h"
#include "EnergyDistribution.h"

std::unordered_map<double, EnergyDistribution*> EnergyDistribution::s_allDistributions;

EnergyDistribution::EnergyDistribution()
	: TH1D()
{
}

EnergyDistribution::~EnergyDistribution()
{
	// needs to be removed from hash map
	for (auto it = s_allDistributions.begin(); it != s_allDistributions.end(); it++)
	{
		if (it->second == this) 
		{
			s_allDistributions.erase(it); 
			break;
		}
	}
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
	eBeamParameter.detuningEnergy = pow(sqrt(labEnergiesParameter.centerLabEnergy) - sqrt(eBeamParameter.coolingEnergy), 2);
	
	SetupLabellingThings();
	std::vector<double> binEdges = SetupBinning();

	std::string histDescription = folder.filename().string() + " " + std::to_string(index);
	SetName(("Energy Distribution " + histDescription).c_str());
	SetTitle(("Energy Distribution " + histDescription).c_str());
	SetBins(binEdges.size() - 1, binEdges.data());
}

void EnergyDistribution::SetupLabellingThings()
{
	if (!eBeamParameter.densityFile.get().empty() && !labEnergiesParameter.energyFile.get().empty())
	{
		folder = eBeamParameter.densityFile.get().parent_path().parent_path();
		subFolder = Form("E_cool %.3feV I_e %.2eA r_ion %.4fm", eBeamParameter.coolingEnergy.get(),
			eBeamParameter.electronCurrent.get(), ionBeamParameter.radius.get());

		std::cout << eBeamParameter.densityFile.get().filename() << std::endl;
		index = std::stoi(eBeamParameter.densityFile.get().filename().string().substr(0, 4));
	}

	if (eBeamParameter.hasGaussianShape) tags += "e-gaus, ";
	if (eBeamParameter.hasNoBending) tags += "no bend, ";
	if (eBeamParameter.hasFixedLongitudinalTemperature) tags += "fixed kT||, ";
	if (labEnergiesParameter.useUniformEnergies) tags += "uniform energy, ";
	if (labEnergiesParameter.useOnlySliceXY) tags += Form("energy sliced %.3f, ", labEnergiesParameter.sliceToFill.get());
	if (eDistParameter.cutOutZValues) tags += Form("z samples %.3f - %.3f, ", eDistParameter.cutOutRange.get().x, eDistParameter.cutOutRange.get().y);
	label = Form("%d: U drift = %.2fV, E_d = %.4f", index, labEnergiesParameter.driftTubeVoltage.get(),
		eBeamParameter.detuningEnergy.get());
}

std::vector<double> EnergyDistribution::SetupBinning()
{
	EnergyDistributionManager* eDistManager = (EnergyDistributionManager*)Module::Get("Energy Distribution Manager");

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
			double difference = std::max(nextBin - binEdges[i], eDistParameter.minBinSize.get());
			binEdges.push_back(binEdges[i] + difference);
			//std::cout << binEdges[binEdges.size() - 1] << " " << nextBin << " " << difference << " \n";
		}
		std::cout << "number bins " << binEdges.size() - 1 << "\n";
	}
	return binEdges;
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
	std::string string = Form("# folder: %s\n", folder.filename().string().c_str()) +
		eBeamParameter.toString() +
		labEnergiesParameter.toString() +
		eDistParameter.toString() +
		ionBeamParameter.toString() +
		mcmcParameter.toString();

	return string;
}

std::string EnergyDistribution::Filename()
{
	std::ostringstream indexSS;
	std::ostringstream eCoolSS;
	eCoolSS << std::fixed << std::setprecision(3) << eBeamParameter.coolingEnergy.get();
	indexSS << std::setw(4) << std::setfill('0') << index;

	std::string string = indexSS.str() + std::string(Form(" E_d %.4feV.asc", eBeamParameter.detuningEnergy.get()));

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
