#include "pch.h"

#include "EnergyDistributionManager.h"
#include "EnergyDistribution.h"
#include "Constants.h"

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

void EnergyDistribution::SetupLabellingThings()
{
	if (!eBeamParameter.densityFile.get().empty() && !labEnergiesParameter.energyFile.get().empty())
	{
		folder = eBeamParameter.densityFile.get().parent_path().parent_path();
		subFolder = Form("E_cool %.3feV I_e %.2eA r_ion %.4fm", eBeamParameter.coolingEnergy.get(),
			eBeamParameter.electronCurrent.get(), ionBeamParameter.radius.get());

		index = std::stoi(eBeamParameter.densityFile.get().filename().string().substr(0, 4));
	}

	if (eBeamParameter.hasGaussianShape) tags += "e-gaus, ";
	if (eBeamParameter.hasCylindricalShape) tags += "e-cylinder, ";
	if (eBeamParameter.hasNoBending) tags += "no bend, ";
	if (eBeamParameter.hasFixedLongitudinalTemperature) tags += "fixed kT||, ";
	if (labEnergiesParameter.useUniformEnergies) tags += "uniform energy, ";
	if (labEnergiesParameter.useOnlySliceXY) tags += Form("energy sliced %.3f, ", labEnergiesParameter.sliceToFill.get());
	if (eDistParameter.cutOutZValues) tags += Form("z samples %.3f - %.3f, ", eDistParameter.cutOutRange.get().x, eDistParameter.cutOutRange.get().y);
	label = Form("%d: U drift = %.2fV, E_d = %.4f", index, labEnergiesParameter.driftTubeVoltage.get(),
		eBeamParameter.detuningEnergy.get());

	std::string histDescription = folder.filename().string() + " " + std::to_string(index);
	SetName(("Energy Distribution " + histDescription).c_str());
	SetTitle(("Energy Distribution " + histDescription).c_str());
}

void EnergyDistribution::SetupBinning(const BinningSettings& binSettings)
{
	int numberBins;
	std::vector<double> binEdges;

	double min = std::max(binSettings.energyRange[0], 1e-9f);
	double max = binSettings.energyRange[1];

	double firstPeak = eBeamParameter.detuningEnergy;
	double secondPeak = pow(sqrt(CSR::energyOutsideDriftTube) - sqrt(eBeamParameter.coolingEnergy), 2);
	// estimate half of the width
	double estimatedPeakWidth1 = 3 * 0.16; //3 * sqrt(pow((eBeamParameter.transverse_kT * log(2)), 2) + 16 * log(2) * eBeamParameter.longitudinal_kT * firstPeak);
	double estimatedPeakWidth2 = 1; //5 * sqrt(pow((eBeamParameter.transverse_kT * log(2)), 2) + 16 * log(2) * eBeamParameter.longitudinal_kT * secondPeak);

	if (binSettings.constantBinSize)
	{
		double step = binSettings.normalStepSize;
		double smallStep = binSettings.peakStepSize;
		numberBins = (max - min) / step + 4 * estimatedPeakWidth2 / smallStep;
		//std::cout << "guess of bin number: " << numberBins << std::endl;
		//std::cout << "estimate peak positions: " << firstPeak << ", " << secondPeak << std::endl;
		//std::cout << "estimate peak widths: " << estimatedPeakWidth1 << ", " << estimatedPeakWidth2 << std::endl;
		binEdges.reserve(numberBins);
		binEdges.push_back(min);

		while(binEdges.back() < max)
		{
			double lastEdge = binEdges.back();
			if (binSettings.increasePeakResolution)
			{
				// check if the next bin edge would be inside a peak or steps over a peak
				if (std::abs(lastEdge + step - firstPeak) < estimatedPeakWidth1 || (lastEdge - firstPeak) / (lastEdge + step - firstPeak) < 0)
				{
					// put one edge at the start of the peak
					binEdges.push_back(std::max(firstPeak - estimatedPeakWidth1, min));

					// and then add more bins until the end of the peak
					while (binEdges.back() < firstPeak + estimatedPeakWidth1)
					{
						binEdges.push_back(binEdges.back() + smallStep);
					}
				}
				// same for the second peak
				else if (std::abs(lastEdge + step - secondPeak) < estimatedPeakWidth2 || (lastEdge - secondPeak) / (lastEdge + step - secondPeak) < 0)
				{
					binEdges.push_back(std::max(secondPeak - estimatedPeakWidth2, min));

					while (binEdges.back() < secondPeak + estimatedPeakWidth2)
					{
						binEdges.push_back(binEdges.back() + smallStep);
					}
				}
				// otherwise just add a constant value to the previous edge
				else
				{
					binEdges.push_back(lastEdge + step);
				}
			}
			else
			{
				binEdges.push_back(lastEdge + step);
			}
		}
	}
	else if (binSettings.factorBinning)
	{
		numberBins = log10(max / min) * binSettings.binsPerDecade;
		double factor = TMath::Power((max / min), (1.0 / numberBins));
		
		binEdges.reserve(numberBins + 1);
		binEdges.push_back(min);

		for (int i = 0; i < numberBins; i++)
		{
			binEdges.push_back(binEdges.back() * factor);
		}
	}
	binCenters.reserve(binEdges.size() - 1);
	binValues.reserve(binEdges.size() - 1);
	binValuesNormalised.reserve(binEdges.size() - 1);

	std::cout << "number bins " << binEdges.size() - 1 << "\n";
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
	if (isAnalytical)
	{
		return analyticalParameter.toString();
	}

	bool excludeOptionals = !(folder.filename().string() == "Test");
	std::string string = Form("# folder: %s\n", folder.filename().string().c_str()) +
		eBeamParameter.toString(excludeOptionals) +
		labEnergiesParameter.toString(excludeOptionals) +
		eDistParameter.toString(excludeOptionals) +
		ionBeamParameter.toString(excludeOptionals) +
		mcmcParameter.toString(excludeOptionals) +
		analyticalParameter.toString(excludeOptionals);

	return string;
}

std::string EnergyDistribution::Filename()
{
	std::ostringstream indexSS;
	std::ostringstream eCoolSS;
	eCoolSS << std::fixed << std::setprecision(3) << eBeamParameter.coolingEnergy.get();
	indexSS << std::setw(4) << std::setfill('0') << index;

	std::string string = indexSS.str() + std::string(Form(" E_d %.4feV", eBeamParameter.detuningEnergy.get()));

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
