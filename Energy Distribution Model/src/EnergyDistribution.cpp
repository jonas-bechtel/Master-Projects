#include "pch.h"

#include "EnergyDistribution.h"
#include "AnalyticalDistribution.h"
#include "Constants.h"

std::unordered_map<double, EnergyDistribution*> EnergyDistribution::s_allDistributions;

EnergyDistribution::EnergyDistribution()
	: TH1D()
{
	//std::cout << "calling Energy Distribution default Constructor" << std::endl;
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
	//std::cout << "calling Energy Distribution destructor" << std::endl;
}

EnergyDistribution::EnergyDistribution(EnergyDistribution&& other)
	: TH1D(std::move(other))
{
	other.Reset();

	collisionEnergies	= std::move(other.collisionEnergies);
	binCenters			= std::move(other.binCenters);
	binValues			= std::move(other.binValues);
	binValuesNormalised = std::move(other.binValuesNormalised);
	fitX				= std::move(other.fitX);
	fitY				= std::move(other.fitY);

	mcmcParameter = other.mcmcParameter;
	eBeamParameter = other.eBeamParameter;
	ionBeamParameter = other.ionBeamParameter;
	labEnergiesParameter = other.labEnergiesParameter;
	analyticalParameter = other.analyticalParameter;
	simplifyParams = other.simplifyParams;

	label = std::move(other.label);
	tags = std::move(other.tags);
	folder = std::move(other.folder);
	subFolder = std::move(other.subFolder);
	index = std::move(other.index);

	rateCoefficient = other.rateCoefficient;
	psi = std::move(other.psi);

	plotted = other.plotted;
	showNormalisedByWidth = other.showNormalisedByWidth;

	other.ResetDefaultValues();

	//std::cout << "calling Energy Distribution Move Constructor" << std::endl;
}

EnergyDistribution& EnergyDistribution::operator=(EnergyDistribution&& other)
{
	if (this == &other) return *this;
	
	TH1D::operator=(std::move(other));

	collisionEnergies = std::move(other.collisionEnergies);
	binCenters = std::move(other.binCenters);
	binValues = std::move(other.binValues);
	binValuesNormalised = std::move(other.binValuesNormalised);
	fitX = std::move(other.fitX);
	fitY = std::move(other.fitY);

	mcmcParameter = other.mcmcParameter;
	eBeamParameter = other.eBeamParameter;
	ionBeamParameter = other.ionBeamParameter;
	labEnergiesParameter = other.labEnergiesParameter;
	analyticalParameter = other.analyticalParameter;
	simplifyParams = other.simplifyParams;

	label = std::move(other.label);
	tags = std::move(other.tags);
	folder = std::move(other.folder);
	subFolder = std::move(other.subFolder);
	index = std::move(other.index);

	rateCoefficient = other.rateCoefficient;
	psi = std::move(other.psi);

	plotted = other.plotted;
	showNormalisedByWidth = other.showNormalisedByWidth;

	other.ResetDefaultValues();

	//std::cout << "calling Energy Distribution Move assignment operator" << std::endl;

	return *this;
}

void EnergyDistribution::ResetDefaultValues()
{
	//std::cout << "binvalues size: " << binValues.size() << std::endl;
	folder = "Test";
	index = 0;
	rateCoefficient = 0;
	plotted = false;
	showNormalisedByWidth = true;
}

void EnergyDistribution::SetupLabellingThings()
{
	if (!eBeamParameter.densityFile.get().empty() && !labEnergiesParameter.energyFile.get().empty())
	{
		folder = eBeamParameter.densityFile.get().parent_path().parent_path();
		subFolder = Form("E_cool %.3feV I_e %.2eA", eBeamParameter.coolingEnergy.get(),
			eBeamParameter.electronCurrent.get());

		index = std::stoi(eBeamParameter.densityFile.get().filename().string().substr(0, 4));
	}

	if (simplifyParams.gaussianElectronBeam) tags += "e-gaus, ";
	if (simplifyParams.cylindricalElectronBeam) tags += "e-cylinder, ";
	if (simplifyParams.noElectronBeamBend) tags += "no bend, ";
	if (simplifyParams.fixedLongitudinalTemperature) tags += "fixed kT||, ";
	if (simplifyParams.uniformLabEnergies) tags += "uniform energy, ";
	if (simplifyParams.sliceLabEnergies) tags += Form("energy sliced %.3f, ", simplifyParams.sliceToFill.get());
	if (simplifyParams.cutOutZValues) tags += Form("z samples %.3f - %.3f, ", simplifyParams.cutOutRange.get().x, simplifyParams.cutOutRange.get().y);
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
	//std::cout << "delta E1: " << sqrt(pow((eBeamParameter.transverse_kT * log(2)), 2) + 16 * log(2) * eBeamParameter.longitudinal_kT * firstPeak) << std::endl;
	// estimate half of the width
	double estimatedPeakWidth1 = 2 * sqrt(pow((eBeamParameter.transverse_kT * log(2)), 2) + 16 * log(2) * eBeamParameter.longitudinal_kT_estimate * firstPeak);
	double estimatedPeakWidth2 = 0.7; //5 * sqrt(pow((eBeamParameter.transverse_kT * log(2)), 2) + 16 * log(2) * eBeamParameter.longitudinal_kT * secondPeak);

	// in case of constant bins
	double step = binSettings.normalStepSize;
	double smallStep = binSettings.peakStepSize;
	numberBins = (max - min) / step + 4 * estimatedPeakWidth2 / smallStep;

	// in case of factor binning
	int numberNormalBins = log10(max / min) * binSettings.binsPerDecade;
	int numberPeakBins = binSettings.binsAtPeak * binSettings.increasePeakResolution;
	numberBins = numberNormalBins + 2 * numberPeakBins;
	double normalFactor = TMath::Power((max / min), (1.0 / numberNormalBins));
	double peakFactor1 = TMath::Power(((firstPeak + estimatedPeakWidth1) / std::max(firstPeak - estimatedPeakWidth1, min)), (1.0 / numberPeakBins));
	double peakFactor2 = TMath::Power(((secondPeak + estimatedPeakWidth2) / std::max(secondPeak - estimatedPeakWidth2, min)), (2.0 / numberPeakBins));

	//std::cout << "peak factor 1: " << peakFactor1 << std::endl;
	//std::cout << "peak factor 2: " << peakFactor2 << std::endl;
	//std::cout << "normalFactor: " << normalFactor << std::endl;
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
			double propsedNextEdge = 0;
			if (binSettings.constantBinSize)
				propsedNextEdge = lastEdge + step;

			if (binSettings.factorBinning)
				propsedNextEdge = lastEdge * normalFactor;

			// check if the next bin edge would be inside a peak or steps over a peak
			if (std::abs(propsedNextEdge - firstPeak) < estimatedPeakWidth1 || (lastEdge - firstPeak) / (propsedNextEdge - firstPeak) < 0)
			{
				// put one edge at the start of the peak
				binEdges.push_back(std::max(firstPeak - estimatedPeakWidth1, min));
				//std::cout << "found first peak\n";
				// and then add more bins until the end of the peak
				while (binEdges.back() < firstPeak + estimatedPeakWidth1)
				{
					if (binSettings.constantBinSize)
						binEdges.push_back(binEdges.back() + smallStep);

					if (binSettings.factorBinning)
						binEdges.push_back(binEdges.back() * peakFactor1);
				}
			}
			// same for the second peak
			else if (std::abs(propsedNextEdge - secondPeak) < estimatedPeakWidth2 || (lastEdge - secondPeak) / (propsedNextEdge - secondPeak) < 0)
			{
				binEdges.push_back(std::max(secondPeak - estimatedPeakWidth2, min));
				//std::cout << "found second peak\n";
				while (binEdges.back() < secondPeak + estimatedPeakWidth2)
				{
					if (binSettings.constantBinSize)
						binEdges.push_back(binEdges.back() + smallStep);

					if (binSettings.factorBinning)
						binEdges.push_back(binEdges.back() * peakFactor2);
				}
			}
			// otherwise just add a constant value / multiply by constant factor to the previous edge
			else
			{
				binEdges.push_back(propsedNextEdge);
			}
		}
		else
		{
			if(binSettings.constantBinSize) 
				binEdges.push_back(lastEdge + step);

			if(binSettings.factorBinning)	
				binEdges.push_back(lastEdge * normalFactor);
		}
	}
	
	binCenters.reserve(binEdges.size() - 1);
	binValues.reserve(binEdges.size() - 1);
	binValuesNormalised.reserve(binEdges.size() - 1);

	std::cout << "number bins " << binEdges.size() - 1 << "\n";
	SetBins(binEdges.size() - 1, binEdges.data());
}

void EnergyDistribution::FillVectorsFromHist()
{
	// write non normalised values to the vector
	for (int i = 1; i <= GetNbinsX(); i++)
	{
		binCenters.push_back(GetBinCenter(i));
		binValues.push_back(GetBinContent(i));
	}

	// normalise the distribution to the width and to one
	double numberEntries = GetEntries();

	for (int i = 1; i <= GetNbinsX(); i++)
	{
		SetBinContent(i, GetBinContent(i) / (GetBinWidth(i) * numberEntries));
	}

	// write normalised values to the vector
	for (int i = 1; i <= GetNbinsX(); i++)
	{
		binValuesNormalised.push_back(GetBinContent(i));
	}
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

void EnergyDistribution::FitAnalyticalToPeak()
{
	double detuningEnergy = eBeamParameter.detuningEnergy;
	double kt_trans = eBeamParameter.transverse_kT;
	double kT_long = eBeamParameter.longitudinal_kT_estimate;

	double peakWidthGuess = sqrt(pow((kt_trans * log(2)), 2) + 16 * log(2) * kT_long * detuningEnergy);
	double energyMin = eBeamParameter.detuningEnergy - 1 * peakWidthGuess;
	double energyMax = eBeamParameter.detuningEnergy + 1 * peakWidthGuess;
	const int nParameter = 4;

	//std::cout << "min: " << energyMin << " max: " << energyMax << std::endl;	
	// first parameter is a scaling factor
	double initialGuess[nParameter] = { 1, detuningEnergy, kt_trans, kT_long };

	TF1* fitFunction = new TF1("fit function", AnalyticalEnergyDistributionFit, energyMin, energyMax, nParameter);
	fitFunction->SetParameters(initialGuess);
	//fitFunction->FixParameter(1, detuningEnergy);
	fitFunction->FixParameter(2, kt_trans);

	TFitResultPtr result = Fit(fitFunction, "QRN0");
	double maxValue = fitFunction->GetMaximum();
	double energyOfMax = fitFunction->GetMaximumX();
	double xLeft = fitFunction->GetX(maxValue / 2, 0, energyOfMax);
	double xRight = fitFunction->GetX(maxValue / 2, energyOfMax, 1000);

	double* fitParameter = fitFunction->GetParameters();
	analyticalParameter.detuningEnergy = fitParameter[1];
	//analyticalParameter.transverseTemperature = fitParameter[2];
	analyticalParameter.longitudinalTemperature = fitParameter[3];
	analyticalParameter.FWHM = xRight - xLeft;
	analyticalParameter.effectiveLength = fitParameter[0] * CSR::overlapLength;
	analyticalParameter.scalingFactor = fitParameter[0];

	// Fill distribution Fit data with these parameters
	double energyStep = 0.002;// (energyMax - energyMin) / 200;
	int i = 0;
	while (true)
	{
		double energy = 0 + i * energyStep;
		double value = AnalyticalEnergyDistributionFit(&energy, fitParameter);
		if (value > 1e-1)
		{
			fitX.push_back(energy);
			fitY.push_back(value);
		}
		else if (fitX.size() > 0)
		{
			break;
		}
		// emergency stop
		if (energy > 100) break;

		i++;
	}

	std::cout << "fit array size: " << fitX.size() << std::endl;
	fitFunction->Delete();
}

std::string EnergyDistribution::String()
{
	std::string string = Form("# folder: %s\n", (folder.filename().string() + subFolder.filename().string()).c_str()) +
		eBeamParameter.toString() +
		labEnergiesParameter.toString() +
		ionBeamParameter.toString() +
		mcmcParameter.toString() +
		analyticalParameter.toString();

	if (folder.filename().string() == "Test")
	{
		string += simplifyParams.toString();
	}

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
