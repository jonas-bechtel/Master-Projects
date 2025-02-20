#include "pch.h"

#include "EnergyDistribution.h"
#include "AnalyticalDistribution.h"
#include "Constants.h"

double GetFWHM(TF1* function)
{
	double maxValue = function->GetMaximum();
	double XofMax = function->GetMaximumX();
	double xLeft = function->GetX(maxValue / 2, std::max(0.0, XofMax - 1), XofMax);
	double xRight = function->GetX(maxValue / 2, XofMax, XofMax + 1);

	return xRight - xLeft;
}

EnergyDistribution::EnergyDistribution()
	: TH1D()
{
	//std::cout << "calling Energy Distribution default Constructor" << std::endl;
}

EnergyDistribution::~EnergyDistribution()
{
	Delete();
	//std::cout << "calling Energy Distribution destructor" << std::endl;
}

EnergyDistribution::EnergyDistribution(EnergyDistribution&& other) noexcept
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
	outputParameter = other.outputParameter;
	simplifyParams = other.simplifyParams;

	label = std::move(other.label);
	tags = std::move(other.tags);
	//folder = std::move(other.folder);
	//subFolder = std::move(other.subFolder);
	index = std::move(other.index);

	psi = std::move(other.psi);

	plotted = other.plotted;
	showNormalisedByWidth = other.showNormalisedByWidth;

	other.ResetDefaultValues();

	//std::cout << "calling Energy Distribution Move Constructor" << std::endl;
}

EnergyDistribution& EnergyDistribution::operator=(EnergyDistribution&& other) noexcept
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
	outputParameter = other.outputParameter;
	simplifyParams = other.simplifyParams;

	label = std::move(other.label);
	tags = std::move(other.tags);
	//folder = std::move(other.folder);
	//subFolder = std::move(other.subFolder);
	index = std::move(other.index);

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
	//folder = "Test";
	index = 0;
	plotted = false;
	showNormalisedByWidth = true;
}

void EnergyDistribution::SetupLabellingThings()
{
	if (!eBeamParameter.densityFile.get().empty() && !labEnergiesParameter.energyFile.get().empty())
	{
		//folder = eBeamParameter.densityFile.get().parent_path().parent_path();
		//subFolder = Form("E_cool %.3feV I_e %.2eA", eBeamParameter.coolingEnergy.get(),
		//	eBeamParameter.electronCurrent.get());

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

	//std::string histDescription = folder.filename().string() + " " + std::to_string(index);
	//SetName(("Energy Distribution " + histDescription).c_str());
	//SetTitle(("Energy Distribution " + histDescription).c_str());
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
	peakFactor1 = std::min(peakFactor1, normalFactor);
	peakFactor2 = std::min(peakFactor2, normalFactor);

	//std::cout << "peak factor 1: " << peakFactor1 << std::endl;
	//std::cout << "peak factor 2: " << peakFactor2 << std::endl;
	//std::cout << "normalFactor: " << normalFactor << std::endl;
	//std::cout << "guess of bin number: " << numberBins << std::endl;
	//std::cout << "estimate peak positions: " << firstPeak << ", " << secondPeak << std::endl;
	//std::cout << "estimate peak widths: " << estimatedPeakWidth1 << ", " << estimatedPeakWidth2 << std::endl;
	binEdges.reserve(numberBins);
	binEdges.push_back(0);
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
				if(firstPeak - estimatedPeakWidth1 > min)
					binEdges.push_back(firstPeak - estimatedPeakWidth1);
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
		//std::cout << GetBinCenter(i) << " : " << GetBinContent(i) << " : " << GetBinWidth(i) << std::endl;
	}

	// normalise the distribution to the width and to one
	double numberEntries = GetEntries();

	for (int i = 1; i <= GetNbinsX(); i++)
	{
		SetBinContent(i, GetBinContent(i) / (GetBinWidth(i) * numberEntries));
		//std::cout << GetBinCenter(i) << " : " << GetBinContent(i) << std::endl;
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

void EnergyDistribution::CalculateFWHM()
{
	auto maxValueIt = std::max_element(binValuesNormalised.begin(), binValuesNormalised.end());
	int index = std::distance(binValuesNormalised.begin(), maxValueIt);
	double maxValue = *maxValueIt;
	double energyOfMaxValue = binCenters.at(index);
	//std::cout << "maxvalue: " << maxValue << " index: " << index << std::endl;
	if (std::abs(energyOfMaxValue - eBeamParameter.detuningEnergy) > 1)
	{
		std::cout << "maximum value position and detuning energy do not match: " << energyOfMaxValue << " != " << eBeamParameter.detuningEnergy << std::endl;
	}

	double energyRight;
	for (double energy = energyOfMaxValue; true; energy += 0.0001)
	{
		if (Interpolate(energy) < maxValue / 2)
		{
			energyRight = energy;
			break;
		}
	}
	double energyLeft;
	for (double energy = energyOfMaxValue; true; energy -= 0.0001)
	{
		if (energy < 0)
		{
			std::cout << "left side value did not go below half maximum" << std::endl;
			energyLeft = 0;
			break;
		}
		if (Interpolate(energy) < maxValue / 2)
		{
			energyLeft = energy;
			break;
		}
	}

	outputParameter.FWHM = energyRight - energyLeft;
}

void EnergyDistribution::FitAnalyticalToPeak(const PeakFitSettings& settings)
{
	double detuningEnergy = eBeamParameter.detuningEnergy;
	double kt_trans = eBeamParameter.transverse_kT;
	double kT_long = eBeamParameter.longitudinal_kT_estimate;

	constexpr int nParameter = 4;
	
	// first parameter is a scaling factor
	double initialGuess[nParameter] = { 1, detuningEnergy, kt_trans, kT_long };

	TF1* fitFunction = new TF1("fit function", AnalyticalEnergyDistributionFit, settings.initialRange[0], settings.initialRange[1], nParameter);
	fitFunction->SetParameters(initialGuess);
	fitFunction->FixParameter(1, detuningEnergy);
	fitFunction->FixParameter(2, kt_trans);
	fitFunction->FixParameter(3, kT_long);

	for (int i = 0; i < settings.fitRounds; i++)
	{
		if (settings.adjustRange[i])
		{
			//double* parameter = fitFunction->GetParameters();
			//double deltaE = sqrt(pow((parameter[2] * log(2)), 2) + 16 * log(2) * parameter[2] * parameter[1]);
			//double fwhm = GetFWHM(fitFunction);
			//std::cout << "fwhm gues: " << fwhm << std::endl;
			///double peakWidthGuess = std::max(0.001, fwhm / 2);
			double maxValue = fitFunction->GetMaximum();
			double XofMax = fitFunction->GetMaximumX();
			double energyMin = fitFunction->GetX(maxValue / 8, std::max(0.0, XofMax - 1), XofMax);   //std::max(0.0, fitFunction->GetMaximumX() - 1 * peakWidthGuess);
			double energyMax = fitFunction->GetX(maxValue / 8, XofMax, XofMax + 1);     //std::max(0.002, fitFunction->GetMaximumX() + 1 * peakWidthGuess);
			fitFunction->SetRange(energyMin, energyMax);

			//std::cout << "Fit round " << i << ": peakWidthGuess: " << peakWidthGuess << ", energyMin : " << energyMin << ", energyMax : " << energyMax << std::endl;
		}
		
		double* parameter = fitFunction->GetParameters();
		//std::cout << parameter[1] << ", " << parameter[2] << ", " << parameter[3] << std::endl;
		settings.freeDetuningEnergy[i] ? fitFunction->ReleaseParameter(1) : fitFunction->FixParameter(1, parameter[1]);
		settings.freeKT_trans[i] ? fitFunction->ReleaseParameter(2) : fitFunction->FixParameter(2, parameter[2]);
		settings.freekT_long[i] ? fitFunction->ReleaseParameter(3) : fitFunction->FixParameter(3, parameter[3]);
		Fit(fitFunction, "QRN0");
	}
	
	
	//std::cout << "maxValue: " << maxValue << " maxValue / 2: " << maxValue / 2 << " energyOfMax: " << energyOfMax << std::endl;
	//std::cout << "xLeft: " << xLeft << " xRight: " << xRight << std::endl;
	//std::cout << "f(xLeft) = " << fitFunction->Eval(xLeft) << " f(xRight) = " << fitFunction->Eval(xRight) << std::endl;

	// set all fit results
	double* fitParameter = fitFunction->GetParameters();
	double rangeMin, rangeMax;
	fitFunction->GetRange(rangeMin, rangeMax);
	outputParameter.fitRange = { (float)rangeMin, (float)rangeMax };
	outputParameter.fitDetuningEnergy = fitParameter[1];
	outputParameter.fitTransverseTemperature = fitParameter[2];
	outputParameter.fitLongitudinalTemperature = fitParameter[3];
	outputParameter.fitFWHM = GetFWHM(fitFunction);
	outputParameter.fitScalingFactor = fitParameter[0];
	outputParameter.effectiveLength = fitParameter[0] * CSR::overlapLength;

	// Fill distribution Fit data with these parameters
	if (settings.showLimitedFit)
	{
		double energyStep = (rangeMax - rangeMin) / 200;
		int i = 0;
		while (i < 1e5)
		{
			double energy = rangeMin + i * energyStep;
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
	}
	else
	{
		double energyStep = (rangeMax - rangeMin) / 20000;

		for (double energy = rangeMin; energy <= rangeMax; energy += energyStep)
		{
			double value = AnalyticalEnergyDistributionFit(&energy, fitParameter);
			fitX.push_back(energy);
			fitY.push_back(value);
		}
	}

	std::cout << "fit array size: " << fitX.size() << std::endl;
	fitFunction->Delete();
}

void EnergyDistribution::CalculatePsisFromBinning(TH1D* crossSection)
{
	psi.clear();
	int nBins = crossSection->GetNbinsX();
	//std::cout << nBins << std::endl;
	//distribution.psi.reserve(nBins);
	psi.resize(nBins);

	//std::cout << "distribution: " << distribution.index << std::endl;
	for (double energy : collisionEnergies)
	{
		int bin = crossSection->FindBin(energy);
		if (bin == 0)
		{
			std::cout << "crossection hist too small, got underflow when putting in " << energy << std::endl;
			continue;
		}
		//if (bin < 10) std::cout << "bin " << bin << ": " << energy << std::endl;
		double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
		if (bin - 1 >= psi.size())
		{
			//std::cout << "want to access " << bin - 1 << " but size is " << distribution.psi.size() << std::endl;
			//std::cout << energy << std::endl;
		}
		psi[bin - 1] += velocity;
	}
	for (int i = 0; i < psi.size(); i++)
	{
		psi[i] /= collisionEnergies.size();
		//std::cout << "Psi_" << i << ": " << distribution.psi[i] << "\t" << crossSectionFit->GetBinLowEdge(i+1)
		//	 << " - " << crossSectionFit->GetBinLowEdge(i+2) << "\n";
	}
}

double EnergyDistribution::CalculateTestRateCoefficient()
{
	if (collisionEnergies.empty()) return 0.0;
	double rc = 0.0;
	for (const double collisionEnergy : collisionEnergies)
	{
		double crossSectionValue = 1 / collisionEnergy;
		double collosionVelocity = TMath::Sqrt(2 * collisionEnergy * TMath::Qe() / PhysicalConstants::electronMass);
		rc += crossSectionValue * collosionVelocity;
	}
	rc /= collisionEnergies.size();

	return rc;
}

std::string EnergyDistribution::String() const
{
	std::string string = //Form("# folder: %s\n", (folder.filename().string() + subFolder.filename().string()).c_str()) +
		eBeamParameter.toString() +
		labEnergiesParameter.toString() +
		ionBeamParameter.toString() +
		mcmcParameter.toString() +
		outputParameter.toString();

	//if (folder.filename().string() == "Test")
	//{
	//	string += simplifyParams.toString();
	//}

	return string;
}

std::string EnergyDistribution::Filename() const
{
	std::ostringstream indexSS;
	std::ostringstream eCoolSS;
	eCoolSS << std::fixed << std::setprecision(3) << eBeamParameter.coolingEnergy.get();
	indexSS << std::setw(4) << std::setfill('0') << index;

	std::string string = indexSS.str() + std::string(Form(" E_d %.4feV", eBeamParameter.detuningEnergy.get()));

	return string;
}

void EnergyDistributionSet::AddDistribution(EnergyDistribution&& distribution)
{
	info.AddDistributionValues(distribution);

	// will call move Constructor
	distributions.emplace_back(std::move(distribution));
	EnergyDistribution& justMoved = distributions.back();
	//std::cout << "E_d: " << justMoved.eBeamParameter.detuningEnergy << std::endl;
	EdToDistMap[justMoved.eBeamParameter.detuningEnergy] = &justMoved;
}

EnergyDistribution* EnergyDistributionSet::FindByEd(double detuningEnergy)
{
	if (EdToDistMap.find(detuningEnergy) == EdToDistMap.end())
	{
		std::cout << "no energy distribution with E_d = " << detuningEnergy << " was found\n";
		return nullptr;
	}
	return EdToDistMap.at(detuningEnergy);
}

void EnergyDistributionSet::SetAllPlotted(bool plotted)
{
	for (EnergyDistribution& eDist : distributions)
	{
		eDist.plotted = plotted;
	}
}

void EnergyDistributionSet::SetAllShowNormalised(bool showNormalised)
{
	for (EnergyDistribution& eDist : distributions)
	{
		eDist.showNormalisedByWidth = showNormalised;
	}
}

void EnergyDistributionSet::CalculatePsisFromBinning(TH1D* crossSection)
{
	for (EnergyDistribution& eDist: distributions)
	{
		eDist.CalculatePsisFromBinning(crossSection);
	}
}

