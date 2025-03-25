#include "pch.h"
#include "EnergyDistribution.h"

#include "ElectronBeam.h"
#include "IonBeam.h"
#include "LabEnergies.h"
#include "MCMC.h"

#include "AnalyticalDistribution.h"
#include "Constants.h"
#include "FileUtils.h"

RNG_engine EnergyDistribution::generator = RNG_engine();

std::normal_distribution<double> EnergyDistribution::longitudinalNormalDistribution = std::normal_distribution<double>();
std::normal_distribution<double> EnergyDistribution::transverseNormalDistribution = std::normal_distribution<double>();

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

	label = std::move(other.label);
	tags = std::move(other.tags);
	
	index = std::move(other.index);

	psi = std::move(other.psi);

	showPlot = other.showPlot;
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

	label = std::move(other.label);
	tags = std::move(other.tags);
	//folder = std::move(other.folder);
	//subFolder = std::move(other.subFolder);
	index = std::move(other.index);

	psi = std::move(other.psi);

	cfData = std::move(other.cfData);

	showPlot = other.showPlot;
	showNormalisedByWidth = other.showNormalisedByWidth;

	other.ResetDefaultValues();

	//std::cout << "calling Energy Distribution Move assignment operator" << std::endl;

	return *this;
}

void EnergyDistribution::CopyParameters()
{
	mcmcParameter = MCMC::GetParameters();
	eBeamParameter = ElectronBeam::GetParameters();
	ionBeamParameter = IonBeam::GetParameters();
	labEnergiesParameter = LabEnergy::GetParameters();
}

void EnergyDistribution::ResetDefaultValues()
{
	index = 0;
	showPlot = false;
	showNormalisedByWidth = true;
}

void EnergyDistribution::SetupLabellingThings()
{
	if (!eBeamParameter.densityFile.get().empty() && !labEnergiesParameter.energyFile.get().empty())
	{
		index = std::stoi(eBeamParameter.densityFile.get().filename().string().substr(0, 4));
	}

	tags += MCMC::GetTags();
	tags += ElectronBeam::GetTags();
	tags += LabEnergy::GetTags();
	tags += IonBeam::GetTags();
	label = Form("%d: U drift = %.2fV, E_d = %.4f", index, labEnergiesParameter.driftTubeVoltage.get(),
		eBeamParameter.detuningEnergy.get());
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

	TF1* fitFunction = new TF1("fit function", AnalyticalDistribution::FitFunction, settings.initialRange[0], settings.initialRange[1], nParameter);
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
			double value = AnalyticalDistribution::FitFunction(&energy, fitParameter);
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
			double value = AnalyticalDistribution::FitFunction(&energy, fitParameter);
			fitX.push_back(energy);
			fitY.push_back(value);
		}
	}

	std::cout << "fit array size: " << fitX.size() << std::endl;
	fitFunction->Delete();
}

void EnergyDistribution::Generate(std::filesystem::path descriptionFile, int index, const BinningSettings& binSettings, const PeakFitSettings& fitSettings)
{
	// get all necessary modules
	std::filesystem::path folder = descriptionFile.parent_path();

	// get 3 parameters: U drift tube, electron current, center E lab if index is in file
	std::array<float, 3> additionalParameter = FileUtils::GetParametersFromDescriptionFileAtIndex(descriptionFile, index);

	// if they are not found the index is not in the file
	if (!additionalParameter[0])
	{
		std::cout << "index " << index << " is not in the file " << descriptionFile.filename() << std::endl;
		return;
	}
	
	// set read electron current and center lab energy
	LabEnergy::SetDriftTubeVoltage(additionalParameter[0]);
	LabEnergy::SetCenterEnergy(additionalParameter[2]);
	ElectronBeam::SetElectronCurrent(additionalParameter[1]);
	ElectronBeam::CalculateEstimateLongkT();
	ElectronBeam::CalculateDetuningEnergy();
	ElectronBeam::CalculateDetuningVelocity();

	// full procedure to generate one energy distribution 
	// 1. setup necessary distributions
	std::filesystem::path densityfile = FileUtils::FindFileWithIndex(folder / "e-densities", index);
	if (densityfile.empty())
	{
		std::cout << "density file not found: " << densityfile.filename() << std::endl;
		return;
	}
	ElectronBeam::SetupDistribution(densityfile);

	std::filesystem::path energyfile = FileUtils::FindFileWithIndex(folder / "lab-energies", index);
	if (energyfile.empty())
	{
		std::cout << "energy file not found: " << energyfile.filename() << std::endl;
		return;
	}
	LabEnergy::SetupDistribution(energyfile);

	IonBeam::CreateFromReference(ElectronBeam::Get());

	TH3D* ionElectronBeam = (TH3D*)ElectronBeam::Get()->Clone("ion electron beam");
	ionElectronBeam->Multiply(IonBeam::Get());
	MCMC::SetTargetDist(ionElectronBeam);

	// 2. sample from this distribution
	MCMC::GenerateSamples();

	// final setup of current distribution
	CopyParameters();
	SetupLabellingThings();
	SetupBinning(binSettings);

	// 3. generate energy distribution
	std::vector<Point3D> positionSamples = MCMC::GetSamples();
	collisionEnergies.reserve(positionSamples.size());

	if (positionSamples.empty())
	{
		std::cout << "no sampled positions were given\n";
		return;
	}

	for (const Point3D& point : positionSamples)
	{
		double x = point.x;
		double y = point.y;
		double z = point.z;

		// calculate velocity magnitude from lab energy given as matrix (TH3D) from outside
		double labEnergy = LabEnergy::Get(x, y, z);
		double electronVelocityMagnitude = TMath::Sqrt(2 * labEnergy * TMath::Qe() / PhysicalConstants::electronMass);

		// determine direction of velocity based on beam trajectory function
		TVector3 longitudinalDirection = ElectronBeam::GetDirection(point.z);
		TVector3 transverseDirection = longitudinalDirection.Orthogonal();

		// add random values to velocity in transverse and longitudinal directions:
		// - calculate longitudinal kT, transverse kT is fixed
		double long_kT = ElectronBeam::GetLongitudinal_kT(labEnergy);
		double trans_kT = ElectronBeam::GetTransverse_kT();

		// - use kT to calculate sigmas of gaussians
		double longSigma = TMath::Sqrt(long_kT * TMath::Qe() / PhysicalConstants::electronMass);
		double transSigma = TMath::Sqrt(trans_kT * TMath::Qe() / PhysicalConstants::electronMass);

		// - sample from gaussians with these sigmas and add that to the electron velocity
		longitudinalNormalDistribution = std::normal_distribution<double>(0, longSigma);
		transverseNormalDistribution = std::normal_distribution<double>(0, transSigma);
		double longitudinalAddition = longitudinalNormalDistribution(generator);
		double transverseAdditionX = transverseNormalDistribution(generator);
		double transverseAdditionY = transverseNormalDistribution(generator);

		// we need a vector that is never in line with the longitudinalDirection
		TVector3 helpVector = TVector3(1, 0, 0);
		TVector3 transverseDirection1 = longitudinalDirection.Cross(helpVector);
		TVector3 transverseDirection2 = longitudinalDirection.Cross(transverseDirection1);

		TVector3 finalElectronVelocity = (electronVelocityMagnitude + longitudinalAddition) * longitudinalDirection
			+ transverseAdditionX * transverseDirection1
			+ transverseAdditionY * transverseDirection2;

		// calculate collision velocity vector and magnitude using a fixed ion beam velocity
		double ionVelocityMagnitude = TMath::Sqrt(2 * eBeamParameter.coolingEnergy * TMath::Qe() / PhysicalConstants::electronMass); // calc from cooling energy;
		TVector3 ionVelocityDirection = IonBeam::GetDirection();
		TVector3 ionVelocity = ionVelocityDirection * ionVelocityMagnitude;

		TVector3 collisionVelocity = ionVelocity - finalElectronVelocity;
		double collisionVelocityMagnitude = collisionVelocity.Mag();

		// calculate collision energy [eV] and put it in a histogram
		double collisionEnergy = 0.5 * PhysicalConstants::electronMass * pow(collisionVelocityMagnitude, 2) / TMath::Qe();
		Fill(collisionEnergy);
		collisionEnergies.push_back(collisionEnergy);
	}

	FillVectorsFromHist();
	RemoveEdgeZeros();
	CalculateFWHM();
	FitAnalyticalToPeak(fitSettings);
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

void EnergyDistribution::Plot(bool showMarkers, bool showFit) const
{
	if (!showPlot)
	{
		return;
	}

	if (showMarkers) ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);

	if (showNormalisedByWidth)
	{
		ImPlot::PlotLine(label.c_str(), binCenters.data(), binValuesNormalised.data(), binValuesNormalised.size());

		ImVec4 color = ImPlot::GetLastItemColor();

		if (showFit)
		{
			color.x *= 2;
			color.y *= 2;
			color.z *= 2;

			// Plot the second line with a lighter color and dashed
			ImPlot::SetNextLineStyle(ImVec4(color.x, color.y, color.z, color.w));
			ImPlot::PlotLine("##", fitX.data(), fitY.data(), fitX.size(), ImPlotLineFlags_Segments);
		}
	}
	else
	{
		ImPlot::PlotLine(label.c_str(), binCenters.data(), binValues.data(), binCenters.size());
	}
}

void EnergyDistribution::SetPlot(bool plot)
{
	showPlot = plot;
}

void EnergyDistribution::SetNormalised(bool normalised)
{
	showNormalisedByWidth = normalised;
}

double EnergyDistribution::GetDetuningEnergy()
{
	return eBeamParameter.detuningEnergy;
}

void EnergyDistribution::ShowListItem()
{
	std::string labelTags = label;
	if (!tags.empty())
	{
		labelTags += "\n";
		labelTags += tags;
	}

	// Render each item as selectable
	if (ImGui::Selectable(labelTags.c_str(), showPlot, ImGuiSelectableFlags_AllowItemOverlap))
	{
		showPlot = !showPlot;
	}

	if (ImGui::BeginDragDropSource(ImGuiDragDropFlags_None))
	{
		ImGui::SetDragDropPayload("Analytical_Pars", &outputParameter, sizeof(OutputParameters));
		ImGui::Text("dragging stuff");
		ImGui::EndDragDropSource();
	}

	if (ImGui::BeginItemTooltip())
	{
		ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
		ImGui::TextUnformatted(GetHeaderString().c_str());
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}

	ImGui::SameLine();
	ImGui::Checkbox("normalised", &showNormalisedByWidth);
}

void EnergyDistribution::SaveSamples(std::filesystem::path folder) const
{
	std::filesystem::path file = folder / (Filename() + ".samples");
	std::ofstream outfile(file);

	if (!outfile.is_open())
	{
		std::cerr << "Error opening file" << std::endl;
		return;
	}

	outfile << GetHeaderString();

	outfile << "# sampled collision energy values\n";
	for (double energy : collisionEnergies)
	{
		outfile << energy << "\n";
	}

	outfile.close();
}

void EnergyDistribution::SaveHist(std::filesystem::path folder) const
{
	std::filesystem::path file = folder / (Filename() + ".asc");
	std::ofstream outfile(file);

	if (!outfile.is_open()) {
		std::cerr << "Error opening file" << std::endl;
		return;
	}
	outfile << GetHeaderString();

	outfile << "# bin center [eV]\tbin value\tbin value normalised by bin width\n";
	for (int i = 0; i < binCenters.size(); i++)
	{
		outfile << binCenters[i] << "\t";
		outfile << binValues[i] << "\t";
		outfile << binValuesNormalised[i] << "\n";
	}

	outfile.close();
}

void EnergyDistribution::Load(std::filesystem::path& file, bool loadSamples)
{
	// load the .asc file with the histogram data
	std::ifstream histFile(file);

	// Check if the file was successfully opened
	if (!histFile.is_open())
	{
		std::cerr << "Error: Could not open the file " << file << std::endl;
		return;
	}

	std::string header = FileUtils::GetHeaderFromFile(histFile);
	
	eBeamParameter.fromString(header);
	ionBeamParameter.fromString(header);
	labEnergiesParameter.fromString(header);
	mcmcParameter.fromString(header);
	outputParameter.fromString(header);

	SetupLabellingThings();

	std::string line;
	while (std::getline(histFile, line))
	{
		std::vector<std::string> tokens = FileUtils::SplitLine(line, "\t");

		binCenters.push_back(std::stod(tokens[0]));
		binValues.push_back(std::stod(tokens[1]));
		binValuesNormalised.push_back(std::stod(tokens[2]));
	}
	std::cout << "loaded file: " << file.filename();

	// see if .samples file exist with collision energy data
	if (loadSamples)
	{
		std::filesystem::path samplesFilename = file.replace_extension(".samples");
		if (std::filesystem::exists(samplesFilename))
		{
			std::ifstream samplesFile(samplesFilename);

			// Check if the file was successfully opened
			if (!samplesFile.is_open())
			{
				std::cerr << "Error: Could not open the file " << samplesFilename << std::endl;
				return;
			}
			// get header to get rid of it
			FileUtils::GetHeaderFromFile(samplesFile);

			collisionEnergies.reserve(mcmcParameter.numberSamples);

			while (std::getline(samplesFile, line))
			{
				collisionEnergies.push_back(std::stod(line));
			}
			std::cout << "\tsamples file found";
		}
	}
	std::cout << std::endl;
}

std::string EnergyDistribution::GetHeaderString() const
{
	std::string string = 
		eBeamParameter.toString() +
		labEnergiesParameter.toString() +
		ionBeamParameter.toString() +
		mcmcParameter.toString() +
		outputParameter.toString();

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

void PeakFitSettings::ShowWindow(bool& show)
{
	if (!show)
	{
		return;
	}
	if (ImGui::Begin("Peak Fit Settings", &show, ImGuiWindowFlags_NoDocking))
	{
		ImGui::BeginDisabled(adjustRange[0]);
		ImGui::SetNextItemWidth(100.0f);
		ImGui::InputDouble("##", &initialRange[0], 0.0, 0.0, "%.4e");
		ImGui::SameLine();
		ImGui::SetNextItemWidth(100.0f);
		ImGui::InputDouble("initial Range", &initialRange[1], 0.0, 0.0, "%.4e");
		ImGui::EndDisabled();

		//ImGui::SetNextItemWidth(160.0f);
		ImGui::SliderInt("number fit rounds", &fitRounds, 1, 4);
		for (int i = 0; i < fitRounds; i++)
		{
			ImGui::PushID(i);
			ImGui::BeginGroup();
			ImGui::Text(("fit " + std::to_string(i + 1)).c_str());
			ImGui::Checkbox("free E_d", &freeDetuningEnergy[i]);
			ImGui::Checkbox("free kT_long", &freekT_long[i]);
			ImGui::Checkbox("free kT_trans", &freeKT_trans[i]);
			ImGui::Checkbox("adjust range", &adjustRange[i]);
			ImGui::EndGroup();
			ImGui::SameLine();
		}
	}
	ImGui::End();
}

void BinningSettings::ShowWindow(bool& show)
{
	if (!show)
	{
		return;
	}
	if (ImGui::Begin("Binning settings", &show, ImGuiWindowFlags_NoDocking))
	{
		ImGui::SetNextItemWidth(150.0f);
		ImGui::InputFloat2("energy range", energyRange, "%.1e");
		ImGui::SameLine();
		ImGui::Checkbox("more bins at peaks", &increasePeakResolution);

		if (ImGui::Checkbox("factor binning", &factorBinning))
		{
			constantBinSize = !factorBinning;
		}
		ImGui::SameLine();
		ImGui::BeginDisabled(!factorBinning);
		ImGui::SetNextItemWidth(50.0f);
		ImGui::InputInt("bins per decade", &binsPerDecade, 0);
		ImGui::SameLine();
		ImGui::BeginDisabled(!increasePeakResolution);
		ImGui::SetNextItemWidth(50.0f);
		ImGui::InputInt("bins at peak", &binsAtPeak, 0);
		ImGui::EndDisabled();
		ImGui::EndDisabled();

		if (ImGui::Checkbox("constant bin size", &constantBinSize))
		{
			factorBinning = !constantBinSize;
		}
		ImGui::BeginDisabled(!constantBinSize);
		ImGui::SameLine();
		ImGui::SetNextItemWidth(80.0f);
		ImGui::InputDouble("step size", &normalStepSize, 0, 0, "%.6f");
		ImGui::SameLine();
		ImGui::BeginDisabled(!increasePeakResolution);
		ImGui::SetNextItemWidth(80.0f);
		ImGui::InputDouble("peak step size", &peakStepSize, 0, 0, "%.6f");
		ImGui::EndDisabled();
		ImGui::EndDisabled();

	}
	ImGui::End();
}
