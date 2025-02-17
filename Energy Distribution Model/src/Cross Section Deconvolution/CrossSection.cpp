#include "pch.h"

#include "CrossSection.h"
#include "Constants.h"
#include "RateCoefficient.h"
#include "EnergyDistribution.h"

CrossSection::CrossSection()
{

}

void CrossSection::SetValues(double* newValues, bool square)
{
	for (int i = 0; i < values.size(); i++)
	{
		if (square)
		{
			hist->SetBinContent(i + 1, newValues[i] * newValues[i]);
			values[i] = newValues[i] * newValues[i];
		}
		else
		{
			hist->SetBinContent(i + 1, newValues[i]);
			values[i] = newValues[i];
		}
		
	}
}

void CrossSection::SetupBinning(CrossSectionBinningSettings binSettings, const RateCoefficient& rc)
{
	// first edge needs to be 0
	std::vector<double> binEdges;
	double minEnergy = 0;    
	double maxEnergy = 100;  //just a guess, not fixed

	// binning like in the paper 
	if (binSettings.scheme == PaperBinning)
	{
		//EnergyDistribution& representativeEnergyDist = energyDistributionList.back();
		double kT_trans = 0.002; // representativeEnergyDist.eBeamParameter.transverse_kT;
		double kT_long = 0.00047; // representativeEnergyDist->eBeamParameter.longitudinal_kT; //eBeam->GetLongitudinal_kT(labEnergiesParameter.centerLabEnergy);

		binEdges.push_back(0);
		binEdges.push_back(kT_trans / 20);
		int binFactor = 1;

		for (int i = 1; binEdges[i] < kT_trans; i++)
		{
			binEdges.push_back(binFactor * 2 * binEdges.back());
		}
		while (binEdges.back() < maxEnergy)
		{
			double previousEdge = binEdges.back();
			double delta_E = sqrt(pow((kT_trans * log(2)), 2) + 16 * log(2) * kT_long * previousEdge);
			binEdges.push_back(previousEdge + binFactor * delta_E);

			//std::cout << "delta E " << delta_E << ", last edge:	" << previousEdge << ", ratio: " << previousEdge /delta_E << "\n";
			if (previousEdge / delta_E > binSettings.maxRatio) break;
		}
		double lastEdgeSoFar = binEdges.back();

		// add edges so rc points are in bin center
		for (size_t i = rc.detuningEnergies.size() - 1; i > 0; i--)
		{
			if (rc.detuningEnergies.at(i) < lastEdgeSoFar) continue;
			
			double Ed_i = rc.detuningEnergies.at(i);
			double Ed_after_i = rc.detuningEnergies.at(i - 1);
			//std::cout << "Ed[i] = " << Ed_i << " Ed[i + 1] = " << Ed_after_i << " i = " << i << std::endl;
			binEdges.push_back((Ed_i + Ed_after_i) / 2);
		}
		// add one last edge
		binEdges.push_back(binEdges.back() + 2 * (rc.detuningEnergies.front() - binEdges.back()));
	}
	// bin width increses by a constant factor
	if (binSettings.scheme == FactorBinning)
	{
		double min = minEnergy;
		double factor = TMath::Power((maxEnergy / min), (1.0 / binSettings.numberBins));

		binEdges.reserve(binSettings.numberBins + 1);
		binEdges.push_back(min);
		for (int i = 0; i < binSettings.numberBins; i++)
		{
			binEdges.push_back(binEdges[i] * factor);
		}
	}
	if (binSettings.scheme == PaperFactorMix)
	{
		//EnergyDistribution& representativeEnergyDist = energyDistributionList.back();
		double kT_trans = 0.002; // representativeEnergyDist.eBeamParameter.transverse_kT;
		//double kT_long = 0.00047; // representativeEnergyDist->eBeamParameter.longitudinal_kT; //eBeam->GetLongitudinal_kT(labEnergiesParameter.centerLabEnergy);

		binEdges.push_back(0);
		binEdges.push_back(kT_trans / 20);

		for (int i = 1; binEdges[i] < kT_trans; i++)
		{
			binEdges.push_back(2 * binEdges[i]);
			//std::cout << binEdges[i + 1] << "\n";
		}

		double factor = TMath::Power((maxEnergy / binEdges.back()), (1.0 / binSettings.numberBins));

		for (int i = 0; i < binSettings.numberBins; i++)
		{
			binEdges.push_back(binEdges.back() * factor);
		}
	}

	//for (double edge : binEdges)
	//{
	//	std::cout << edge << "\n";
	//}
	for (int i = 1; i < binEdges.size(); i++)
	{
		std::cout << "bin center: " << (binEdges[i] + binEdges[i - 1]) / 2 << "\tbin width:	" << (binEdges[i] - binEdges[i - 1]) << std::endl;
	}
	std::cout << "number cross section bins: " << binEdges.size() - 1 << "\n";

	hist = new TH1D("cross section fit", "cross section fit", binEdges.size() - 1, binEdges.data());
}

void CrossSection::SetupInitialGuess(const RateCoefficient& rc, bool squareRoot)
{
	values.clear();
	values.reserve(hist->GetNbinsX());
	errors.clear();
	errors.resize(hist->GetNbinsX());

	for (int i = 1; i <= hist->GetNbinsX(); i++)
	{
		double energy = hist->GetBinCenter(i);
		double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
		double alpha = rc.graph->Eval(energy);
		if (squareRoot)
		{
			values.push_back(sqrt(alpha / velocity));
		}
		else
		{
			values.push_back(alpha / velocity);
		}
		
		initialGuess.push_back(alpha / velocity);
		energies.push_back(energy);
	}
}

void CrossSection::FillWithOneOverE(int scale)
{
	// setup binning
	std::vector<double> binEdges;
	double minEnergy = 1e-4;
	double maxEnergy = 100;
	int numberBins = 1000;

	double factor = TMath::Power((maxEnergy / minEnergy), (1.0 / numberBins));

	binEdges.reserve(numberBins + 2);

	binEdges.push_back(0);
	binEdges.push_back(minEnergy);
	for (int i = 0; i < numberBins; i++)
	{
		binEdges.push_back(binEdges.back() * factor);
	}
	hist = new TH1D("cross section 1/E", "cross section 1/E", binEdges.size() - 1, binEdges.data());

	// fill hist and lists
	values.reserve(hist->GetNbinsX());
	energies.reserve(hist->GetNbinsX());

	for (int i = 1; i <= hist->GetNbinsX(); i++)
	{
		double energy = hist->GetBinCenter(i);
		double value = scale / energy;
		
		energies.push_back(energy);
		values.push_back(value);
		hist->SetBinContent(i, value);
	}
}
