#include "pch.h"

#include "CrossSection.h"
#include "Constants.h"
#include "RateCoefficient.h"
#include "EnergyDistribution.h"

CrossSection::CrossSection()
{

}

void CrossSection::SetValues(double* newValues)
{
	for (int i = 0; i < values.size(); i++)
	{
		values[i] = newValues[i];
	}
}

void CrossSection::SetupBinning(CrossSectionBinningSettings binSettings)
{
	std::vector<double> binEdges;
	double maxEnergy = 100;  //just a gues, not fixed
	double minEnergy = 0;    //energyDistributionList.back().eBeamParameter.detuningEnergy / 10;

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
			//std::cout << binEdges.back() << "\n";
		}
		while (binEdges[binEdges.size() - 1] < maxEnergy)
		{
			double previousEdge = binEdges[binEdges.size() - 1];
			double delta_E = sqrt(pow((kT_trans * log(2)), 2) + 16 * log(2) * kT_long * previousEdge);
			binEdges.push_back(previousEdge + binFactor * delta_E);

			//std::cout << "delta E " << delta_E << "\n";
			//std::cout << binEdges[binEdges.size()] << "\n";
			//std::cout << (binEdges[binEdges.size()] < maxEnergy) << "\n";
		}
	}
	// bin width increses by a constant factor
	if (binSettings.scheme == FactorBinning)
	{
		float min = minEnergy;
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

	if (binSettings.scheme == Paper_FWHM)
	{

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

void CrossSection::SetupInitialGuess(const RateCoefficient& rc)
{
	values.clear();
	values.reserve(hist->GetNbinsX());

	for (int i = 1; i <= hist->GetNbinsX(); i++)
	{
		double energy = hist->GetBinCenter(i);
		double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
		double alpha = rc.graph->Eval(energy);
		values.push_back(alpha / velocity);
		initialGuess.push_back(alpha / velocity);
		energies.push_back(energy);
	}
}
