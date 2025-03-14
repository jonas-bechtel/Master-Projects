#include "pch.h"
#include "CoolingForceData.h"

void CoolingForceData::SetupHistograms(TH3D* reference)
{
	delete forceX;
	delete forceY;
	delete forceZ;
	delete positionSamples;

	forceX = (TH3D*)reference->Clone("cooling force X");
	forceY = (TH3D*)reference->Clone("cooling force Y");
	forceZ = (TH3D*)reference->Clone("cooling force Z");
	positionSamples = (TH3D*)reference->Clone("position samples");

	forceX->Reset();
	forceY->Reset();
	forceZ->Reset();
	positionSamples->Reset();

	forceX->SetTitle("cooling force X");
	forceY->SetTitle("cooling force Y");
	forceZ->SetTitle("cooling force Z");
	positionSamples->SetTitle("position samples");
}

void CoolingForceData::FillData()
{
	forceXIntegral = CalculateIntegral(forceX);
	forceYIntegral = CalculateIntegral(forceY);
	forceZIntegral = CalculateIntegral(forceZ);

	forceXValue = forceXIntegral / 0.8;
	forceYValue = forceYIntegral / 0.8;
	forceZValue = forceZIntegral / 0.8;

	xAxis.reserve(forceZ->GetNbinsX());
	yAxis.reserve(forceZ->GetNbinsY());
	zAxis.reserve(forceZ->GetNbinsZ());

	forceZProjectionX.reserve(forceZ->GetNbinsX());
	forceZProjectionY.reserve(forceZ->GetNbinsY());
	forceZProjectionZ.reserve(forceZ->GetNbinsZ());

	for (int i = 1; i <= forceZ->GetNbinsX(); i++)
	{
		xAxis.push_back(forceZ->GetXaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= forceZ->GetNbinsY(); i++)
	{
		yAxis.push_back(forceZ->GetYaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= forceZ->GetNbinsZ(); i++)
	{
		zAxis.push_back(forceZ->GetZaxis()->GetBinCenter(i));
	}

	TH1D* projectionX = forceZ->ProjectionX();
	TH1D* projectionY = forceZ->ProjectionY();
	TH1D* projectionZ = forceZ->ProjectionZ();

	for (int i = 1; i <= projectionX->GetNbinsX(); i++)
	{
		forceZProjectionX.push_back(projectionX->GetBinContent(i));
	}
	for (int i = 1; i <= projectionY->GetNbinsX(); i++)
	{
		forceZProjectionY.push_back(projectionY->GetBinContent(i));
	}
	for (int i = 1; i <= projectionZ->GetNbinsX(); i++)
	{
		forceZProjectionZ.push_back(projectionZ->GetBinContent(i));
	}

	delete projectionX;
	delete projectionY;
	delete projectionZ;
}

double CoolingForceData::CalculateIntegral(TH3D* hist)
{
	int nBinsZ = hist->GetZaxis()->GetNbins();
	double zMin = hist->GetZaxis()->GetXmin();
	double zMax = hist->GetZaxis()->GetXmax();

	// Create a 1D histogram for weighted averages along Z
	TH1D* H_weightedAvgZ = new TH1D("H_weightedAvgZ", "Weighted Average per Z", nBinsZ, zMin, zMax);

	// Loop over Z bins
	for (int iZ = 1; iZ <= nBinsZ; iZ++)
	{
		double weightedSum = 0.0;
		double weightSum = 0.0;

		// Loop over all (x, y) bins for this fixed Z bin
		for (int iX = 1; iX <= hist->GetXaxis()->GetNbins(); iX++)
		{
			for (int iY = 1; iY <= hist->GetYaxis()->GetNbins(); iY++)
			{
				int bin = hist->GetBin(iX, iY, iZ);

				double value = hist->GetBinContent(bin);
				double weight = positionSamples->GetBinContent(bin);

				weightedSum += value * weight;
				weightSum += weight;
			}
		}

		// Compute weighted average, handle division by zero
		double avg = (weightSum > 0) ? (weightedSum / weightSum) : 0.0;

		// Store result in the 1D histogram
		H_weightedAvgZ->SetBinContent(iZ, avg);
	}
	double result = H_weightedAvgZ->Integral("width");
	delete H_weightedAvgZ;

	return result;
}
