#pragma once
#include "Module.h"

struct mbrcData
{
	double E_det;
	double rateCoef;
	double RateCoefError;

	mbrcData(double Ed, double rc, double error) : E_det(Ed), rateCoef(rc), RateCoefError(error) {}
};

class RateCoefficient : public Module
{
public:
	RateCoefficient();

private:
	void ShowUI() override;

	void FillVectors(std::vector<mbrcData>& data);

	void PlotRateCoefficient();

private:
	TGraph* rateCoefficientsGraph = new TGraph();

	std::vector<double> detuningEnergies;
	std::vector<double> rateCoefficients;
	std::vector<double> rateCoefficientsError;

	std::vector<double> rateCoefficientsFit;

	// plot parameters
	bool logX = true;
	bool logY = true;
};

