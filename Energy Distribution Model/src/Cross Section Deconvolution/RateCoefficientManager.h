#pragma once
#include "Module.h"

class RateCoefficientManager : public CrossSectionDeconvolutionModule
{
public:
	RateCoefficientManager();

private:
	void ShowUI() override;

	void AddRateCoefficientToList(RateCoefficient& rc);
	void RemoveRateCoefficient(int index);

	void PlotRateCoefficient();

private:

	// plot parameters
	bool logX = true;
	bool logY = true;
};

