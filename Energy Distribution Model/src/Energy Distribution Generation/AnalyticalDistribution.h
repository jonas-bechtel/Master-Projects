#pragma once

namespace AnalyticalDistribution
{
	// functions for the analytical model
	double FitFunction(double* x, double* params);
	double Function(double* x, double* params);
	double Function(double Ecm, double Ed, double Ttr, double Tlong);
	double ComplexErrorFunction(double* x, double* par);
	double DawsonIntegral(double* x, double* par);
	double ExpDiff(double* x, double* par);

	void ShowWindow(bool& show);
	void Update();
	void Plot();
}
