#pragma once

// functions for the analytical model
double AnalyticalEnergyDistributionFit(double* x, double* params);
double AnalyticalEnergyDistribution(double* x, double* params);
double AnalyticalEnergyDistribution(double Ecm, double Ed, double Ttr, double Tlong);
double ComplexErrorFunction(double* x, double* par);
double DawsonIntegral(double* x, double* par);
double ExpDiff(double* x, double* par);
