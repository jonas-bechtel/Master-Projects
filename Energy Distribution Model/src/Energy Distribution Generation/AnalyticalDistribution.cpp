#include "pch.h"

#include "AnalyticalDistribution.h"
#include "ParameterImplementations.h"

namespace AnalyticalDistribution
{
	// Analytical Parameter
	static float scale = 1.0f;
	static float E_d = 10.0f;
	static float kT_long = 1e-4f;
	static float kT_trans = 0.002f;
	static float energyRange[2] = { 0.0f, 1.0f };

	// plot data
	static float energies[200];
	static float values[200];
	static float color[3] = { 1.0f, 0.0f, 0.0f };


	double FitFunction(double* x, double* params)
	{
		double factor = params[0];
		return factor * Function(x, &params[1]);
	}

	double Function(double* x, double* par)
	{
		// Flaterned energy distribution, expects transverse mean velocity = 0!
		// Original version written by Andreas Wolf (edistr.C) returns the same result normalised by 1/sqrt(ee)!!!
		//
		// x[0] energy in CM (the variable of the distribution) [eV]
		// par[0] detuning energy (average CM energy) [eV]
		// par[1] transversal electron temperature [eV]
		// par[2] longitudinal electron temperature [eV]

		Double_t ee = x[0];			//energy in CM (variable of the distribution) [eV]
		Double_t edet = par[0]; 	//detuning energy (average CM energy) [eV]
		Double_t tperp = par[1];	//transversal temperature [eV]
		Double_t tpar = par[2];		//longitudinal temperature [eV]
		//=============================
		Double_t rpar[1];
		rpar[0] = 0.;
		Double_t rx[1];
		//=============================
		const Double_t sqrtpi = 1.77245385090551588;  // TMath::Sqrt(TMath::Pi())
		const Double_t eps = 1e-5;
		//cout<<"edis: "<<ee<<" "<<edet<<" "<<tperp<<" "<<tpar<<" "<<endl;
		ee = TMath::Max(ee, 1e-20);
		tperp = TMath::Max(tperp, 1e-20);
		tpar = TMath::Max(tpar, 1e-20);
		//=============================
		Double_t eta = 1.0 - tpar / tperp;
		Double_t ksi = TMath::Abs(eta);

		if (ksi > eps)
		{
			Double_t aa = TMath::Sqrt(ksi * ee / tpar);
			Double_t bb = TMath::Sqrt(edet / ksi / tpar);
			//=============================
			Double_t exp1 = TMath::Sqrt(ee) - TMath::Sqrt(edet);
			exp1 = exp1 * exp1 / tpar;
			Double_t exp2 = TMath::Sqrt(ee) + TMath::Sqrt(edet);
			exp2 = exp2 * exp2 / tpar;
			//=============================
			if (eta > 0)
			{
				// Tpar < Tperp
				Double_t f = TMath::Exp(-exp1);
				if (bb >= aa)
				{
					rx[0] = bb - aa; f = f * ComplexErrorFunction(rx, rpar);
				}
				else
				{
					rx[0] = aa - bb; Double_t cc = (ee - edet / ksi) / tperp;
					f = 2. * TMath::Exp(-cc) - f * ComplexErrorFunction(rx, rpar);
				}
				rx[0] = bb + aa; f = f - TMath::Exp(-exp2) * ComplexErrorFunction(rx, rpar);
				f = f / (2. * tperp * TMath::Sqrt(tpar) * aa);  //factor for f(E)/sqrt(E)

				return f * sqrt(ee);
			}
			else
			{
				// Tpar > Tperp
				Double_t f = TMath::Exp(-exp1);
				rx[0] = bb + aa; f = f * DawsonIntegral(rx, rpar);
				rx[0] = bb - aa; f = f - TMath::Exp(-exp2) * DawsonIntegral(rx, rpar);
				f = f / (sqrtpi * tperp * TMath::Sqrt(tpar) * aa);  //factor for f(E)/sqrt(E)

				return f * sqrt(ee);
			}
		}
		else
		{
			// Tpar = Tperp
			Double_t aa = 2.0 * TMath::Sqrt(ee * edet) / tperp;
			Double_t f = 1 / (sqrtpi * tperp * TMath::Sqrt(tperp));  //factor for f(E)/sqrt(E)

			if (aa > 1)
			{
				Double_t exp1 = (TMath::Sqrt(ee) - TMath::Sqrt(edet));
				exp1 = exp1 * exp1 / tperp;
				Double_t exp2 = (TMath::Sqrt(ee) + TMath::Sqrt(edet));
				exp2 = exp2 * exp2 / tperp;
				f = f * (TMath::Exp(-exp1) - TMath::Exp(-exp2)) / aa;
			}
			else
			{
				rx[0] = aa;
				f = f * TMath::Exp(-(ee + edet) / tperp) * ExpDiff(rx, rpar);
			}
			return f * sqrt(ee);
		}
	}

	double Function(double Ecm, double Ed, double Ttr, double Tlong)
	{
		// The same as previous, different parameter representation
		//
		// Ecm energy in CM (the variable of the distribution) [eV]
		// Ed detuning energy (average CM energy) [eV]
		// Ttr transversal electron temperature [eV]
		// Tlong longitudinal electron temperature [eV]

		double x[1];
		x[0] = Ecm;

		double par[3];
		par[0] = Ed;
		par[1] = Ttr;
		par[2] = Tlong;

		return Function(x, par);
	}

	double ComplexErrorFunction(double* x, double* par)
	{
		// Special implementation of the complex errorfunction:
		// Returns real part of complex errorfunction (Voigt-Faddeeva function) of the argument representing complex part of the parameter
		//
		// x[0] - complex part of the parameter
		Double_t yy = x[0];
		//=============================
		Double_t factor1 = 2.82842712474619029; // 2.*TMath::Sqrt(2.)
		Double_t factor2 = 2.50662827463100024; // TMath::Sqrt(TMath::TwoPi())
		Double_t f = factor2 * TMath::Voigt(0., 1., factor1 * yy, 5);
		return f;
	}

	double DawsonIntegral(double* xarg, double* par)
	{
		// Returns Dawson's integral for any real x.
		//        (From Numerical Recipes ch. 6.10)
		Int_t i, n0;
#define NMAX 6        // #define NMAX 6 
		const Double_t H = 0.4;      // #define H 0.4 
		const Double_t A1 = 2.0 / 3.0; // #define A1 (2.0/3.0) 
		const Double_t A2 = 0.4;   // #define A2 0.4 
		const Double_t A3 = 2.0 / 7.0; // #define A3 (2.0/7.0) 
		Double_t d1, d2, e1, e2, sum, x2, xp, xx, ans;
		Double_t x = xarg[0];
		static Double_t c[NMAX + 1];
		static Int_t init = 0; //  Flag is 0 if we need to initialize, else 1.

		if (init == 0)
		{
			init = 1;
			for (i = 1; i <= NMAX; i++)
				c[i] = TMath::Exp(-H * H * (Double_t)((2 * i - 1) * (2 * i - 1)));
		}
		if (TMath::Abs(x) < 0.2)
		{   // Use series expansion. 
			x2 = x * x;
			ans = x * (1.0 - A1 * x2 * (1.0 - A2 * x2 * (1.0 - A3 * x2)));
		}
		else
		{
			//Use sampling theorem representation. 
			xx = TMath::Abs(x);
			//    n0=2*(Int_t)(0.5*xx/H+0.5); 
			n0 = 2 * TMath::Nint(0.5 * xx / H);
			xp = xx - (Double_t)n0 * H;
			e1 = TMath::Exp(2.0 * xp * H);
			e2 = e1 * e1;
			d1 = (Double_t)n0 + 1.0;
			d2 = d1 - 2.0;
			sum = 0.0;
			//    for (i=1;i<=NMAX;i++,d1+=2.0,d2-=2.0,e1*=e2) {
			for (i = 1; i <= NMAX; i++, d1 += 2.0, d2 -= 2.0, e1 *= e2)
			{
				sum += c[i] * (e1 / d1 + 1.0 / (d2 * e1));
			}

			ans = 0.5641895835 * TMath::Sign(TMath::Exp(-xp * xp), x) * sum;
			//Constant is 1/sqrt(pi) 
		}
		return ans;
	}

	double ExpDiff(double* x, double* par)
	{
		Double_t yy = x[0];
		const Double_t A = 1. / 3.;
		//=============================
		Double_t f;
		if (TMath::Abs(yy) > 1e-6)
		{
			f = (TMath::Exp(yy) - TMath::Exp(-yy)) / yy;
		}
		else
		{
			f = 2. + A * yy * yy;
		}
		return f;
	}

	void ShowWindow(bool& show)
	{
		if (!show)
		{
			return;
		}
		if (ImGui::Begin("Analytical parameters", &show, ImGuiWindowFlags_NoDocking))
		{
			ImGui::BeginGroup();
			bool changed = false;
			changed |= ImGui::SliderFloat2("range", energyRange, 1e-6f, 100.0f, "%.6f", ImGuiSliderFlags_Logarithmic);
			changed |= ImGui::SliderFloat("scale", &scale, 0.0f, 2.0f);
			changed |= ImGui::SliderFloat("E_d", &E_d, 0.0f, 100.0f, "%.8f", ImGuiSliderFlags_Logarithmic);
			changed |= ImGui::SliderFloat("kT long", &kT_long, 1e-6f, 0.1f, "%.6f", ImGuiSliderFlags_Logarithmic);
			changed |= ImGui::SliderFloat("kT trans", &kT_trans, 0.0f, 0.1f, "%.6f", ImGuiSliderFlags_Logarithmic);
			ImGui::ColorEdit3("color", color);
			ImGui::EndGroup();
			if (changed)
			{
				Update();
			}

			if (ImGui::BeginDragDropTarget())
			{
				if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("Analytical_Pars"))
				{
					OutputParameters& parameter = *(OutputParameters*)payload->Data;
					energyRange[0] = std::max(1e-6f, parameter.fitRange.get().x);
					energyRange[1] = parameter.fitRange.get().y;
					scale = parameter.fitScalingFactor;
					E_d = parameter.fitDetuningEnergy;
					kT_long = parameter.fitLongitudinalTemperature;
					kT_trans = parameter.fitTransverseTemperature;
					Update();
				}
				ImGui::EndDragDropTarget();
			}

		}
		ImGui::End();
	}

	void Update()
	{
		float step = (energyRange[1] - energyRange[0]) / 199;
		for (int i = 0; i < 200; i++)
		{
			energies[i] = energyRange[0] + i * step;
			values[i] = scale * Function(energies[i], E_d, kT_trans, kT_long);
		}
	}

	void Plot()
	{
		ImPlot::PushStyleColor(ImPlotCol_Line, { color[0], color[1], color[2], 1.0 });
		ImPlot::PlotLine("analytical", energies, values, 200, ImPlotLineFlags_Segments);
		ImPlot::PopStyleColor();
	}

}
