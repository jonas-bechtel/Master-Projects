#pragma once

namespace PhysicalConstants {

	// electron mass [kg]
	const double electronMass = 9.10938291e-31;
	// electron mass [amu]
	const double electronMassAMU = 5.4857990946e-4;
	// permitivity of vacuum [Coulomb / (V * m)]
	const double epsilon_0 = 8.854187817e-12;
	// atomic mass unit [kg]
	const double atomicMassUnit = 1.66053886e-27;
		
}

namespace CSR
{
	// taken from Daniel Paul PHD thesis Appendix C equation C.29 [m]
	const double overlapLength = 1.13614;

	// is a guess from looking at the lab energy files, used for estimating the second peak [eV]
	const double energyOutsideDriftTube = 44.0;
}