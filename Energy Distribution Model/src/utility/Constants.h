#pragma once

namespace PhysicalConstants {

	// electron mass [kg]
	constexpr double electronMass = 9.10938291e-31;
	// electron mass [amu]
	constexpr double electronMassAMU = 5.4857990946e-4;
	// permitivity of vacuum [Coulomb / (V * m)]
	constexpr double epsilon_0 = 8.854187817e-12;
	// atomic mass unit [kg]
	constexpr double atomicMassUnit = 1.66053886e-27;
		
}

namespace CSR
{
	// taken from Daniel Paul PHD thesis Appendix C equation C.29 [m]
	constexpr double overlapLength = 1.13614;

	// drift Tube length where cooling force is doing something [m]
	constexpr double coolerLength = 0.8;
}