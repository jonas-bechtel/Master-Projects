#include "pch.h"
#include "RateCoefficient.h"

RateCoefficient::RateCoefficient()
{
}

int RateCoefficient::GetIndexOfDetuningEnergy(double Ed)
{
    for (int i = 0; i < detuningEnergies.size(); i++)
    {
        if (detuningEnergies.at(i) == Ed)
        {
            //std::cout << detuningEnergies.at(i) << " == " << Ed << std::endl;
            return i;
        }
    }
    return -1;
}




