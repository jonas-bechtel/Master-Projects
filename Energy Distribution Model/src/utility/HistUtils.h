#pragma once
#include "Point3D.h"
namespace HistUtils
{
	double GetValueAtPosition(TH3D* hist, const Point3D& point, bool interpolate = true);
}

