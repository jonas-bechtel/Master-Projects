#include "pch.h"
#include "HistUtils.h"

double HistUtils::GetValueAtPosition(TH3D* hist, const Point3D& point, bool interpolate)
{
	// check if outside
	if (point.x < hist->GetXaxis()->GetXmin() || point.x > hist->GetXaxis()->GetXmax() ||
		point.y < hist->GetYaxis()->GetXmin() || point.y > hist->GetYaxis()->GetXmax() ||
		point.z < hist->GetZaxis()->GetXmin() || point.z > hist->GetZaxis()->GetXmax())
	{
		return 0;
	}

	if (interpolate)
	{
		int numberBinsX = hist->GetXaxis()->GetNbins();
		int numberBinsY = hist->GetYaxis()->GetNbins();
		int numberBinsZ = hist->GetZaxis()->GetNbins();

		double x_min = hist->GetXaxis()->GetBinCenter(1);
		double x_max = hist->GetXaxis()->GetBinCenter(numberBinsX);
		double y_min = hist->GetYaxis()->GetBinCenter(1);
		double y_max = hist->GetYaxis()->GetBinCenter(numberBinsY);
		double z_min = hist->GetZaxis()->GetBinCenter(1);
		double z_max = hist->GetZaxis()->GetBinCenter(numberBinsZ);

		double x = std::clamp(point.x, x_min, x_max - 1e-8);
		double y = std::clamp(point.y, y_min, y_max - 1e-8);
		double z = std::clamp(point.z, z_min, z_max - 1e-8);

		return hist->Interpolate(x, y, z);
	}
	else
	{
		return hist->GetBinContent(hist->FindBin(point.x, point.y, point.z));
	}
}
