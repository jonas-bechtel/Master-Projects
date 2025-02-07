#pragma once

struct HeatMapData
{
	std::vector<double> values;
	int nRows;
	int nCols;
	double minValue;
	double maxValue;
	ImPlotPoint bottomLeft;
	ImPlotPoint topRight;

	void FromTH3D(TH3D* hist, float zSliceValue)
	{
		nRows = hist->GetNbinsX();
		nCols = hist->GetNbinsY();

		values.clear();
		values.reserve(nRows * nCols);

		int z_bin = hist->GetZaxis()->FindBin(zSliceValue);

		for (int x_bin = 1; x_bin <= nRows; x_bin++)
		{
			for (int y_bin = 1; y_bin <= nCols; y_bin++)
			{
				double content = hist->GetBinContent(x_bin, y_bin, z_bin);
				values.push_back(content);
			}
		}

		minValue = *std::min_element(values.begin(), values.end());
		maxValue = *std::max_element(values.begin(), values.end());

		bottomLeft = { hist->GetXaxis()->GetBinLowEdge(1), hist->GetYaxis()->GetBinLowEdge(1) };
		topRight = { hist->GetXaxis()->GetBinUpEdge(nRows), hist->GetYaxis()->GetBinUpEdge(nCols) };
	}
};
