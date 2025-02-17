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
		nRows = hist->GetNbinsY();
		nCols = hist->GetNbinsX();
		//std::cout << nRows << ", " << nCols << std::endl;
		values.clear();
		values.reserve(nRows * nCols);

		int z_bin = hist->GetZaxis()->FindBin(zSliceValue);
		
		for (int y_bin = nRows; y_bin >= 1; y_bin--)
		{
			for (int x_bin = 1; x_bin <= nCols; x_bin++)
			{
				double content = hist->GetBinContent(x_bin, y_bin, z_bin);
				values.push_back(content);
			}
		}

		minValue = *std::min_element(values.begin(), values.end());
		maxValue = *std::max_element(values.begin(), values.end());

		bottomLeft = { hist->GetXaxis()->GetBinLowEdge(1), hist->GetYaxis()->GetBinLowEdge(1) };
		topRight = { hist->GetXaxis()->GetXmax(), hist->GetYaxis()->GetXmax()};
	}
};
