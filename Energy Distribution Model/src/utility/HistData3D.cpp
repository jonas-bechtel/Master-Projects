#include "pch.h"
#include "HistData3D.h"

#include "ElectronBeam.h"

int HistData3D::rebinningFactors[3] = { 10, 10, 10 };

HistData3D::HistData3D(TH3D* hist)
{
	if(!hist)
	{
		fullHistogram = nullptr;
		return;
	}

	delete fullHistogram;
	fullHistogram = hist;

	xAxis.reserve(fullHistogram->GetNbinsX());
	yAxis.reserve(fullHistogram->GetNbinsY());
	zAxis.reserve(fullHistogram->GetNbinsZ());

	projectionValuesX.reserve(fullHistogram->GetNbinsX());
	projectionValuesY.reserve(fullHistogram->GetNbinsY());
	projectionValuesZ.reserve(fullHistogram->GetNbinsZ());

	centerValue.reserve(fullHistogram->GetNbinsZ());
	outsideValue.reserve(fullHistogram->GetNbinsZ());

	for (int i = 1; i <= fullHistogram->GetNbinsX(); i++)
	{
		xAxis.push_back(fullHistogram->GetXaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= fullHistogram->GetNbinsY(); i++)
	{
		yAxis.push_back(fullHistogram->GetYaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= fullHistogram->GetNbinsZ(); i++)
	{
		zAxis.push_back(fullHistogram->GetZaxis()->GetBinCenter(i));
	}

	UpdateData();
}

HistData3D::HistData3D(HistData3D&& other) noexcept
{
	fullHistogram = other.fullHistogram;
	other.fullHistogram = nullptr;

	xAxis = std::move(other.xAxis);
	yAxis = std::move(other.yAxis);
	zAxis = std::move(other.zAxis);

	projectionValuesX = std::move(other.projectionValuesX);
	projectionValuesY = std::move(other.projectionValuesY);
	projectionValuesZ = std::move(other.projectionValuesZ);

	centerValue = std::move(other.centerValue);
	outsideValue = std::move(other.outsideValue);

	slice = std::move(other.slice);

	label = std::move(other.label);
}

HistData3D& HistData3D::operator=(HistData3D&& other) noexcept
{
	if (this == &other) return *this;

	delete fullHistogram;
	fullHistogram = other.fullHistogram;
	other.fullHistogram = nullptr;

	xAxis = std::move(other.xAxis);
	yAxis = std::move(other.yAxis);
	zAxis = std::move(other.zAxis);

	projectionValuesX = std::move(other.projectionValuesX);
	projectionValuesY = std::move(other.projectionValuesY);
	projectionValuesZ = std::move(other.projectionValuesZ);

	centerValue = std::move(other.centerValue);
	outsideValue = std::move(other.outsideValue);

	slice = std::move(other.slice);

	label = std::move(other.label);

	return *this;
}

HistData3D::~HistData3D()
{
	delete fullHistogramSmall;
	delete fullHistogram;
}

TH3D* HistData3D::GetHist() const
{
	return fullHistogram;
}

std::string HistData3D::GetLabel()
{
	return label;
}

void HistData3D::SetLabel(std::string str)
{
	label = str;
}

void HistData3D::UpdateSlice(float zValue)
{
	slice.FromTH3D(fullHistogram, zValue);
}

void HistData3D::UpdateData()
{
	TH1D* projectionX = fullHistogram->ProjectionX();
	TH1D* projectionY = fullHistogram->ProjectionY();
	TH1D* projectionZ = fullHistogram->ProjectionZ();
	projectionX->Scale(1.0 / projectionX->Integral());
	projectionY->Scale(1.0 / projectionY->Integral());
	projectionZ->Scale(1.0 / projectionZ->Integral());

	for (int i = 1; i <= projectionX->GetNbinsX(); i++)
	{
		projectionValuesX.push_back(projectionX->GetBinContent(i));
	}
	for (int i = 1; i <= projectionY->GetNbinsX(); i++)
	{
		projectionValuesY.push_back(projectionY->GetBinContent(i));
	}
	for (int i = 1; i <= projectionZ->GetNbinsX(); i++)
	{
		projectionValuesZ.push_back(projectionZ->GetBinContent(i));
	}

	delete projectionX;
	delete projectionY;
	delete projectionZ;

	int binInCenterX = fullHistogram->GetNbinsX() / 2;
	for (int i = 1; i <= fullHistogram->GetNbinsZ(); i++)
	{
		double zValue = fullHistogram->GetZaxis()->GetBinCenter(i);
		int binInCenterY = fullHistogram->GetYaxis()->FindBin(ElectronBeam::Trajectory(zValue));

		double energyValueIn = fullHistogram->GetBinContent(binInCenterX, binInCenterY, i);
		double energyValueOut = fullHistogram->GetBinContent(1, 1, i);

		centerValue.push_back(energyValueIn);
		outsideValue.push_back(energyValueOut);
	}

	//// quick hack to save the center values to a file
	//std::ofstream outFile("potential_values.txt");
	//for (size_t i = 0; i < centerValue.size(); i++)
	//{
	//	outFile << zAxis[i] << "\t" << centerValue[i] << "\t" << outsideValue[i] << "\n";
	//}
	//outFile.close();
}

void HistData3D::Plot3D(ROOTCanvas* canvas, int pos)
{
	if (!canvas->IsShown() || !fullHistogram) return;

	canvas->cd(pos);

	delete fullHistogramSmall;
	fullHistogramSmall = (TH3D*)fullHistogram->Rebin3D(rebinningFactors[0],
		rebinningFactors[1],
		rebinningFactors[2], "dist small");

	if (!fullHistogramSmall) return;

	fullHistogramSmall->GetXaxis()->SetTitle("x-axis");
	fullHistogramSmall->GetYaxis()->SetTitle("y-axis");
	fullHistogramSmall->GetZaxis()->SetTitle("z-axis");
	fullHistogramSmall->Draw("BOX2 COLZ");
}

void HistData3D::PlotSlice() const
{
	slice.Plot(label);
}

void HistData3D::PlotProjectionX(ImPlotLineFlags_ flags) const
{
	ImPlot::PlotLine(label.c_str(), xAxis.data(), projectionValuesX.data(), xAxis.size());
}

void HistData3D::PlotProjectionY(ImPlotLineFlags_ flags) const
{
	ImPlot::PlotLine(label.c_str(), yAxis.data(), projectionValuesY.data(), yAxis.size());
}

void HistData3D::PlotProjectionZ(ImPlotLineFlags_ flags) const
{
	ImPlot::PlotLine(label.c_str(), zAxis.data(), projectionValuesZ.data(), zAxis.size());
}

void HistData3D::PlotInsideOutsideValue() const
{
	ImPlot::PlotLine(label.c_str(), zAxis.data(), centerValue.data(), zAxis.size());
	ImPlot::PlotLine(label.c_str(), zAxis.data(), outsideValue.data(), zAxis.size(), ImPlotLineFlags_Segments);
}

bool HistData3D::ShowRebinningFactorsInput()
{
	ImGui::SetNextItemWidth(120.0f);
	return ImGui::InputInt3("Rebinning factors", rebinningFactors);
}

