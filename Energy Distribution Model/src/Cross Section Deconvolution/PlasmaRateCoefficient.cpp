#include "pch.h"
#include "PlasmaRateCoefficient.h"
#include "Constants.h"
#include "BoltzmannDistribution.h"

#include "FileUtils.h"

std::string PlasmaRateCoefficient::GetLabel()
{
	return label;
}

void PlasmaRateCoefficient::SetLabel(std::string label)
{
	this->label = label;
}

void PlasmaRateCoefficient::Convolve(const CrossSection& cs)
{
	label = "plasma from " + cs.GetLabel();

	int numberValue = 10000;
	double T_start = 1;
	double T_end = 50000;
	double step = (T_end - T_start) / numberValue;

	temperatures.clear();
	values.clear();
	temperatures.reserve(numberValue);
	values.reserve(numberValue);

	for (double T = 1; T <= T_end; T += step)
	{
		temperatures.push_back(T);
		values.push_back(0);
		for (int i = 1; i <= cs.hist->GetNbinsX(); i++)
		{
			double energy = cs.hist->GetBinCenter(i);
			//std::cout << "energy: " << energy << std::endl;
			double csValue = cs.hist->GetBinContent(i) * cs.hist->GetBinWidth(i);
			//std::cout << "csValue: " << csValue << std::endl;
			double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
			//std::cout << "velocity: " << velocity << std::endl;
			double f_pl = BoltzmannDistribution::Function(energy, T);
			//std::cout << "f_pl: " << f_pl << std::endl;
			values.back() += csValue * velocity * f_pl;
		}
		errors.push_back(0);
	}
}

void PlasmaRateCoefficient::Plot(bool showMarkers) const
{
	if (showMarkers) ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
	ImPlot::PlotLine(label.c_str(), temperatures.data(), values.data(), values.size());
	ImPlot::PlotErrorBars(label.c_str(), temperatures.data(), values.data(), errors.data(), errors.size());

}

void PlasmaRateCoefficient::Load(std::filesystem::path file)
{
	std::ifstream infile(file);

	// Check if the file was successfully opened
	if (!infile.is_open())
	{
		std::cerr << "Error: Could not open the file " << file << std::endl;
		return;
	}

	std::string line;
	// skip first line
	std::getline(infile, line);

	while (std::getline(infile, line))
	{
		std::vector<std::string> tokens = FileUtils::SplitLine(line, "\t");
		temperatures.push_back(std::stod(tokens[0]));
		values.push_back(std::stod(tokens[1]));
		errors.push_back(std::stod(tokens[2]));
	}
}

void PlasmaRateCoefficient::Save() const
{
	// set the output filepath
	std::filesystem::path filepath = FileUtils::GetPlasmaRateFolder() / (label + ".dat");

	// Create the directories if they don't exist
	if (!std::filesystem::exists(filepath.parent_path()))
	{
		std::filesystem::create_directories(filepath.parent_path());
	}

	std::ofstream outfile(filepath);

	if (!outfile.is_open())
	{
		std::cerr << "Error opening file" << std::endl;
		return;
	}

	outfile << "# Temperature [K]\tPlasma rate coefficient\terror\n";

	for (int i = 0; i < temperatures.size(); i++)
	{
		outfile << temperatures[i] << "\t" << values[i] << "\t" << errors[i] << "\n";
	}

	outfile.close();
}
