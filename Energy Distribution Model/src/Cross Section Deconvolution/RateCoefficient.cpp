#include "pch.h"
#include "RateCoefficient.h"
#include "EnergyDistributionSet.h"
#include "CrossSection.h"

#include "FileUtils.h"

RateCoefficient::RateCoefficient()
{
}

int RateCoefficient::GetIndexOfDetuningEnergy(double Ed) const
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

void RateCoefficient::VaryGraphValues()
{
	// Create a random number generator and normal distribution
	static std::mt19937 rng(std::random_device{}());

	for (int i = 0; i < graph->GetN(); i++)
	{
		// Use error.at(i) as the standard deviation if available, otherwise set to 1.0
		double mean = value.at(i);
		double stddev = (i < error.size()) ? error.at(i) : 1.0;
		std::normal_distribution<double> dist(mean, stddev);

		double variedValue = dist(rng);
		graph->SetPointY(i, variedValue);
	}
}

void RateCoefficient::ResetGraphValues()
{
	for (int i = 0; i < graph->GetN(); i++)
	{
		graph->SetPointY(i, value.at(i));
	}
}

void RateCoefficient::SetLabel(std::string label)
{
	this->label = label;
}

std::string RateCoefficient::GetLabel()
{
	return label;
}

void RateCoefficient::Convolve(const CrossSection& cs, EnergyDistributionSet& set)
{
	energyDistriubtionSetFolder = set.Label();
	crossSectionFile = cs.label;
	measured = false;

	set.CalculatePsisFromBinning(cs.hist);

	for (const EnergyDistribution& eDist : set.distributions)
	{
		detuningEnergies.push_back(eDist.eBeamParameter.detuningEnergy.get());
	}

	value.clear();
	error.clear();
	value.resize(set.distributions.size());
	error.resize(set.distributions.size());
	psiSubfunctions.resize(cs.hist->GetNbinsX());
	for (int i = 0; i < cs.hist->GetNbinsX(); i++)
	{
		psiSubfunctions[i].resize(set.distributions.size());
	}

	int j = 0;
	for (const EnergyDistribution& eDist : set.distributions)
	{
		for (int i = 0; i < eDist.psi.size(); i++)
		{
			double subValue = eDist.psi[i] * cs.values[i];
			value[j] += subValue;
			psiSubfunctions[i][j] = subValue;
		}
		j++;
	}
	graph->Clear();
	for (int i = 0; i < detuningEnergies.size(); i++)
	{
		graph->SetPoint(i, detuningEnergies.at(i), value.at(i));
		graph->SetPointError(i, 0, error.at(i));
	}
}

void RateCoefficient::Plot(bool showMarkers) const
{
	if (showMarkers) ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
	ImPlot::PlotLine(label.c_str(), detuningEnergies.data(), value.data(), value.size());
	ImPlot::PlotErrorBars(label.c_str(), detuningEnergies.data(), value.data(), error.data(), error.size());
}

void RateCoefficient::PlotSubfunctions() const
{
	ImPlot::PushColormap("Jet");
	int j = 0;
	int numberSubFunc = psiSubfunctions.size();
	for (const std::vector<double>& subFunction : psiSubfunctions)
	{
		ImPlot::PushStyleColor(ImPlotCol_Line, ImPlot::SampleColormap(j / (float)numberSubFunc));
		ImPlot::PlotLine(label.c_str(), detuningEnergies.data(), subFunction.data(), subFunction.size());
		ImPlot::PopStyleColor();
		j++;
	}
	ImPlot::PopColormap();
}

void RateCoefficient::Clear()
{
	graph->Clear();
	detuningEnergies.clear();
	value.clear();
	error.clear();
	psiSubfunctions.clear();
	measured = false;
	label = "mbrc";
	energyDistriubtionSetFolder.clear();
	crossSectionFile.clear();
}

void RateCoefficient::Load(std::filesystem::path& filename)
{
	std::ifstream file(filename);

	// Check if the file was successfully opened
	if (!file.is_open())
	{
		std::cerr << "Error: Could not open the file " << filename << std::endl;
		return;
	}

	Clear();

	std::string header = FileUtils::GetHeaderFromFile(file);

	std::string line;
	while (std::getline(file, line))
	{
		std::vector<std::string> tokens = FileUtils::SplitLine(line, "\t");

		if (tokens.size() < 4)
		{
			std::cerr << "Error: Invalid line format in file " << filename << std::endl;
			Clear();
			return;
		}

		detuningEnergies.push_back(std::stod(tokens[0]));
		value.push_back(std::stod(tokens[1]));
		error.push_back(std::stod(tokens[3]));
	}

	// sort values in ascending order of detuning energy
	SortValuesByDetuningEnergy();

	// put values into graph
	for (size_t i = 0; i < detuningEnergies.size(); i++)
	{
		graph->SetPoint(i, detuningEnergies.at(i), value.at(i));
		graph->SetPointError(i, 0, error.at(i));
	}

	measured = true;
	label = filename.filename().string();
}

void RateCoefficient::Save() const
{
	// set the output filepath
	std::filesystem::path file = FileUtils::GetRateCoefficientFolder() / (label + ".dat");

	// Create the directories if they don't exist
	if (!std::filesystem::exists(file.parent_path()))
	{
		std::filesystem::create_directories(file.parent_path());
	}

	std::ofstream outfile(file);

	if (!outfile.is_open())
	{
		std::cerr << "Error opening file" << std::endl;
		return;
	}

	outfile << "# Ed [eV]\tRate Coefficient\terror\tstatistical error\n";

	for (int i = 0; i < detuningEnergies.size(); i++)
	{
		outfile << detuningEnergies[i] << "\t" << value[i] << "\t" << 0.0 << "\t" << error[i] << "\n";
	}

	outfile.close();
}

void RateCoefficient::SortValuesByDetuningEnergy()
{
	// Create a vector of indices
	std::vector<size_t> idx(detuningEnergies.size());
	for (size_t i = 0; i < idx.size(); ++i) 
		idx[i] = i;

	// Sort indices based on detuningEnergies
	std::sort(idx.begin(), idx.end(), [this](size_t a, size_t b) {
		return detuningEnergies[a] < detuningEnergies[b];
		});

	// Create sorted copies
	std::vector<double> sortedDetuningEnergies, sortedValue, sortedError;
	sortedDetuningEnergies.reserve(idx.size());
	sortedValue.reserve(idx.size());
	sortedError.reserve(idx.size());

	for (size_t i : idx) {
		sortedDetuningEnergies.push_back(detuningEnergies[i]);
		sortedValue.push_back(value[i]);
		sortedError.push_back(error[i]);
	}

	// Assign sorted vectors back
	detuningEnergies = std::move(sortedDetuningEnergies);
	value = std::move(sortedValue);
	error = std::move(sortedError);
}


