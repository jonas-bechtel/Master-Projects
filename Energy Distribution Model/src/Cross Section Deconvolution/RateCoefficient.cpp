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

double RateCoefficient::ConvolveFit(double Ed, double* csBins, const EnergyDistributionSet& set, bool squareCS) const
{
	double sum = 0;

	// find correct distribution
	int index = GetIndexOfDetuningEnergy(Ed);
	if (index < 0)
	{
		std::cout << "did not find index of detuning energy " << Ed << std::endl;
		return 0.0;
	}

	if (index >= set.distributions.size())
	{
		std::cout << "index " << index << " is out of bounds for energy distribution set with size " << set.distributions.size() << std::endl;
		return 0.0;
	}
	const EnergyDistribution& distribution = set.distributions.at(index);

	for (int i = 0; i < distribution.psi.size(); i++)
	{
		if(squareCS)
			sum += distribution.psi[i] * csBins[i] * csBins[i];

		else
			sum += distribution.psi[i] * csBins[i];
	}
	//std::cout << "sum " << sum << "\n";
	return sum;
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

void RateCoefficient::Load(std::filesystem::path& filename)
{
	std::ifstream file(filename);

	// Check if the file was successfully opened
	if (!file.is_open())
	{
		std::cerr << "Error: Could not open the file " << filename << std::endl;
		return;
	}

	std::string line;
	// skip first line

	std::getline(file, line);
	int i = 0;
	while (std::getline(file, line))
	{
		std::vector<std::string> tokens = FileUtils::SplitLine(line, "\t");
		detuningEnergies.push_back(std::stod(tokens[0]));
		value.push_back(std::stod(tokens[1]));
		error.push_back(std::stod(tokens[3]));
		graph->SetPoint(i, std::stod(tokens[0]), std::stod(tokens[1]));
		graph->SetPointError(i, 0, std::stod(tokens[3]));
		i++;
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




