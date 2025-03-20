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
		graph->AddPoint(detuningEnergies.at(i), value.at(i));
	}
}

double RateCoefficient::ConvolveFit(double Ed, double* csBins, const EnergyDistributionSet& set) const
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
		//std::cout << "no set for detuning energy " << detuningEnergy << std::endl;
		return 0.0;
	}
	const EnergyDistribution& distribution = set.distributions.at(index);

	for (int i = 0; i < distribution.psi.size(); i++)
	{
		sum += distribution.psi[i] * csBins[i] * csBins[i];
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
	while (std::getline(file, line))
	{
		std::vector<std::string> tokens = FileUtils::SplitLine(line, "\t");
		detuningEnergies.push_back(std::stod(tokens[0]));
		value.push_back(std::stod(tokens[1]));
		error.push_back(std::stod(tokens[3]));
		graph->AddPoint(std::stod(tokens[0]), std::stod(tokens[1]));
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

	outfile << "# Ed [eV]\tRelative rate\terror\tstatistical error\n";

	for (int i = 0; i < detuningEnergies.size(); i++)
	{
		outfile << detuningEnergies[i] << "\t" << value[i] << "\t" << 0.0 << "\t" << error[i] << "\n";
	}

	outfile.close();
}




