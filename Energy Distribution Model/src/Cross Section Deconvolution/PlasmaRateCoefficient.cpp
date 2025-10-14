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

void PlasmaRateCoefficient::ConvolveFromErrorIterationArray(const CrossSection& cs)
{
	label = "plasma from " + cs.GetLabel();

	double factor = TMath::Power((endTemperature / startTemperature), (1.0 / (numberValues)));

	int errorIterations = cs.valueArray.size() / cs.hist->GetNbinsX();
	//std::cout << "error iterations: " << errorIterations << std::endl;
	if (errorIterations < 1)
	{
		std::cerr << "Error: Not enough data in valueArray to perform convolution." << std::endl;
		return;
	}

	temperatures.clear();
	values.clear();
	temperatures.reserve(numberValues);
	values.reserve(numberValues);

	std::vector<double> plasmaValueArray;
	plasmaValueArray.reserve(numberValues * errorIterations);

	for (double T = startTemperature; T <= endTemperature; T *= factor)
	{
		temperatures.push_back(T);
		for (int j = 0; j < errorIterations; j++)
		{
			plasmaValueArray.push_back(0);
			for (int i = 0; i < cs.energies.size(); i++)
			{
				double energy = cs.energies.at(i);
				//std::cout << "energy: " << energy << std::endl;
				double csValue = cs.valueArray.at(i * errorIterations + j) * cs.hist->GetBinWidth(i + 1);
				//std::cout << "csValue: " << csValue << std::endl;
				double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
				velocity *= 100; // convert to cm/s
				//std::cout << "velocity: " << velocity << std::endl;
				double f_pl = BoltzmannDistribution::Function(energy, T);
				//std::cout << "f_pl: " << f_pl << std::endl;
				plasmaValueArray.back() += csValue * velocity * f_pl;
			}
			
		}
	}
	//std::cout << "plasmaValueArray size: " << plasmaValueArray.size() << std::endl;
	
	// Calculate mean and errors
	for (int j = 0; j < numberValues; j++)
	{
		double mean = 0;
		double error = 0;

		for (int i = 0; i < errorIterations; i++)
		{
			//std::cout << "i: " << i << std::endl;
			double value = plasmaValueArray.at(i + j * errorIterations);
			mean += value;
		}
		mean /= errorIterations;

		for (int i = 0; i < errorIterations; i++)
		{
			error += pow(plasmaValueArray[i + j * errorIterations] - mean, 2);
		}
		error = sqrt(error / (errorIterations - 1));

		values.push_back(mean);
		errors.push_back(error);
		//std::cout << "j: " << j << std::endl;
		//std::cout << "T: " << T_start * pow(factor, j) << ", mean: " << mean << ", error: " << error << std::endl;
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

void PlasmaRateCoefficient::ShowConvolutionParamterInputs()
{
	ImGui::PushItemWidth(100.0f);
	ImGui::InputInt("number of values", &numberValues);
	ImGui::InputDouble("start T [K]", &startTemperature);
	ImGui::InputDouble("end T [K]", &endTemperature);
	ImGui::PopItemWidth();
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
