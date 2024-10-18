#pragma once
#include "EnergyDistributionModel.h"

#include <filesystem>
#include <vector>

#include <TH3D.h>

class FileHandler
{
public:
	static FileHandler& GetInstance() 
	{
		static FileHandler instance;
		return instance;
	}
	std::filesystem::path OpenFileExplorer(std::filesystem::path startPath = "data\\");
	std::filesystem::path FindFileWithIndex(std::filesystem::path folder, int index);

	TH3D* LoadMatrixFile(std::filesystem::path filename);
	//std::vector<EnergyDistribution> LoadEnergiesFile(std::filesystem::path filename);
	std::array<float, 3> GetParamtersFromDescriptionFileAtIndex(std::filesystem::path descriptionFile, int index);

	void SaveEnergyDistributionToFile(EnergyDistribution energyDistribution);

private:
	std::vector<std::string> SplitLine(std::string& string, std::string delimiter) const;
	TH3D* CreateTH3DfromHeader(std::ifstream& file) const;
	std::vector<double> CalculateBinEdges(const std::vector<double>& binCenters) const;

private:
	std::filesystem::path dataFolder = "data/";
	std::filesystem::path outputFolder = "output/";
	int headerSize = 9;
	std::string xDelimiter = "\t";
	std::string zDelimiter = ";";
	
	mutable int matrixSize[3];// = { 100, 100, 100 };
};

