#pragma once

#include "EnergyDistributionManager.h"

class FileHandler
{
public:
	static FileHandler& GetInstance() 
	{
		static FileHandler instance;
		return instance;
	}
	std::filesystem::path SelectFile(const std::filesystem::path& startPath = "data\\");
	std::vector<std::filesystem::path> SelectFiles(const std::filesystem::path& startPath = "data\\");
	std::filesystem::path FindFileWithIndex(const std::filesystem::path& folder, int index);

	TH3D* LoadMatrixFile(const std::filesystem::path& filename);
	//EnergyDistribution* LoadEnergyDistributionSamples(const std::filesystem::path& filename);
	EnergyDistribution* LoadEnergyDistribution(std::filesystem::path& filename, bool loadSamples);

	int GetMaxIndex(std::filesystem::path energiesFile);
	std::array<float, 3> GetParamtersFromDescriptionFileAtIndex(const std::filesystem::path& descriptionFile, int index);

	void SaveEnergyDistributionHistToFile(EnergyDistribution* energyDistribution);
	void SaveEnergyDistributionSamplesToFile(EnergyDistribution* energyDistribution);

private:
	std::string GetHeaderFromFile(std::ifstream& file) const;
	std::vector<std::string> SplitLine(std::string& string, const std::string& delimiter) const;
	TH3D* CreateTH3DfromHeader(std::ifstream& file) const;
	std::vector<double> CalculateBinEdges(const std::vector<double>& binCenters) const;
	EnergyDistribution* CreateEnergyDistFromHeader(std::string& header);

private:
	std::filesystem::path dataFolder = "data/";
	std::filesystem::path outputFolder = "output/";
	int headerSize = 9;
	std::string xDelimiter = "\t";
	std::string zDelimiter = ";";
	
	// size will be overwritten when reading files
	mutable int matrixSize[3] = { 100, 100, 100 };
};

