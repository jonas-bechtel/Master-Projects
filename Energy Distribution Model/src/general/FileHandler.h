#pragma once

#include "EnergyDistributionManager.h"

struct mbrcData;

class FileHandler
{
public:
	static FileHandler& GetInstance() 
	{
		static FileHandler instance;
		return instance;
	}
	std::filesystem::path SelectFile(const std::filesystem::path& startPath = "data\\", const std::vector<const char*>& filterPatterns = { "*.asc" });
	std::filesystem::path SelectFolder(const std::filesystem::path& startPath = "data\\");
	std::vector<std::filesystem::path> SelectFiles(const std::filesystem::path& startPath = "data\\");
	std::filesystem::path FindFileWithIndex(const std::filesystem::path& folder, int index);

	TH3D* LoadMatrixFile(const std::filesystem::path& filename);
	EnergyDistribution LoadEnergyDistribution(std::filesystem::path& filename, bool loadSamples);
	EnergyDistributionSet LoadEnergyDistributionSet(std::filesystem::path& folder);
	RateCoefficient LoadRateCoefficients(std::filesystem::path& filename);
	CrossSection LoadCrossSection(std::filesystem::path& filename);
	PlasmaRateCoefficient LoadPlasmaRate(std::filesystem::path& filename);

	int GetMaxIndex(std::filesystem::path energiesFile);
	std::array<float, 3> GetParamtersFromDescriptionFileAtIndex(const std::filesystem::path& descriptionFile, int index);

	void SaveEnergyDistributionSetAsHist(EnergyDistributionSet& eDistSet);
	void SaveEnergyDistributionSetAsSamples(EnergyDistributionSet& eDistSet);
	void SaveEnergyDistributionSetInfo(const EnergyDistributionSet& eDistSet);
	void SaveRateCoefficients(RateCoefficient& rc);
	void SaveCrossSection(CrossSection& cs);
	void SavePlasmaRate(PlasmaRateCoefficient& prc);

private:
	std::string GetHeaderFromFile(std::ifstream& file) const;
	std::vector<std::string> SplitLine(std::string& string, const std::string& delimiter) const;
	TH3D* CreateTH3DfromHeader(std::ifstream& file) const;
	std::vector<double> CalculateBinEdges(const std::vector<double>& binCenters, bool uniformDistances = true) const;
	EnergyDistribution CreateEnergyDistFromHeader(std::string& header);

private:
	std::filesystem::path dataFolder = "data/";
	std::filesystem::path outputFolder = "output/";
	int headerSize = 9;
	std::string xDelimiter = "\t";
	std::string zDelimiter = ";";
	
	// size will be overwritten when reading files
	mutable int matrixSize[3] = { 100, 100, 100 };
};

