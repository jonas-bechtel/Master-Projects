#pragma once

#include "EnergyDistributionWindow.h"

namespace FileUtils
{
	std::filesystem::path GetDataFolder();
	std::filesystem::path GetMeasuredRateCoefficientFolder();
	std::filesystem::path GetOutputFolder();
	std::filesystem::path GetPlasmaRateFolder();
	std::filesystem::path GetCrossSectionFolder();
	std::filesystem::path GetRateCoefficientFolder();
	std::filesystem::path GetEnergyDistSetFolder();
	std::filesystem::path GetCoolingForceCurveFolder();
	std::filesystem::path GetNumericalCoolingForceCurveFolder();
	std::filesystem::path GetMeasuredCoolingForceCurveFolder();

	std::filesystem::path SelectFile(const std::filesystem::path& startPath = "data\\", const std::vector<const char*>& filterPatterns = { "*.asc" });
	std::filesystem::path SelectFolder(const std::filesystem::path& startPath = "data\\");
	std::vector<std::filesystem::path> SelectFiles(const std::filesystem::path& startPath = "data\\", const std::vector<const char*>& filterPatterns = { "*.asc" });
	
	TH3D* LoadMatrixFile(const std::filesystem::path& filename);
	TH3D* CreateTH3DfromHeader(std::ifstream& file);

	std::string GetHeaderFromFile(std::ifstream& file);
	std::vector<std::string> SplitLine(std::string& line, const std::string& delimiter);
	std::filesystem::path FindFileWithIndex(const std::filesystem::path& folder, int index);
	int GetMaxIndex(std::filesystem::path energiesFile);
	std::array<float, 3> GetParametersFromDescriptionFileAtIndex(const std::filesystem::path& descriptionFile, int index);
	std::vector<double> CalculateBinEdges(const std::vector<double>& binCenters, bool uniformDistances = true);
}


