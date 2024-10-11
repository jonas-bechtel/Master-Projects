#pragma once

#include <filesystem>
#include <vector>

#include <TH3D.h>

class FileLoader
{
public:
	static FileLoader& getInstance() 
	{
		static FileLoader instance;
		return instance;
	}
	TH3D* LoadMatrixFile(std::filesystem::path filename);
	std::filesystem::path openFileExplorer(std::filesystem::path startPath);

private:
	std::vector<std::string> SplitLine(std::string& string, std::string& delimiter) const;
	TH3D* CreateTH3DfromHeader(std::ifstream& file) const;
	std::vector<double> CalculateBinEdges(const std::vector<double>& binCenters) const;

private:
	std::filesystem::path dataFilepath = "data/";
	int headerSize = 9;
	std::string xDelimiter = "\t";
	std::string zDelimiter = ";";
	
	mutable int matrixSize[3];// = { 100, 100, 100 };
};

