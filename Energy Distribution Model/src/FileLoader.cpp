#include "FileLoader.h"

#include "tinyfiledialogs.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

TH3D* FileLoader::LoadMatrixFile(std::filesystem::path filename)
{
    std::ifstream file(dataFilepath / filename);

    // Check if the file was successfully opened
    if (!file.is_open()) 
    {
        std::cerr << "Error: Could not open the file " << dataFilepath / filename << std::endl;
        return nullptr;
    }

    TH3D* dataMatrix = CreateTH3DfromHeader(file);

    std::string line;
    int currentY = 0;
    int currentZ = 0;

    // Read the file line by line
    while (std::getline(file, line)) 
    {
        if (line.find(zDelimiter) != std::string::npos)
        {
            //std::cout << currentZ << "\n";
            currentZ++;
            currentY = 0;
            continue;
        }
        std::vector<std::string> tokens = SplitLine(line, xDelimiter);
        //std::cout << tokens.size() << "\n";
        for (int x = 0; x < tokens.size(); x++)
        {
            dataMatrix->SetBinContent(x, currentY, currentZ, std::stod(tokens[x]));
            //std::cout << "position " << x << " " << currentY << " " << currentZ << " value: " << tokens[x] << "\n";
        }
        currentY++;
    }

    // Close the file after reading
    file.close();
    
    std::cout << "loaded file: " << filename << "\n";
    return dataMatrix;
}

std::filesystem::path FileLoader::openFileExplorer(std::filesystem::path startPath = "data\\")
{
    const char* filterPatterns[] = { "*.asc" };
    const char* filePath = tinyfd_openFileDialog(
        "Choose a file",               // Dialog title
        startPath.string().c_str(),    // Default path or file
        1,                             // Number of filters
        filterPatterns,                // Filter patterns (NULL for any file type)
        NULL,                          // Filter description (optional)
        0                              // Allow multiple selection (0 = false)
    );
    if (!filePath)
    {
        return std::filesystem::path();
    }
    return std::filesystem::path(filePath);
}

std::vector<std::string> FileLoader::SplitLine(std::string& string, std::string& delimiter) const
{
    std::vector<std::string> tokens;
    tokens.reserve(matrixSize[0]);

    size_t pos = 0;
    std::string token;

    while ((pos = string.find(delimiter)) != std::string::npos)
    {
        token = string.substr(0, pos);
        tokens.push_back(token);
        string.erase(0, pos + delimiter.length());
    }
    tokens.push_back(string);
    
    return tokens;
}

TH3D* FileLoader::CreateTH3DfromHeader(std::ifstream& file) const
{
    std::string line;
    std::vector<double> xNodes, yNodes, zNodes;
    int lineCount = 0;
    
    while (lineCount < headerSize) 
    {
        std::getline(file, line);
        lineCount++;

        // comment lines starting with '#'
        if (line[0] == '#') 
        {
            // Check for specific comments to guide parsing
            if (line.find("#dim sizes x y z") != std::string::npos) 
            {
                // The next line will contain the dimension sizes
                file >> matrixSize[0] >> matrixSize[1] >> matrixSize[2];
                std::getline(file, line);  // To consume the end of the line
                lineCount++;

                xNodes.reserve(matrixSize[0]);
                yNodes.reserve(matrixSize[1]);
                zNodes.reserve(matrixSize[2]);
                //std::cout << matrixSize[0] << "\n";
            }
            else if (line.find("#x node positions") != std::string::npos)
            {
                // Read x node positions
                std::getline(file, line);
                lineCount++;
                std::stringstream ss(line);

                double value;
                while (ss >> value) 
                {
                    xNodes.push_back(value);
                }
            }
            else if (line.find("#y node positions") != std::string::npos) 
            {
                // Read y node positions
                std::getline(file, line);
                lineCount++;
                std::stringstream ss(line);
                double value;
                while (ss >> value) 
                {
                    yNodes.push_back(value);
                }
            }
            else if (line.find("#z node positions") != std::string::npos) 
            {
                // Read z node positions
                std::getline(file, line);
                lineCount++;
                std::stringstream ss(line);
                double value;
                while (ss >> value) 
                {
                    zNodes.push_back(value);
                }
            }
        }
    }
    std::vector<double> xBinEdges = CalculateBinEdges(xNodes);
    std::vector<double> yBinEdges = CalculateBinEdges(yNodes);
    std::vector<double> zBinEdges = CalculateBinEdges(zNodes);

    TH3D* dataMatrix = new TH3D("name", "title", matrixSize[0], xBinEdges.data(),
                                                 matrixSize[1], yBinEdges.data(),
                                                 matrixSize[2], zBinEdges.data());
    dataMatrix->SetXTitle("X-axis");
    dataMatrix->SetYTitle("Y-axis");
    dataMatrix->SetZTitle("Z-axis");
    
    return dataMatrix;
}

std::vector<double> FileLoader::CalculateBinEdges(const std::vector<double>& binCenters) const
{
    std::vector<double> binEdges;
    if (binCenters.empty()) 
    {
        return binEdges;
    }

    int nBins = binCenters.size();
    binEdges.reserve(nBins + 1);

    // Compute first edge by extrapolation
    double firstEdge = binCenters[0] - (binCenters[1] - binCenters[0]) / 2.0;
    binEdges.push_back(firstEdge);

    // Compute middle edges as averages of adjacent bin centers
    for (int i = 0; i < nBins - 1; ++i) 
    {
        double edge = (binCenters[i] + binCenters[i + 1]) / 2.0;
        binEdges.push_back(edge);
    }

    // Compute last edge by extrapolation
    double lastEdge = binCenters[nBins - 1] + (binCenters[nBins - 1] - binCenters[nBins - 2]) / 2.0;
    binEdges.push_back(lastEdge);

    return binEdges;
}
