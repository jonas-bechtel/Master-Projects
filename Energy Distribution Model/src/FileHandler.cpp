#include "FileHandler.h"

#include "tinyfiledialogs.h"

#include <iostream>
#include <array>
#include <fstream>
#include <string>
#include <sstream>

TH3D* FileHandler::LoadMatrixFile(std::filesystem::path filename)
{
    std::ifstream file(filename);

    // Check if the file was successfully opened
    if (!file.is_open()) 
    {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
        return nullptr;
    }

    TH3D* dataMatrix = CreateTH3DfromHeader(file);

    std::string line;
    int currentY = 1;
    int currentZ = 1;

    // Read the file line by line
    while (std::getline(file, line)) 
    {
        if (line.find(zDelimiter) != std::string::npos)
        {
            currentZ++;
            currentY = 1;
            continue;
        }
        std::vector<std::string> tokens = SplitLine(line, xDelimiter);

        for (int x = 1; x <= tokens.size(); x++)
        {
            dataMatrix->SetBinContent(x, currentY, currentZ, std::stod(tokens[x-1]));
        }
        currentY++;
    }

    // Close the file after reading
    file.close();
    
    std::cout << "loaded file: " << filename.parent_path().parent_path().filename() /
                                    filename.parent_path().filename() /
                                    filename.filename() << "\n";
    return dataMatrix;
}

//std::vector<EnergyDistribution> FileHandler::LoadEnergiesFile(std::filesystem::path filename)
//{
//    std::ifstream file(filename);
//
//    std::vector<EnergyDistribution> energyDistributions;
//
//    // Check if the file was successfully opened
//    if (!file.is_open())
//    {
//        std::cerr << "Error: Could not open the file " << filename << std::endl;
//        return energyDistributions;
//    }
//
//    std::filesystem::path baseDirectory = filename.parent_path();
//    std::filesystem::path densityDirectory = baseDirectory / "e-densities";
//    std::filesystem::path energyDirectory = baseDirectory / "lab-energies";
//
//    std::string line;
//
//    // trow away first line
//    std::getline(file, line);
//
//    // Read the file line by line
//    while (std::getline(file, line))
//    {
//        std::vector<std::string> tokens = SplitLine(line, "\t");
//        EnergyDistribution distribution;
//        distribution.index = std::stoi(tokens[0]);
//        distribution.driftTubeVoltage = std::stod(tokens[1]);
//        distribution.centerLabEnergy = std::stod(tokens[3]);
//
//        // look through density dir to find matching index
//        for (const auto& file : std::filesystem::directory_iterator(densityDirectory))
//        {
//            // Check if it's a regular file
//            if (std::filesystem::is_regular_file(file.status()))
//            {
//                // compare first 4 characters to index
//                if (std::stoi(file.path().filename().string().substr(0, 4)) == distribution.index)
//                {
//                    distribution.densityFile = file.path();
//                    break;
//                }
//            }
//        }
//        // look through energy dir to find matching index
//        for (const auto& file : std::filesystem::directory_iterator(energyDirectory))
//        {
//            // Check if it's a regular file
//            if (std::filesystem::is_regular_file(file.status()))
//            {
//                // compare first 4 characters to index
//                if (std::stoi(file.path().filename().string().substr(0, 4)) == distribution.index)
//                {
//                    distribution.energyFile = file.path();
//                    break;
//                }
//            }
//        }
//        energyDistributions.push_back(distribution);
//    }
//
//    return energyDistributions;
//}

int FileHandler::GetMaxIndex(std::filesystem::path energiesFile)
{
    std::ifstream file(energiesFile);

    // Check if the file was successfully opened
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open the file " << energiesFile << std::endl;
    }

    int maxIndex = 0;
    std::string line;

    // throw away first line
    std::getline(file, line);

    while (std::getline(file, line))
    {
        maxIndex = std::max(maxIndex, std::stoi(line.substr(0, 4)));
    }

    return maxIndex;
}

std::array<float, 3> FileHandler::GetParamtersFromDescriptionFileAtIndex(std::filesystem::path descriptionFile, int index)
{
    std::array<float, 3> parameter = { 0, 0, 0 };

    std::ifstream file(descriptionFile);

    // Check if the file was successfully opened
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open the file " << descriptionFile << std::endl;
        return parameter;
    }

    std::string line;

    // throw away first line
    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::vector<std::string> tokens = SplitLine(line, "\t");
        if (std::stoi(tokens[0]) == index)
        {
            parameter[0] = std::stof(tokens[1]);
            parameter[1] = std::stof(tokens[2]);
            parameter[2] = std::stof(tokens[3]);
            return parameter;
        }
    }

    std::cout << "Index " << index << " not found in file: " << descriptionFile << std::endl;
    return parameter;
}

std::filesystem::path FileHandler::OpenFileExplorer(std::filesystem::path startPath)
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

std::filesystem::path FileHandler::FindFileWithIndex(std::filesystem::path folder, int index)
{
    // look through dir to find matching index
    for (const auto& file : std::filesystem::directory_iterator(folder))
    {
        // Check if it's a regular file
        if (std::filesystem::is_regular_file(file.status()))
        {
            // compare first 4 characters to index
            if (std::stoi(file.path().filename().string().substr(0, 4)) == index)
            {
                return file.path();
            }
        }
    }
    std::cout << "No file with index: " << index << " was found in " << folder << std::endl;
    return std::filesystem::path();
}

void FileHandler::SaveEnergyDistributionToFile(EnergyDistribution energyDistribution)
{
    std::filesystem::path file = outputFolder.string() / energyDistribution.folder.filename() /
        std::filesystem::path(energyDistribution.Filename());

    std::filesystem::path dir = std::filesystem::path(file).parent_path();

    // Create the directories if they don't exist
    if (!std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
    }

    std::ofstream outfile(file);
    
    if (!outfile.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return ;
    }
    outfile << energyDistribution.String();

    outfile << "# bin center [eV]\tbin value\n";
    for (int i = 1; i <= energyDistribution->GetNbinsX(); i++)
    {
        outfile << energyDistribution->GetBinCenter(i) << "\t";
        outfile << energyDistribution->GetBinContent(i) << "\n";
    }
    
    outfile.close();
}

std::vector<std::string> FileHandler::SplitLine(std::string& string, std::string delimiter) const
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

TH3D* FileHandler::CreateTH3DfromHeader(std::ifstream& file) const
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

std::vector<double> FileHandler::CalculateBinEdges(const std::vector<double>& binCenters) const
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
