#include "pch.h"

#include "tinyfiledialogs.h"

#include "RateCoefficient.h"
#include "CrossSection.h"
#include "PlasmaRateCoefficient.h"
#include "FileUtils.h"


namespace FileUtils
{
    static std::filesystem::path dataFolder = "data\\";
    static std::filesystem::path measuredRateCoefficientFolder = dataFolder / "Rate Coefficients\\";

    static std::filesystem::path outputFolder = "output\\";
    static std::filesystem::path plasmaRateFolder = outputFolder / "Plasma Rate Coefficients\\";
    static std::filesystem::path crossSectionFolder = outputFolder / "Cross Sections\\";
    static std::filesystem::path rateCoefficientFolder = outputFolder / "Rate Coefficients\\";
    static std::filesystem::path energyDistSetFolder = outputFolder / "Energy Distribution Sets\\";
    static std::filesystem::path coolingForceCurveFolder = outputFolder / "Cooling Force Curves\\";
    static std::filesystem::path numericalCoolingForceCurveFolder = coolingForceCurveFolder / "Numerical\\";
    static int headerSize = 9;
    static std::string xDelimiter = "\t";
    static std::string zDelimiter = ";";

    // size will be overwritten when reading files
    static int matrixSize[3] = { 100, 100, 100 };


    std::filesystem::path GetDataFolder()
    {
        return dataFolder;
    }

    std::filesystem::path GetMeasuredRateCoefficientFolder()
    {
        return measuredRateCoefficientFolder;
    }

    std::filesystem::path GetOutputFolder()
    {
        return outputFolder;
    }

    std::filesystem::path GetPlasmaRateFolder()
    {
        return plasmaRateFolder;
    }

    std::filesystem::path GetCrossSectionFolder()
    {
        return crossSectionFolder;
    }

    std::filesystem::path GetRateCoefficientFolder()
    {
        return rateCoefficientFolder;
    }

    std::filesystem::path GetEnergyDistSetFolder()
    {
        return energyDistSetFolder;
    }

    std::filesystem::path GetCoolingForceCurveFolder()
    {
        return coolingForceCurveFolder;
    }

    std::filesystem::path GetNumericalCoolingForceCurveFolder()
    {
        return numericalCoolingForceCurveFolder;
    }

    std::filesystem::path SelectFile(const std::filesystem::path& startPath, const std::vector<const char*>& filterPatterns)
    {
        //const char* filterPatterns[] = { "*.asc" };
        const char* filePath = tinyfd_openFileDialog(
            "Choose a file",               // Dialog title
            startPath.string().c_str(),    // Default path or file
            1,                             // Number of filters
            filterPatterns.data(),                // Filter patterns (NULL for any file type)
            NULL,                          // Filter description (optional)
            0                              // Allow multiple selection (0 = false)
        );
        if (!filePath)
        {
            return std::filesystem::path();
        }
        return std::filesystem::path(filePath);
    }

    std::filesystem::path SelectFolder(const std::filesystem::path& startPath)
    {
        std::filesystem::path projectDir = std::filesystem::current_path();
        std::filesystem::path fullStartPath = projectDir / startPath;

        const char* folder = tinyfd_selectFolderDialog("Choose a folder", fullStartPath.string().c_str());

        if (!folder)
        {
            return std::filesystem::path();
        }
        return std::filesystem::path(folder);
    }

    std::vector<std::filesystem::path> SelectFiles(const std::filesystem::path& startPath, const std::vector<const char*>& filterPatterns)
    {
        const char* filePaths = tinyfd_openFileDialog(
            "Choose a file",               // Dialog title
            startPath.string().c_str(),    // Default path or file
            1,                             // Number of filters
            filterPatterns.data(),         // Filter patterns (NULL for any file type)
            NULL,                          // Filter description (optional)
            1                              // Allow multiple selection (0 = false)
        );
        if (!filePaths)
        {
            return std::vector<std::filesystem::path>();
        }
        std::string paths1 = std::string(filePaths);
        std::vector<std::string> paths = SplitLine(paths1, "|");
        std::vector<std::filesystem::path> splitPaths(paths.begin(), paths.end());
        return splitPaths;
    }

    TH3D* LoadMatrixFile(const std::filesystem::path& filename)
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
                dataMatrix->SetBinContent(x, currentY, currentZ, std::stod(tokens[x - 1]));
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

    TH3D* CreateTH3DfromHeader(std::ifstream& file)
    {
        std::string line;
        std::vector<double> xNodes, yNodes, zNodes;

        while (file.peek() == '#')
        {
            std::getline(file, line);

            // comment lines starting with '#'
            if (line[0] == '#')
            {
                // Check for specific comments to guide parsing
                if (line.find("#dim sizes x y z") != std::string::npos)
                {
                    // The next line will contain the dimension sizes
                    file >> matrixSize[0] >> matrixSize[1] >> matrixSize[2];
                    std::getline(file, line);  // To consume the end of the line

                    xNodes.reserve(matrixSize[0]);
                    yNodes.reserve(matrixSize[1]);
                    zNodes.reserve(matrixSize[2]);
                    //std::cout << matrixSize[0] << "\n";
                }
                else if (line.find("#x node positions") != std::string::npos)
                {
                    // Read x node positions
                    std::getline(file, line);
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

    std::string GetHeaderFromFile(std::ifstream& file)
    {
        std::string line;
        std::string header = "";

        while (std::getline(file, line))
        {
            header += line + "\n";

            if (!(file.peek() == '#'))
                break;
        }

        return header;
    }

    int GetMaxIndex(std::filesystem::path energiesFile)
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

    std::array<float, 3> GetParametersFromDescriptionFileAtIndex(const std::filesystem::path& descriptionFile, int index)
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

    std::filesystem::path FindFileWithIndex(const std::filesystem::path& folder, int index)
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

    std::vector<std::string> SplitLine(std::string& string, const std::string& delimiter)
    {
        std::vector<std::string> tokens;

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

    std::vector<double> CalculateBinEdges(const std::vector<double>& binCenters, bool uniformDistances)
    {
        std::vector<double> binEdges;
        if (binCenters.empty())
        {
            return binEdges;
        }

        int nBins = binCenters.size();
        binEdges.reserve(nBins + 1);

        if (uniformDistances)
        {
            double firstEdge = binCenters[0] - (binCenters[1] - binCenters[0]) / 2.0;
            binEdges.push_back(firstEdge);

            // Compute middle edges as averages of adjacent bin centers
            for (int i = 0; i < nBins - 1; i++)
            {
                double edge = (binCenters[i] + binCenters[i + 1]) / 2.0;
                binEdges.push_back(edge);
            }

            // Compute last edge by extrapolation
            double lastEdge = binCenters[nBins - 1] + (binCenters[nBins - 1] - binCenters[nBins - 2]) / 2.0;
            binEdges.push_back(lastEdge);
        }
        else
        {
            // first edge is assumed to be 0
            binEdges.push_back(0);

            for (int i = 0; i < nBins; i++)
            {
                double edge = binEdges.back() + 2 * (binCenters[i] - binEdges.back());
                binEdges.push_back(edge);
            }
        }

        //for (const auto edge : binEdges)
        //{
        //    std::cout << "edge: " << edge << std::endl;
        //}

        return binEdges;
    }
}





