#include "pch.h"

#include "FileHandler.h"

#include "tinyfiledialogs.h"

#include "RateCoefficient.h"
#include "CrossSection.h"


TH3D* FileHandler::LoadMatrixFile(const std::filesystem::path& filename)
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

EnergyDistribution FileHandler::LoadEnergyDistribution(std::filesystem::path& filename, bool loadSamples)
{
    // load the .asc file with the histogram data
    std::ifstream histFile(filename);

    // Check if the file was successfully opened
    if (!histFile.is_open())
    {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
        return EnergyDistribution();
    }

    std::string header = GetHeaderFromFile(histFile);
    EnergyDistribution energyDist = CreateEnergyDistFromHeader(header);

    std::string line;
    while (std::getline(histFile, line))
    {
        std::vector<std::string> tokens = SplitLine(line, xDelimiter);

        energyDist.binCenters.push_back(std::stod(tokens[0]));
        energyDist.binValues.push_back(std::stod(tokens[1]));
        energyDist.binValuesNormalised.push_back(std::stod(tokens[2]));
    }
    std::cout << "loaded file: " << filename.filename();

    // see if .samples file exist with collision energy data
    if (loadSamples)
    {
        std::filesystem::path samplesFilename = filename.replace_extension(".samples");
        if (std::filesystem::exists(samplesFilename))
        {
            std::ifstream samplesFile(samplesFilename);

            // Check if the file was successfully opened
            if (!samplesFile.is_open())
            {
                std::cerr << "Error: Could not open the file " << filename << std::endl;
                return energyDist;
            }

            header = GetHeaderFromFile(samplesFile);
            energyDist.collisionEnergies.reserve(energyDist.mcmcParameter.numberSamples);

            while (std::getline(samplesFile, line))
            {
                energyDist.collisionEnergies.push_back(std::stod(line));
            }
            std::cout << "\tsamples file found";
        }
    }
    std::cout << std::endl;
    
    return energyDist;
}

EnergyDistributionSet FileHandler::LoadEnergyDistributionSet(std::filesystem::path& folder)
{
    EnergyDistributionSet set = EnergyDistributionSet();

    if (!std::filesystem::exists(folder) || !std::filesystem::is_directory(folder)) 
    {
        std::cerr << "Invalid directory path!" << std::endl;
        return set;
    }

    // Iterate through the directory
    for (const auto& entry : std::filesystem::directory_iterator(folder)) 
    {
        if (entry.is_regular_file() && entry.path().extension() == ".asc")
        {
            std::filesystem::path file = entry.path();
            set.distributions.emplace_back(std::move(LoadEnergyDistribution(file, true)));
        }
    }
    set.folder = folder.parent_path().filename();
    set.subFolder = folder.filename();
    return set;
}

RateCoefficient FileHandler::LoadRateCoefficients(std::filesystem::path& filename)
{
    std::ifstream file(filename);
    RateCoefficient rc = RateCoefficient();

    // Check if the file was successfully opened
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
        return rc;
    }

    std::string line;
    // skip first line
    std::getline(file, line);
    while (std::getline(file, line))
    {
        std::vector<std::string> tokens = SplitLine(line, "\t");
        rc.detuningEnergies.push_back(std::stod(tokens[0]));
        rc.value.push_back(std::stod(tokens[1]));
        rc.error.push_back(std::stod(tokens[3]));
        rc.graph->AddPoint(std::stod(tokens[0]), std::stod(tokens[1]));
    }
    return rc;
}

CrossSection FileHandler::LoadCrossSection(std::filesystem::path& filename)
{
    std::ifstream file(filename);
    CrossSection cs = CrossSection();

    // Check if the file was successfully opened
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
        return cs;
    }

    std::string line;
    // skip first line
    std::getline(file, line);
    while (std::getline(file, line))
    {
        std::vector<std::string> tokens = SplitLine(line, "\t");
        cs.energies.push_back(std::stod(tokens[0]));
        cs.values.push_back(std::stod(tokens[1]));
        cs.errors.push_back(std::stod(tokens[2]));
    }
    return cs;
}

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

std::array<float, 3> FileHandler::GetParamtersFromDescriptionFileAtIndex(const std::filesystem::path& descriptionFile, int index)
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

std::filesystem::path FileHandler::SelectFile(const std::filesystem::path& startPath, const std::vector<const char*>& filterPatterns)
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

std::filesystem::path FileHandler::SelectFolder(const std::filesystem::path& startPath)
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

std::vector<std::filesystem::path> FileHandler::SelectFiles(const std::filesystem::path& startPath)
{
    const char* filterPatterns[] = { "*.asc" };
    const char* filePaths = tinyfd_openFileDialog(
        "Choose a file",               // Dialog title
        startPath.string().c_str(),    // Default path or file
        1,                             // Number of filters
        filterPatterns,                // Filter patterns (NULL for any file type)
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

std::filesystem::path FileHandler::FindFileWithIndex(const std::filesystem::path& folder, int index)
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

void FileHandler::SaveEnergyDistributionSetAsHist(EnergyDistributionSet& eDistSet)
{
    // set the output filepath
    std::filesystem::path folder = outputFolder.string() /        // general output folder
        std::filesystem::path("Energy Distribution Sets") /
        eDistSet.folder.filename() /                  // folder of corresponfding desription file
        eDistSet.subFolder.filename();              // subfolder with specific parameters

    // extract the directory 
    //std::filesystem::path dir = std::filesystem::path(file).parent_path();

    // Create the directories if they don't exist
    if (!std::filesystem::exists(folder)) 
    {
        std::filesystem::create_directories(folder);
    }
    for (const EnergyDistribution& edist : eDistSet.distributions)
    {
        std::filesystem::path file = folder / std::filesystem::path(edist.Filename() + ".asc");
        std::ofstream outfile(file);

        if (!outfile.is_open()) {
            std::cerr << "Error opening file" << std::endl;
            return;
        }
        outfile << edist.String();

        outfile << "# bin center [eV]\tbin value\tbin value normalised by bin width\n";
        for (int i = 0; i < edist.binCenters.size(); i++)
        {
            outfile << edist.binCenters[i] << "\t";
            outfile << edist.binValues[i] << "\t";
            outfile << edist.binValuesNormalised[i] << "\n";
        }

        outfile.close();
    }
}

void FileHandler::SaveEnergyDistributionSetAsSamples(EnergyDistributionSet& eDistSet)
{
    // set the output filepath
    std::filesystem::path folder = outputFolder.string() /                    // general output folder
        std::filesystem::path("Energy Distribution Sets") /
        eDistSet.folder.filename() /                             // folder of corresponfding desription file
        eDistSet.subFolder.filename();                          // subfolder with specific parameters

    // extract the directory 
    //std::filesystem::path dir = std::filesystem::path(folder).parent_path();

    // Create the directories if they don't exist
    if (!std::filesystem::exists(folder)) 
    {
        std::filesystem::create_directories(folder);
    }

    for (const EnergyDistribution& edist : eDistSet.distributions)
    {
        std::filesystem::path file = folder / std::filesystem::path(edist.Filename() + ".samples");
        std::ofstream outfile(file);

        if (!outfile.is_open())
        {
            std::cerr << "Error opening file" << std::endl;
            return;
        }

        outfile << edist.String();

        outfile << "# sampled collision energy values\n";
        for (double energy : edist.collisionEnergies)
        {
            outfile << energy << "\n";
        }

        outfile.close();
    }
}

void FileHandler::SaveRateCoefficients(RateCoefficient& rc)
{
    // set the output filepath
    std::filesystem::path file = outputFolder.string() /                    // general output folder
        std::filesystem::path("Rate Coefficients") /
        rc.file;
       
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

    for (int i = 0; i < rc.detuningEnergies.size(); i++)
    {
        outfile << rc.detuningEnergies[i] << "\t" << rc.value[i] << "\t" << 0.0 << "\t" << rc.error[i] << "\n";
    }

    outfile.close();
}

void FileHandler::SaveCrossSection(CrossSection& cs)
{
    // set the output filepath
    std::filesystem::path file = outputFolder.string() /                    // general output folder
        std::filesystem::path("Cross Sections") /
        cs.file;

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

    outfile << "# Energy [eV]\tCross Section Value\terror\n";

    for (int i = 0; i < cs.energies.size(); i++)
    {
        outfile << cs.energies[i] << "\t" << cs.values[i]  << "\t" << cs.errors[i] << "\n";
    }

    outfile.close();
}

std::string FileHandler::GetHeaderFromFile(std::ifstream& file) const
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

std::vector<std::string> FileHandler::SplitLine(std::string& string, const std::string& delimiter) const
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

TH3D* FileHandler::CreateTH3DfromHeader(std::ifstream& file) const
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

EnergyDistribution FileHandler::CreateEnergyDistFromHeader(std::string& header)
{
    EnergyDistribution energyDist;

    energyDist.eBeamParameter.fromString(header);
    energyDist.ionBeamParameter.fromString(header);
    energyDist.labEnergiesParameter.fromString(header);
    energyDist.mcmcParameter.fromString(header);
    energyDist.analyticalParameter.fromString(header);
    energyDist.simplifyParams.fromString(header);
              
    energyDist.SetupLabellingThings();

    return energyDist;
}
