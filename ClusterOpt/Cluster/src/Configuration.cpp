#include "../include/Configuration.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <regex>

Configuration::SystemConfig Configuration::getDefaultConfig() {
    SystemConfig config;

    // General parameters
    config.totalAtoms = 38;
    config.numElementA = 38;
    config.numElementB = 0;
    config.elementA = "Pt";
    config.elementB = "Co";

    // Algorithm selection
    config.algorithm = AlgorithmType::IDE;

    // IDE parameters
    config.ideParams.populationSize = 40;
    config.ideParams.maxGenerations = 200;
    config.ideParams.exchangeInterval = 20;
    config.ideParams.useLocalSearch = true;
    config.ideParams.localSearchFrequency = 1;
    config.ideParams.useMultiPopulation = true;

    // SaNSDE parameters
    config.sansdeParams.populationSize = 40;
    config.sansdeParams.maxGenerations = 200;
    config.sansdeParams.learningPeriod = 20;
    config.sansdeParams.F_min = 0.1;
    config.sansdeParams.F_max = 1.0;
    config.sansdeParams.CR_min = 0.0;
    config.sansdeParams.CR_max = 1.0;
    config.sansdeParams.p_min = 0.1;
    config.sansdeParams.neighborhoodSizeMin = 2;
    config.sansdeParams.neighborhoodSizeMax = 15;
    config.sansdeParams.memorySize = 100;
    config.sansdeParams.useLocalSearch = true;
    config.sansdeParams.localSearchFrequency = 1;

    // File and output parameters
    config.potentialFile = "data/gupta_PtCo.txt";
    config.outputDirectory = "results";
    config.saveIntermediates = true;
    config.saveFrequency = 10;
    config.numRuns = 1;
    config.runAllCompositions = false;

    // Advanced options
    config.verbose = true;
    config.randomSeed = -1;
    config.convergenceTolerance = 1e-6;
    config.stallGenerations = 50;

    return config;
}

bool Configuration::loadFromFile(const std::string& filename, SystemConfig& config) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Cannot open configuration file: " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key;
        if (std::getline(iss, key, '=')) {
            std::string value;
            if (std::getline(iss, value)) {
                // Trim whitespace
                key.erase(0, key.find_first_not_of(" \t"));
                key.erase(key.find_last_not_of(" \t") + 1);
                value.erase(0, value.find_first_not_of(" \t"));
                value.erase(value.find_last_not_of(" \t") + 1);

                // General parameters
                if (key == "totalAtoms") config.totalAtoms = std::stoi(value);
                else if (key == "numElementA") config.numElementA = std::stoi(value);
                else if (key == "numElementB") config.numElementB = std::stoi(value);
                else if (key == "elementA") config.elementA = value;
                else if (key == "elementB") config.elementB = value;

                // Algorithm selection
                else if (key == "algorithm") {
                    if (value == "IDE" || value == "ide") config.algorithm = AlgorithmType::IDE;
                    else if (value == "SaNSDE" || value == "sansde") config.algorithm = AlgorithmType::SaNSDE;
                    else if (value == "PSO" || value == "pso") config.algorithm = AlgorithmType::PSO;
                }

                // IDE parameters
                else if (key == "ide.populationSize") config.ideParams.populationSize = std::stoi(value);
                else if (key == "ide.maxGenerations") config.ideParams.maxGenerations = std::stoi(value);
                else if (key == "ide.exchangeInterval") config.ideParams.exchangeInterval = std::stoi(value);
                else if (key == "ide.useLocalSearch") config.ideParams.useLocalSearch = (value == "true" || value == "1");
                else if (key == "ide.localSearchFrequency") config.ideParams.localSearchFrequency = std::stoi(value);
                else if (key == "ide.useMultiPopulation") config.ideParams.useMultiPopulation = (value == "true" || value == "1");

                // SaNSDE parameters
                else if (key == "sansde.populationSize") config.sansdeParams.populationSize = std::stoi(value);
                else if (key == "sansde.maxGenerations") config.sansdeParams.maxGenerations = std::stoi(value);
                else if (key == "sansde.learningPeriod") config.sansdeParams.learningPeriod = std::stoi(value);
                else if (key == "sansde.F_min") config.sansdeParams.F_min = std::stod(value);
                else if (key == "sansde.F_max") config.sansdeParams.F_max = std::stod(value);
                else if (key == "sansde.CR_min") config.sansdeParams.CR_min = std::stod(value);
                else if (key == "sansde.CR_max") config.sansdeParams.CR_max = std::stod(value);
                else if (key == "sansde.p_min") config.sansdeParams.p_min = std::stod(value);
                else if (key == "sansde.neighborhoodSizeMin") config.sansdeParams.neighborhoodSizeMin = std::stoi(value);
                else if (key == "sansde.neighborhoodSizeMax") config.sansdeParams.neighborhoodSizeMax = std::stoi(value);
                else if (key == "sansde.memorySize") config.sansdeParams.memorySize = std::stoi(value);
                else if (key == "sansde.useLocalSearch") config.sansdeParams.useLocalSearch = (value == "true" || value == "1");
                else if (key == "sansde.localSearchFrequency") config.sansdeParams.localSearchFrequency = std::stoi(value);

                // File and output parameters
                else if (key == "potentialFile") config.potentialFile = value;
                else if (key == "outputDirectory") config.outputDirectory = value;
                else if (key == "saveIntermediates") config.saveIntermediates = (value == "true" || value == "1");
                else if (key == "saveFrequency") config.saveFrequency = std::stoi(value);
                else if (key == "numRuns") config.numRuns = std::stoi(value);
                else if (key == "runAllCompositions") config.runAllCompositions = (value == "true" || value == "1");

                // Advanced options
                else if (key == "verbose") config.verbose = (value == "true" || value == "1");
                else if (key == "randomSeed") config.randomSeed = std::stoi(value);
                else if (key == "convergenceTolerance") config.convergenceTolerance = std::stod(value);
                else if (key == "stallGenerations") config.stallGenerations = std::stoi(value);
            }
        }
    }

    // Auto-calculate numElementB if not specified
    config.numElementB = config.totalAtoms - config.numElementA;

    file.close();
    return true;
}

bool Configuration::saveToFile(const std::string& filename, const SystemConfig& config) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Cannot create configuration file: " << filename << std::endl;
        return false;
    }

    file << "# Binary Alloy Cluster Optimization Configuration File\n";
    file << "# ====================================================\n\n";

    file << "# General Parameters\n";
    file << "totalAtoms=" << config.totalAtoms << "\n";
    file << "numElementA=" << config.numElementA << "\n";
    file << "numElementB=" << config.numElementB << "\n";
    file << "elementA=" << config.elementA << "\n";
    file << "elementB=" << config.elementB << "\n\n";

    file << "# Algorithm Selection (IDE or SaNSDE)\n";
    file << "algorithm=" << (config.algorithm == AlgorithmType::IDE ? "IDE" : "SaNSDE") << "\n\n";

    file << "# IDE Algorithm Parameters\n";
    file << "ide.populationSize=" << config.ideParams.populationSize << "\n";
    file << "ide.maxGenerations=" << config.ideParams.maxGenerations << "\n";
    file << "ide.exchangeInterval=" << config.ideParams.exchangeInterval << "\n";
    file << "ide.useLocalSearch=" << (config.ideParams.useLocalSearch ? "true" : "false") << "\n";
    file << "ide.localSearchFrequency=" << config.ideParams.localSearchFrequency << "\n";
    file << "ide.useMultiPopulation=" << (config.ideParams.useMultiPopulation ? "true" : "false") << "\n\n";

    file << "# SaNSDE Algorithm Parameters\n";
    file << "sansde.populationSize=" << config.sansdeParams.populationSize << "\n";
    file << "sansde.maxGenerations=" << config.sansdeParams.maxGenerations << "\n";
    file << "sansde.learningPeriod=" << config.sansdeParams.learningPeriod << "\n";
    file << "sansde.F_min=" << config.sansdeParams.F_min << "\n";
    file << "sansde.F_max=" << config.sansdeParams.F_max << "\n";
    file << "sansde.CR_min=" << config.sansdeParams.CR_min << "\n";
    file << "sansde.CR_max=" << config.sansdeParams.CR_max << "\n";
    file << "sansde.p_min=" << config.sansdeParams.p_min << "\n";
    file << "sansde.neighborhoodSizeMin=" << config.sansdeParams.neighborhoodSizeMin << "\n";
    file << "sansde.neighborhoodSizeMax=" << config.sansdeParams.neighborhoodSizeMax << "\n";
    file << "sansde.memorySize=" << config.sansdeParams.memorySize << "\n";
    file << "sansde.useLocalSearch=" << (config.sansdeParams.useLocalSearch ? "true" : "false") << "\n";
    file << "sansde.localSearchFrequency=" << config.sansdeParams.localSearchFrequency << "\n\n";

    file << "# File and Output Parameters\n";
    file << "potentialFile=" << config.potentialFile << "\n";
    file << "outputDirectory=" << config.outputDirectory << "\n";
    file << "saveIntermediates=" << (config.saveIntermediates ? "true" : "false") << "\n";
    file << "saveFrequency=" << config.saveFrequency << "\n\n";

    file << "# Run Control\n";
    file << "numRuns=" << config.numRuns << "\n";
    file << "runAllCompositions=" << (config.runAllCompositions ? "true" : "false") << "\n\n";

    file << "# Advanced Options\n";
    file << "verbose=" << (config.verbose ? "true" : "false") << "\n";
    file << "randomSeed=" << config.randomSeed << "\n";
    file << "convergenceTolerance=" << config.convergenceTolerance << "\n";
    file << "stallGenerations=" << config.stallGenerations << "\n";

    file.close();
    return true;
}

void Configuration::printConfig(const SystemConfig& config) {
    std::cout << "\n=== Configuration Settings ===\n";
    std::cout << "Cluster: " << config.elementA << config.numElementA
        << config.elementB << config.numElementB << "\n";
    std::cout << "Algorithm: " << (config.algorithm == AlgorithmType::IDE ? "IDE" : "SaNSDE") << "\n";

    if (config.algorithm == AlgorithmType::IDE) {
        std::cout << "IDE Parameters:\n";
        std::cout << "  Population Size: " << config.ideParams.populationSize << "\n";
        std::cout << "  Max Generations: " << config.ideParams.maxGenerations << "\n";
        std::cout << "  Multi-Population: " << (config.ideParams.useMultiPopulation ? "Yes" : "No") << "\n";
    }
    else {
        std::cout << "SaNSDE Parameters:\n";
        std::cout << "  Population Size: " << config.sansdeParams.populationSize << "\n";
        std::cout << "  Max Generations: " << config.sansdeParams.maxGenerations << "\n";
        std::cout << "  F range: [" << config.sansdeParams.F_min << ", " << config.sansdeParams.F_max << "]\n";
        std::cout << "  CR range: [" << config.sansdeParams.CR_min << ", " << config.sansdeParams.CR_max << "]\n";
    }

    std::cout << "Local Search: " << ((config.algorithm == AlgorithmType::IDE ?
        config.ideParams.useLocalSearch :
        config.sansdeParams.useLocalSearch) ? "Enabled" : "Disabled") << "\n";
    std::cout << "Number of Runs: " << config.numRuns << "\n";
    std::cout << "Run All Compositions: " << (config.runAllCompositions ? "Yes" : "No") << "\n";
    std::cout << "==============================\n\n";
}

std::vector<std::pair<int, int>> Configuration::generateAllCompositions(int totalAtoms) {
    std::vector<std::pair<int, int>> compositions;

    // For binary alloy, generate all combinations from (0, totalAtoms) to (totalAtoms, 0)
    // Skip pure elements if desired
    for (int nA = 0; nA <= totalAtoms; ++nA) {
        int nB = totalAtoms - nA;
        compositions.push_back(std::make_pair(nA, nB));
    }

    return compositions;
}


bool Configuration::loadStructureFromFile(
    const std::string& filename,
    BinaryAlloyCluster& cluster,
    const std::string& elementA,
    const std::string& elementB
) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    std::vector<std::array<double, 3>> positions;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string symbol;
        double x, y, z;

        if (iss >> symbol >> x >> y >> z) {
            positions.push_back({ x, y, z });
        }
    }

    if (positions.size() != static_cast<size_t>(cluster.getNumAtoms())) {
        std::cerr << "Warning: Number of atoms in file (" << positions.size()
            << ") doesn't match cluster size (" << cluster.getNumAtoms() << ")" << std::endl;
        return false;
    }

    for (size_t i = 0; i < positions.size(); ++i) {
        cluster.setAtomPosition(i, positions[i][0], positions[i][1], positions[i][2]);
    }

    std::cout << "Loaded structure from " << filename
        << " with " << positions.size() << " atoms" << std::endl;

    return true;
}

std::vector<BinaryAlloyCluster> Configuration::loadInitialStructures(
    const std::string& directory,
    const std::string& elementA,
    const std::string& elementB
) {
    std::vector<BinaryAlloyCluster> structures;

    if (directory.empty() || !fs::exists(directory)) {
        std::cout << "Initial structures directory not specified or doesn't exist" << std::endl;
        return structures;
    }

    std::cout << "Loading initial structures from: " << directory << std::endl;

    std::regex pattern("tmp[0-9]+\\.txt");

    std::vector<fs::path> matchedFiles;
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();
            if (std::regex_match(filename, pattern)) {
                matchedFiles.push_back(entry.path());
            }
        }
    }

    std::sort(matchedFiles.begin(), matchedFiles.end());

    std::cout << "Found " << matchedFiles.size() << " tmp*.txt files" << std::endl;

    for (const auto& filepath : matchedFiles) {
        std::ifstream checkFile(filepath);
        int atomCount = 0;
        std::string line;
        while (std::getline(checkFile, line)) {
            if (!line.empty() && line[0] != '#') {
                atomCount++;
            }
        }
        checkFile.close();

        if (atomCount == 0) continue;

        int countA = atomCount / 2;
        int countB = atomCount - countA;

        BinaryAlloyCluster cluster(countA, countB, elementA, elementB);

        if (loadStructureFromFile(filepath.string(), cluster, elementA, elementB)) {
            structures.push_back(cluster);
        }
    }

    std::cout << "Successfully loaded " << structures.size() << " initial structures" << std::endl;

    return structures;
}

