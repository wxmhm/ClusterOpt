#include "../include/Configuration.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
Configuration::SystemConfig Configuration::getDefaultConfig() {
    SystemConfig config;

    // General parameters
    config.totalAtoms = 38;
    config.numElementA = 38;
    config.numElementB = 0;
    config.elementA = "Pt";
    config.elementB = "Co";

    // Algorithm selection
    config.algorithm = AlgorithmType::CDE;

    // CDE parameters
    config.cdeParams.populationSize = 40;
    config.cdeParams.maxGenerations = 200;
    config.cdeParams.exchangeInterval = 20;
    config.cdeParams.useLocalSearch = true;
    config.cdeParams.localSearchFrequency = 1;
    config.cdeParams.useMultiPopulation = true;

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
    config.sansdeParams.useThreading = false;

    // File and output parameters
    config.potentialFile = "data/gupta_PtCo.txt";
    config.outputDirectory = "results";
    config.numRuns = 1;
    config.runAllCompositions = false;

    // Advanced options
    config.verbose = true;

    // Potential type
    config.potentialType = PotentialType::Gupta;
    return config;
}

PotentialType Configuration::stringToPotentialType(const std::string& str) {
    if (str == "Gupta" || str == "gupta" || str == "GUPTA") {
        return PotentialType::Gupta;
    } else if (str == "FinnisSinclair" || str == "finnis-sinclair" || str == "FS" || str == "fs") {
        return PotentialType::FinnisSinclair;
    } else if (str == "SuttonChen" || str == "sutton-chen" || str == "SC" || str == "sc") {
        return PotentialType::SuttonChen;
    }
    return PotentialType::Gupta;
}

std::string Configuration::potentialTypeToString(PotentialType type) {
    switch (type) {
        case PotentialType::Gupta:
            return "Gupta";
        case PotentialType::FinnisSinclair:
            return "FinnisSinclair";
        case PotentialType::SuttonChen:
            return "SuttonChen";
        default:
            return "Gupta";
    }
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
                    if (value == "CDE" || value == "cde") config.algorithm = AlgorithmType::CDE;
                    else if (value == "SaNSDE" || value == "sansde") config.algorithm = AlgorithmType::SaNSDE;
                    else if (value == "PSO" || value == "pso") config.algorithm = AlgorithmType::PSO;
                }
                
                // Potential type selection
                else if (key == "potentialType") {
                    config.potentialType = stringToPotentialType(value);
                }

                // CDE parameters
                else if (key == "cde.populationSize") config.cdeParams.populationSize = std::stoi(value);
                else if (key == "cde.maxGenerations") config.cdeParams.maxGenerations = std::stoi(value);
                else if (key == "cde.exchangeInterval") config.cdeParams.exchangeInterval = std::stoi(value);
                else if (key == "cde.useLocalSearch") config.cdeParams.useLocalSearch = (value == "true" || value == "1");
                else if (key == "cde.localSearchFrequency") config.cdeParams.localSearchFrequency = std::stoi(value);
                else if (key == "cde.useMultiPopulation") config.cdeParams.useMultiPopulation = (value == "true" || value == "1");
                else if (key == "cde.useThreading") config.cdeParams.useThreading = (value == "true" || value == "1");

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
                else if (key == "sansde.useThreading") config.sansdeParams.useThreading = (value == "true" || value == "1");

                // File and output parameters
                else if (key == "potentialFile") config.potentialFile = value;
                else if (key == "outputDirectory") config.outputDirectory = value;
                else if (key == "numRuns") config.numRuns = std::stoi(value);
                else if (key == "runAllCompositions") config.runAllCompositions = (value == "true" || value == "1");

                // Advanced options
                else if (key == "verbose") config.verbose = (value == "true" || value == "1");
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

    file << "# Algorithm Selection (CDE, SaNSDE, or PSO)\n";
    file << "algorithm=" << (config.algorithm == AlgorithmType::CDE ? "CDE" : 
                             config.algorithm == AlgorithmType::SaNSDE ? "SaNSDE" : "PSO") << "\n\n";
    
    file << "# Potential Type Selection (Gupta, FinnisSinclair, SuttonChen)\n";
    file << "potentialType=" << potentialTypeToString(config.potentialType) << "\n\n";

    file << "# CDE Algorithm Parameters\n";
    file << "cde.populationSize=" << config.cdeParams.populationSize << "\n";
    file << "cde.maxGenerations=" << config.cdeParams.maxGenerations << "\n";
    file << "cde.exchangeInterval=" << config.cdeParams.exchangeInterval << "\n";
    file << "cde.useLocalSearch=" << (config.cdeParams.useLocalSearch ? "true" : "false") << "\n";
    file << "cde.localSearchFrequency=" << config.cdeParams.localSearchFrequency << "\n";
    file << "cde.useMultiPopulation=" << (config.cdeParams.useMultiPopulation ? "true" : "false") << "\n";
    file << "cde.useThreading=" << (config.cdeParams.useThreading ? "true" : "false") << "\n\n";

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
    file << "sansde.localSearchFrequency=" << config.sansdeParams.localSearchFrequency << "\n";
    file << "sansde.useThreading=" << (config.sansdeParams.useThreading ? "true" : "false") << "\n\n";

    file << "# File and Output Parameters\n";
    file << "potentialFile=" << config.potentialFile << "\n";
    file << "outputDirectory=" << config.outputDirectory << "\n\n";

    file << "# Run Control\n";
    file << "numRuns=" << config.numRuns << "\n";
    file << "runAllCompositions=" << (config.runAllCompositions ? "true" : "false") << "\n\n";

    file << "# Advanced Options\n";
    file << "verbose=" << (config.verbose ? "true" : "false") << "\n";

    file.close();
    return true;
}

void Configuration::printConfig(const SystemConfig& config) {
    std::cout << "\n=== Configuration Settings ===\n";
    std::cout << "Cluster: " << config.elementA << config.numElementA
        << config.elementB << config.numElementB << "\n";
    std::cout << "Algorithm: " << (config.algorithm == AlgorithmType::CDE ? "CDE" : 
                                    config.algorithm == AlgorithmType::SaNSDE ? "SaNSDE" : "PSO") << "\n";
    std::cout << "Potential Type:       " << potentialTypeToString(config.potentialType) << std::endl;

    if (config.algorithm == AlgorithmType::CDE) {
        std::cout << "CDE Parameters:\n";
        std::cout << "  Population Size: " << config.cdeParams.populationSize << "\n";
        std::cout << "  Max Generations: " << config.cdeParams.maxGenerations << "\n";
        std::cout << "  Multi-Population: " << (config.cdeParams.useMultiPopulation ? "Yes" : "No") << "\n";
    }
    else {
        std::cout << "SaNSDE Parameters:\n";
        std::cout << "  Population Size: " << config.sansdeParams.populationSize << "\n";
        std::cout << "  Max Generations: " << config.sansdeParams.maxGenerations << "\n";
        std::cout << "  F range: [" << config.sansdeParams.F_min << ", " << config.sansdeParams.F_max << "]\n";
        std::cout << "  CR range: [" << config.sansdeParams.CR_min << ", " << config.sansdeParams.CR_max << "]\n";
    }

    std::cout << "Local Search: " << ((config.algorithm == AlgorithmType::CDE ?
        config.cdeParams.useLocalSearch :
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

