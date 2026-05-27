#include "../include/Configuration.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
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

    // PSO parameters
    config.psoParams.populationSize = 30;
    config.psoParams.maxGenerations = 200;
    config.psoParams.omegaMax = 0.9;
    config.psoParams.omegaMin = 0.4;
    config.psoParams.c1 = 2.5;
    config.psoParams.c2 = 1.5;
    config.psoParams.vmax = 0.1;
    config.psoParams.psoRatio = 0.6;
    config.psoParams.useLocalSearch = true;
    config.psoParams.localSearchFrequency = 1;
    config.psoParams.useThreading = false;
    config.psoParams.stagnationThreshold = 10;
    config.psoParams.numPerturbations = 20;
    config.psoParams.perturbationScale = 0.05;
    config.psoParams.stagnationLimit = 15;
    config.psoParams.c = 1.49445;
    config.psoParams.cAdaptRate = 0.1;
    config.psoParams.exemplarGap = 7;
    config.psoParams.useRingTopology = false;
    config.psoParams.explorerRatio = 0.7;

    // BH parameters
    config.bhParams.T0 = 1.0;
    config.bhParams.Tfinal = 0.01;
    config.bhParams.coolingRate = 0.97;
    config.bhParams.stepsPerT = 50;
    config.bhParams.probSingle = 0.3;
    config.bhParams.probCollective = 0.4;
    config.bhParams.probSwap = 0.3;
    config.bhParams.singleStepSize = 0.3;
    config.bhParams.collectiveStepSize = 0.5;
    config.bhParams.collectiveSize = 4;
    config.bhParams.numReplicas = 5;
    config.bhParams.exchangeInterval = 10;
    config.bhParams.maxTotalSteps = 5000;
    config.bhParams.useThreading = true;
    config.bhParams.useLocalSearch = true;
    config.bhParams.localSearchFrequency = 1;
    config.bhParams.useSurfaceBias = true;
    config.bhParams.surfaceBiasStrength = 0.7;
    config.bhParams.useStructuralMemory = true;
    config.bhParams.opMemoryResolution = 0.03;
    config.bhParams.opMemoryMaxSize = 200;
    config.bhParams.escapeSteps = 5;
    config.bhParams.useAdaptiveStep = true;
    config.bhParams.targetAcceptRate = 0.5;
    config.bhParams.stepAdaptRate = 0.05;

    // OPG parameters
    config.opgParams.maxSteps = 5000;
    config.opgParams.archiveMaxSize = 100;
    config.opgParams.useThreading = true;
    config.opgParams.useLocalSearch = true;
    config.opgParams.localSearchFrequency = 1;
    config.opgParams.temperature = 0.1;
    config.opgParams.Tmin = 0.01;
    config.opgParams.Tmax = 2.0;
    config.opgParams.smallPerturbScale = 0.15;
    config.opgParams.largePerturbScale = 0.6;
    config.opgParams.collectiveGroupSize = 5;
    config.opgParams.swapProbability = 0.15;
    config.opgParams.stepsBeforeLargePerturb = 30;
    config.opgParams.opDiversityWeight = 0.5;
    config.opgParams.opGridResolution = 0.03;
    config.opgParams.operatorLearnWindow = 50;
    config.opgParams.opRegionBins = 5;
    config.opgParams.thompsonStrength = 2.0;
    config.opgParams.stuckThreshold = 25;
    config.opgParams.heatFactor = 1.5;
    config.opgParams.coolFactor = 0.9;

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

                try {
                    // General parameters
                    if (key == "totalAtoms") config.totalAtoms = std::stoi(value);
                    else if (key == "numElementA") config.numElementA = std::stoi(value);
                    else if (key == "elementA") config.elementA = value;
                    else if (key == "elementB") config.elementB = value;

                    // Algorithm selection
                    else if (key == "algorithm") {
                        if (value == "CDE" || value == "cde") config.algorithm = AlgorithmType::CDE;
                        else if (value == "SaNSDE" || value == "sansde") config.algorithm = AlgorithmType::SaNSDE;
                        else if (value == "PSO" || value == "pso") config.algorithm = AlgorithmType::PSO;
                        else if (value == "BH" || value == "bh") config.algorithm = AlgorithmType::BH;
                        else if (value == "OPG" || value == "opg") config.algorithm = AlgorithmType::OPG;
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
                    else if (key == "cde.numThreads") config.cdeParams.numThreads = std::stoi(value);

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

                    // PSO parameters
                    else if (key == "pso.populationSize") config.psoParams.populationSize = std::stoi(value);
                    else if (key == "pso.maxGenerations") config.psoParams.maxGenerations = std::stoi(value);
                    else if (key == "pso.omegaMax") config.psoParams.omegaMax = std::stod(value);
                    else if (key == "pso.omegaMin") config.psoParams.omegaMin = std::stod(value);
                    else if (key == "pso.c1") config.psoParams.c1 = std::stod(value);
                    else if (key == "pso.c2") config.psoParams.c2 = std::stod(value);
                    else if (key == "pso.vmax") config.psoParams.vmax = std::stod(value);
                    else if (key == "pso.psoRatio") config.psoParams.psoRatio = std::stod(value);
                    else if (key == "pso.useLocalSearch") config.psoParams.useLocalSearch = (value == "true" || value == "1");
                    else if (key == "pso.localSearchFrequency") config.psoParams.localSearchFrequency = std::stoi(value);
                    else if (key == "pso.useThreading") config.psoParams.useThreading = (value == "true" || value == "1");
                    else if (key == "pso.stagnationThreshold") config.psoParams.stagnationThreshold = std::stoi(value);
                    else if (key == "pso.numPerturbations") config.psoParams.numPerturbations = std::stoi(value);
                    else if (key == "pso.perturbationScale") config.psoParams.perturbationScale = std::stod(value);
                    else if (key == "pso.stagnationLimit") config.psoParams.stagnationLimit = std::stoi(value);
                    else if (key == "pso.c") config.psoParams.c = std::stod(value);
                    else if (key == "pso.cAdaptRate") config.psoParams.cAdaptRate = std::stod(value);
                    else if (key == "pso.exemplarGap") config.psoParams.exemplarGap = std::stoi(value);
                    else if (key == "pso.useRingTopology") config.psoParams.useRingTopology = (value == "true" || value == "1");
                    else if (key == "pso.explorerRatio") config.psoParams.explorerRatio = std::stod(value);

                    // BH parameters
                    else if (key == "bh.T0") config.bhParams.T0 = std::stod(value);
                    else if (key == "bh.Tfinal") config.bhParams.Tfinal = std::stod(value);
                    else if (key == "bh.coolingRate") config.bhParams.coolingRate = std::stod(value);
                    else if (key == "bh.stepsPerT") config.bhParams.stepsPerT = std::stoi(value);
                    else if (key == "bh.probSingle") config.bhParams.probSingle = std::stod(value);
                    else if (key == "bh.probCollective") config.bhParams.probCollective = std::stod(value);
                    else if (key == "bh.probSwap") config.bhParams.probSwap = std::stod(value);
                    else if (key == "bh.singleStepSize") config.bhParams.singleStepSize = std::stod(value);
                    else if (key == "bh.collectiveStepSize") config.bhParams.collectiveStepSize = std::stod(value);
                    else if (key == "bh.collectiveSize") config.bhParams.collectiveSize = std::stoi(value);
                    else if (key == "bh.numReplicas") config.bhParams.numReplicas = std::stoi(value);
                    else if (key == "bh.exchangeInterval") config.bhParams.exchangeInterval = std::stoi(value);
                    else if (key == "bh.maxTotalSteps") config.bhParams.maxTotalSteps = std::stoi(value);
                    else if (key == "bh.useThreading") config.bhParams.useThreading = (value == "true" || value == "1");
                    else if (key == "bh.useLocalSearch") config.bhParams.useLocalSearch = (value == "true" || value == "1");
                    else if (key == "bh.localSearchFrequency") config.bhParams.localSearchFrequency = std::stoi(value);
                    else if (key == "bh.useSurfaceBias") config.bhParams.useSurfaceBias = (value == "true" || value == "1");
                    else if (key == "bh.surfaceBiasStrength") config.bhParams.surfaceBiasStrength = std::stod(value);
                    else if (key == "bh.useStructuralMemory") config.bhParams.useStructuralMemory = (value == "true" || value == "1");
                    else if (key == "bh.opMemoryResolution") config.bhParams.opMemoryResolution = std::stod(value);
                    else if (key == "bh.opMemoryMaxSize") config.bhParams.opMemoryMaxSize = std::stoi(value);
                    else if (key == "bh.escapeSteps") config.bhParams.escapeSteps = std::stoi(value);
                    else if (key == "bh.useAdaptiveStep") config.bhParams.useAdaptiveStep = (value == "true" || value == "1");
                    else if (key == "bh.targetAcceptRate") config.bhParams.targetAcceptRate = std::stod(value);
                    else if (key == "bh.stepAdaptRate") config.bhParams.stepAdaptRate = std::stod(value);

                    // OPG parameters
                    else if (key == "opg.maxSteps") config.opgParams.maxSteps = std::stoi(value);
                    else if (key == "opg.archiveMaxSize") config.opgParams.archiveMaxSize = std::stoi(value);
                    else if (key == "opg.useThreading") config.opgParams.useThreading = (value == "true" || value == "1");
                    else if (key == "opg.useLocalSearch") config.opgParams.useLocalSearch = (value == "true" || value == "1");
                    else if (key == "opg.localSearchFrequency") config.opgParams.localSearchFrequency = std::stoi(value);
                    else if (key == "opg.temperature") config.opgParams.temperature = std::stod(value);
                    else if (key == "opg.smallPerturbScale") config.opgParams.smallPerturbScale = std::stod(value);
                    else if (key == "opg.largePerturbScale") config.opgParams.largePerturbScale = std::stod(value);
                    else if (key == "opg.collectiveGroupSize") config.opgParams.collectiveGroupSize = std::stoi(value);
                    else if (key == "opg.swapProbability") config.opgParams.swapProbability = std::stod(value);
                    else if (key == "opg.stepsBeforeLargePerturb") config.opgParams.stepsBeforeLargePerturb = std::stoi(value);
                    else if (key == "opg.opDiversityWeight") config.opgParams.opDiversityWeight = std::stod(value);
                    else if (key == "opg.opGridResolution") config.opgParams.opGridResolution = std::stod(value);
                    else if (key == "opg.Tmin") config.opgParams.Tmin = std::stod(value);
                    else if (key == "opg.Tmax") config.opgParams.Tmax = std::stod(value);
                    else if (key == "opg.operatorLearnWindow") config.opgParams.operatorLearnWindow = std::stoi(value);
                    else if (key == "opg.opRegionBins") config.opgParams.opRegionBins = std::stoi(value);
                    else if (key == "opg.thompsonStrength") config.opgParams.thompsonStrength = std::stod(value);
                    else if (key == "opg.stuckThreshold") config.opgParams.stuckThreshold = std::stoi(value);
                    else if (key == "opg.heatFactor") config.opgParams.heatFactor = std::stod(value);
                    else if (key == "opg.coolFactor") config.opgParams.coolFactor = std::stod(value);

                    // File and output parameters
                    else if (key == "potentialFile") config.potentialFile = value;
                    else if (key == "outputDirectory") config.outputDirectory = value;
                    else if (key == "numRuns") config.numRuns = std::stoi(value);
                    else if (key == "runAllCompositions") config.runAllCompositions = (value == "true" || value == "1");

                    // Advanced options
                    else if (key == "verbose") config.verbose = (value == "true" || value == "1");
                } catch (const std::exception& e) {
                    std::cerr << "Warning: Failed to parse config key '" << key << "': " << e.what() << std::endl;
                }
            }
        }
    }

    // numElementB is always derived from totalAtoms - numElementA
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

    file << "# Algorithm Selection (CDE, SaNSDE, PSO, or BH)\n";
    std::string algStr;
    switch (config.algorithm) {
        case AlgorithmType::CDE: algStr = "CDE"; break;
        case AlgorithmType::SaNSDE: algStr = "SaNSDE"; break;
        case AlgorithmType::PSO: algStr = "PSO"; break;
        case AlgorithmType::BH: algStr = "BH"; break;
        case AlgorithmType::OPG: algStr = "OPG"; break;
    }
    file << "algorithm=" << algStr << "\n\n";
    
    file << "# Potential Type Selection (Gupta, FinnisSinclair, SuttonChen)\n";
    file << "potentialType=" << potentialTypeToString(config.potentialType) << "\n\n";

    file << "# CDE Algorithm Parameters\n";
    file << "cde.populationSize=" << config.cdeParams.populationSize << "\n";
    file << "cde.maxGenerations=" << config.cdeParams.maxGenerations << "\n";
    file << "cde.exchangeInterval=" << config.cdeParams.exchangeInterval << "\n";
    file << "cde.useLocalSearch=" << (config.cdeParams.useLocalSearch ? "true" : "false") << "\n";
    file << "cde.localSearchFrequency=" << config.cdeParams.localSearchFrequency << "\n";
    file << "cde.useMultiPopulation=" << (config.cdeParams.useMultiPopulation ? "true" : "false") << "\n";
    file << "cde.useThreading=" << (config.cdeParams.useThreading ? "true" : "false") << "\n";
    file << "cde.numThreads=" << config.cdeParams.numThreads << "\n\n";

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

    file << "# PSO Algorithm Parameters\n";
    file << "pso.populationSize=" << config.psoParams.populationSize << "\n";
    file << "pso.maxGenerations=" << config.psoParams.maxGenerations << "\n";
    file << "pso.omegaMax=" << config.psoParams.omegaMax << "\n";
    file << "pso.omegaMin=" << config.psoParams.omegaMin << "\n";
    file << "pso.c1=" << config.psoParams.c1 << "\n";
    file << "pso.c2=" << config.psoParams.c2 << "\n";
    file << "pso.vmax=" << config.psoParams.vmax << "\n";
    file << "pso.psoRatio=" << config.psoParams.psoRatio << "\n";
    file << "pso.useLocalSearch=" << (config.psoParams.useLocalSearch ? "true" : "false") << "\n";
    file << "pso.localSearchFrequency=" << config.psoParams.localSearchFrequency << "\n";
    file << "pso.useThreading=" << (config.psoParams.useThreading ? "true" : "false") << "\n";
    file << "pso.stagnationThreshold=" << config.psoParams.stagnationThreshold << "\n";
    file << "pso.numPerturbations=" << config.psoParams.numPerturbations << "\n";
    file << "pso.perturbationScale=" << config.psoParams.perturbationScale << "\n";
    file << "pso.stagnationLimit=" << config.psoParams.stagnationLimit << "\n";
    file << "pso.c=" << config.psoParams.c << "\n";
    file << "pso.cAdaptRate=" << config.psoParams.cAdaptRate << "\n";
    file << "pso.exemplarGap=" << config.psoParams.exemplarGap << "\n";
    file << "pso.useRingTopology=" << (config.psoParams.useRingTopology ? "true" : "false") << "\n";
    file << "pso.explorerRatio=" << config.psoParams.explorerRatio << "\n\n";

    file << "# BH Algorithm Parameters\n";
    file << "bh.T0=" << config.bhParams.T0 << "\n";
    file << "bh.Tfinal=" << config.bhParams.Tfinal << "\n";
    file << "bh.coolingRate=" << config.bhParams.coolingRate << "\n";
    file << "bh.stepsPerT=" << config.bhParams.stepsPerT << "\n";
    file << "bh.probSingle=" << config.bhParams.probSingle << "\n";
    file << "bh.probCollective=" << config.bhParams.probCollective << "\n";
    file << "bh.probSwap=" << config.bhParams.probSwap << "\n";
    file << "bh.singleStepSize=" << config.bhParams.singleStepSize << "\n";
    file << "bh.collectiveStepSize=" << config.bhParams.collectiveStepSize << "\n";
    file << "bh.collectiveSize=" << config.bhParams.collectiveSize << "\n";
    file << "bh.numReplicas=" << config.bhParams.numReplicas << "\n";
    file << "bh.exchangeInterval=" << config.bhParams.exchangeInterval << "\n";
    file << "bh.maxTotalSteps=" << config.bhParams.maxTotalSteps << "\n";
    file << "bh.useThreading=" << (config.bhParams.useThreading ? "true" : "false") << "\n";
    file << "bh.useLocalSearch=" << (config.bhParams.useLocalSearch ? "true" : "false") << "\n";
    file << "bh.localSearchFrequency=" << config.bhParams.localSearchFrequency << "\n\n";

    file << "# OPG Algorithm Parameters\n";
    file << "opg.maxSteps=" << config.opgParams.maxSteps << "\n";
    file << "opg.archiveMaxSize=" << config.opgParams.archiveMaxSize << "\n";
    file << "opg.useThreading=" << (config.opgParams.useThreading ? "true" : "false") << "\n";
    file << "opg.useLocalSearch=" << (config.opgParams.useLocalSearch ? "true" : "false") << "\n";
    file << "opg.localSearchFrequency=" << config.opgParams.localSearchFrequency << "\n";
    file << "opg.temperature=" << config.opgParams.temperature << "\n";
    file << "opg.smallPerturbScale=" << config.opgParams.smallPerturbScale << "\n";
    file << "opg.largePerturbScale=" << config.opgParams.largePerturbScale << "\n";
    file << "opg.collectiveGroupSize=" << config.opgParams.collectiveGroupSize << "\n";
    file << "opg.swapProbability=" << config.opgParams.swapProbability << "\n";
    file << "opg.stepsBeforeLargePerturb=" << config.opgParams.stepsBeforeLargePerturb << "\n";
    file << "opg.opDiversityWeight=" << config.opgParams.opDiversityWeight << "\n";
    file << "opg.opGridResolution=" << config.opgParams.opGridResolution << "\n\n";

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
    std::string algName;
    switch (config.algorithm) {
        case AlgorithmType::CDE: algName = "CDE"; break;
        case AlgorithmType::SaNSDE: algName = "SaNSDE"; break;
        case AlgorithmType::PSO: algName = "PSO"; break;
        case AlgorithmType::BH: algName = "BH"; break;
        case AlgorithmType::OPG: algName = "OPG"; break;
    }
    std::cout << "Algorithm: " << algName << "\n";
    std::cout << "Potential Type:       " << potentialTypeToString(config.potentialType) << std::endl;

    if (config.algorithm == AlgorithmType::CDE) {
        std::cout << "CDE Parameters:\n";
        std::cout << "  Population Size: " << config.cdeParams.populationSize << "\n";
        std::cout << "  Max Generations: " << config.cdeParams.maxGenerations << "\n";
        std::cout << "  Multi-Population: " << (config.cdeParams.useMultiPopulation ? "Yes" : "No") << "\n";
        std::cout << "  Threading: " << (config.cdeParams.useThreading ? "Yes" : "No");
        if (config.cdeParams.useThreading) {
            std::cout << " (" << config.cdeParams.numThreads << " threads)";
        }
        std::cout << "\n";
        std::cout << "Local Search: " << (config.cdeParams.useLocalSearch ? "Enabled" : "Disabled") << "\n";
    }
    else if (config.algorithm == AlgorithmType::SaNSDE) {
        std::cout << "SaNSDE Parameters:\n";
        std::cout << "  Population Size: " << config.sansdeParams.populationSize << "\n";
        std::cout << "  Max Generations: " << config.sansdeParams.maxGenerations << "\n";
        std::cout << "  F range: [" << config.sansdeParams.F_min << ", " << config.sansdeParams.F_max << "]\n";
        std::cout << "  CR range: [" << config.sansdeParams.CR_min << ", " << config.sansdeParams.CR_max << "]\n";
        std::cout << "Local Search: " << (config.sansdeParams.useLocalSearch ? "Enabled" : "Disabled") << "\n";
    }
    else if (config.algorithm == AlgorithmType::PSO) {
        std::cout << "PSO Parameters:\n";
        std::cout << "  Population Size: " << config.psoParams.populationSize << "\n";
        std::cout << "  Max Generations: " << config.psoParams.maxGenerations << "\n";
        std::cout << "  Omega: [" << config.psoParams.omegaMax << " -> " << config.psoParams.omegaMin << "]\n";
        std::cout << "  c1: " << config.psoParams.c1 << ", c2: " << config.psoParams.c2 << "\n";
        std::cout << "  Vmax: " << config.psoParams.vmax << ", PSO Ratio: " << config.psoParams.psoRatio << "\n";
        std::cout << "  Threading: " << (config.psoParams.useThreading ? "Yes" : "No") << "\n";
        std::cout << "  Local Search: " << (config.psoParams.useLocalSearch ? "Enabled" : "Disabled") << "\n";
        std::cout << "  Stagnation Threshold: " << config.psoParams.stagnationThreshold << "\n";
        std::cout << "  Perturbations: " << config.psoParams.numPerturbations
            << " (scale=" << config.psoParams.perturbationScale << ")\n";
        std::cout << "  Velocity Reset Limit: " << config.psoParams.stagnationLimit << "\n";
        std::cout << "  HCLPSO: c=" << config.psoParams.c
            << ", adaptRate=" << config.psoParams.cAdaptRate
            << ", exemplarGap=" << config.psoParams.exemplarGap
            << ", explorerRatio=" << config.psoParams.explorerRatio << "\n";
        std::cout << "  Ring Topology: " << (config.psoParams.useRingTopology ? "Yes" : "No") << "\n";
    }
    else if (config.algorithm == AlgorithmType::BH) {
        std::cout << "BH Parameters:\n";
        std::cout << "  T range: [" << config.bhParams.T0 << " -> " << config.bhParams.Tfinal << "]\n";
        std::cout << "  Cooling Rate: " << config.bhParams.coolingRate
            << ", Steps per T: " << config.bhParams.stepsPerT << "\n";
        std::cout << "  Max Total Steps: " << config.bhParams.maxTotalSteps << "\n";
        std::cout << "  Move Probs: single=" << config.bhParams.probSingle
            << ", collective=" << config.bhParams.probCollective
            << ", swap=" << config.bhParams.probSwap << "\n";
        std::cout << "  Step Sizes: single=" << config.bhParams.singleStepSize
            << ", collective=" << config.bhParams.collectiveStepSize
            << ", collectiveSize=" << config.bhParams.collectiveSize << "\n";
        std::cout << "  Replicas: " << config.bhParams.numReplicas
            << ", exchangeInterval=" << config.bhParams.exchangeInterval << "\n";
        std::cout << "  Threading: " << (config.bhParams.useThreading ? "Yes" : "No") << "\n";
        std::cout << "  Local Search: " << (config.bhParams.useLocalSearch ? "Enabled" : "Disabled")
            << " (freq=" << config.bhParams.localSearchFrequency << ")\n";
    }
    else if (config.algorithm == AlgorithmType::OPG) {
        std::cout << "OPG Parameters:\n";
        std::cout << "  Max Steps: " << config.opgParams.maxSteps << "\n";
        std::cout << "  Archive Max Size: " << config.opgParams.archiveMaxSize
            << ", Resolution: " << config.opgParams.opGridResolution << "\n";
        std::cout << "  Temperature: " << config.opgParams.temperature
            << ", Diversity Weight: " << config.opgParams.opDiversityWeight << "\n";
        std::cout << "  Perturb Scales: small=" << config.opgParams.smallPerturbScale
            << ", large=" << config.opgParams.largePerturbScale << "\n";
        std::cout << "  Collective Group Size: " << config.opgParams.collectiveGroupSize
            << ", Swap Prob: " << config.opgParams.swapProbability << "\n";
        std::cout << "  Steps Before Large Perturb: " << config.opgParams.stepsBeforeLargePerturb << "\n";
        std::cout << "  Threading: " << (config.opgParams.useThreading ? "Yes" : "No") << "\n";
        std::cout << "  Local Search: " << (config.opgParams.useLocalSearch ? "Enabled" : "Disabled")
            << " (freq=" << config.opgParams.localSearchFrequency << ")\n";
    }
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

