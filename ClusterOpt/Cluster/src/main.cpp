#include "../include/Common.h"
#include "../include/BinaryAlloyCluster.h"
#include "../include/GuptaPotential.h"
#include "../include/LocalOptimizer.h"
#include "../include/IDE.h"
#include "../include/SaNSDE.h"
#include "../include/PSO.h"
#include "../include/ResultManager.h"
#include "../include/Configuration.h"
#include <memory>

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -config <file>         Configuration file (default: data/config.txt)" << std::endl;
    std::cout << "  -atoms <n>             Total number of atoms" << std::endl;
    std::cout << "  -composition <nA>      Number of element A atoms" << std::endl;
    std::cout << "  -elements <A> <B>      Element symbols" << std::endl;
    std::cout << "  -algorithm <IDE|SaNSDE|PSO> Select optimization algorithm" << std::endl;
    std::cout << "  -generations <n>       Maximum generations" << std::endl;
    std::cout << "  -runs <n>              Number of independent runs" << std::endl;
    std::cout << "  -allcompositions       Run all possible compositions for given atom count" << std::endl;
    std::cout << "  -saveconfig <file>     Save current configuration to file" << std::endl;
    std::cout << "  -help                  Show this help" << std::endl;
}

BinaryAlloyCluster runSingleOptimization(
    const Configuration::SystemConfig& config,
    GuptaPotential* potential,
    LocalOptimizer* localOptimizer,
    ResultManager* resultManager,
    const BinaryAlloyCluster& initial) {

    BinaryAlloyCluster best = initial;

    if (config.algorithm == AlgorithmType::IDE) {
        IDE optimizer(config.ideParams, potential, localOptimizer);
        best = optimizer.optimize(initial, resultManager);
    }
    else if (config.algorithm == AlgorithmType::SaNSDE) {
        SaNSDE optimizer(config.sansdeParams, potential, localOptimizer);
        best = optimizer.optimize(initial, resultManager);
    }
    else if (config.algorithm == AlgorithmType::PSO) {
        PSO optimizer(potential, localOptimizer);
        best = optimizer.optimize(initial, resultManager);
    }
    else {
        std::cerr << "Unknown algorithm type!" << std::endl;
    }

    return best;
}

void runOptimizationForComposition(
    Configuration::SystemConfig config,
    int nA, int nB) {

    config.numElementA = nA;
    config.numElementB = nB;

    std::cout << "\n========================================" << std::endl;
    std::cout << "Optimizing: " << config.elementA << nA << config.elementB << nB << std::endl;
    std::cout << "Algorithm: " << (config.algorithm == AlgorithmType::IDE ? "IDE" :
        config.algorithm == AlgorithmType::SaNSDE ? "SaNSDE" : "PSO") << std::endl;
    std::cout << "========================================" << std::endl;

    GuptaPotential potential(config.elementA, config.elementB);
    if (!potential.loadParameters(config.potentialFile)) {
        std::cerr << "Warning: Could not load potential parameters from "
            << config.potentialFile << ", using defaults" << std::endl;
    }

    LocalOptimizer localOptimizer(&potential);

    std::string compositionName = config.elementA + std::to_string(nA) +
        config.elementB + std::to_string(nB);
    std::string outputDir = config.outputDirectory + "/" + compositionName;

    ResultManager resultManager(outputDir);
    resultManager.loadHistoricalBest(compositionName);

    if (resultManager.hasHistory()) {
        std::cout << "Historical best energy: "
            << resultManager.getHistoricalBestEnergy() << " eV" << std::endl;
    }

    std::vector<double> bestEnergies;
    BinaryAlloyCluster currentBest(nA, nB, config.elementA, config.elementB);
    double currentBestEnergy = 1e10;

    for (int run = 0; run < config.numRuns; ++run) {
        if (config.verbose) {
            std::cout << "\nRun " << (run + 1) << "/" << config.numRuns << std::endl;
        }

        resultManager.startNewRun(run + 1);

        BinaryAlloyCluster initial(nA, nB, config.elementA, config.elementB);
        initial.randomInitialize(2.75);

        BinaryAlloyCluster best = runSingleOptimization(
            config, &potential, &localOptimizer, &resultManager, initial);

        bestEnergies.push_back(best.getEnergy());

        if (config.verbose) {
            std::cout << "Energy: " << best.getEnergy() << " eV" << std::endl;
        }

        if (best.getEnergy() < currentBestEnergy) {
            currentBestEnergy = best.getEnergy();
            currentBest = best;

            if (config.verbose) {
                std::cout << "New best energy: " << currentBestEnergy << " eV" << std::endl;
            }

            resultManager.saveXYZ(currentBest, "best.xyz");
            resultManager.saveDiamond(currentBest, "best.txt");
        }

        resultManager.closeCurrentRun();
    }

    if (config.numRuns > 1) {
        double avgEnergy = std::accumulate(bestEnergies.begin(), bestEnergies.end(), 0.0)
            / config.numRuns;
        double minEnergy = *std::min_element(bestEnergies.begin(), bestEnergies.end());
        double maxEnergy = *std::max_element(bestEnergies.begin(), bestEnergies.end());

        double variance = 0.0;
        for (double e : bestEnergies) {
            variance += (e - avgEnergy) * (e - avgEnergy);
        }
        variance /= config.numRuns;
        double stdDev = std::sqrt(variance);

        std::cout << "\n=== Statistics for " << compositionName << " ===" << std::endl;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Best:     " << minEnergy << " eV" << std::endl;
        std::cout << "Worst:    " << maxEnergy << " eV" << std::endl;
        std::cout << "Average:  " << avgEnergy << " eV" << std::endl;
        std::cout << "Std Dev:  " << stdDev << " eV" << std::endl;

        double threshold = minEnergy * 1.01;
        int successCount = std::count_if(bestEnergies.begin(), bestEnergies.end(),
            [threshold](double e) { return e <= threshold; });
        double successRate = 100.0 * successCount / config.numRuns;
        std::cout << "Success Rate (within 1% of best): " << successRate << "%" << std::endl;
    }

    if (resultManager.hasHistory() && currentBestEnergy < resultManager.getHistoricalBestEnergy()) {
        std::cout << "\n*** New historical best found! ***" << std::endl;
        std::cout << "Previous: " << resultManager.getHistoricalBestEnergy() << " eV" << std::endl;
        std::cout << "New:      " << currentBestEnergy << " eV" << std::endl;
    }
}

int main(int argc, char* argv[]) {
    ElementInfo::initialize();

    Configuration::SystemConfig config = Configuration::getDefaultConfig();
    std::string configFile = "data/config.txt";
    bool saveConfigRequested = false;
    std::string saveConfigFile = "";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-help" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        }
        else if (arg == "-config" && i + 1 < argc) {
            configFile = argv[++i];
        }
        else if (arg == "-atoms" && i + 1 < argc) {
            config.totalAtoms = std::atoi(argv[++i]);
            config.numElementB = config.totalAtoms - config.numElementA;
        }
        else if (arg == "-composition" && i + 1 < argc) {
            config.numElementA = std::atoi(argv[++i]);
            config.numElementB = config.totalAtoms - config.numElementA;
        }
        else if (arg == "-elements" && i + 2 < argc) {
            config.elementA = argv[++i];
            config.elementB = argv[++i];
        }
        else if (arg == "-algorithm" && i + 1 < argc) {
            std::string alg = argv[++i];
            if (alg == "IDE" || alg == "ide") {
                config.algorithm = AlgorithmType::IDE;
            }
            else if (alg == "SaNSDE" || alg == "sansde") {
                config.algorithm = AlgorithmType::SaNSDE;
            }
            else if (alg == "PSO" || alg == "pso") {
                config.algorithm = AlgorithmType::PSO;
            }
        }
        else if (arg == "-generations" && i + 1 < argc) {
            int gens = std::atoi(argv[++i]);
            config.ideParams.maxGenerations = gens;
            config.sansdeParams.maxGenerations = gens;
        }
        else if (arg == "-runs" && i + 1 < argc) {
            config.numRuns = std::atoi(argv[++i]);
        }
        else if (arg == "-allcompositions") {
            config.runAllCompositions = true;
        }
        else if (arg == "-saveconfig" && i + 1 < argc) {
            saveConfigRequested = true;
            saveConfigFile = argv[++i];
        }
    }

    if (fs::exists(configFile)) {
        if (!Configuration::loadFromFile(configFile, config)) {
            std::cerr << "Error loading configuration file, using defaults" << std::endl;
        }
        else {
            std::cout << "Configuration loaded from: " << configFile << std::endl;
        }
    }
    else {
        std::cout << "Configuration file not found, using defaults" << std::endl;
    }

    if (saveConfigRequested) {
        if (Configuration::saveToFile(saveConfigFile, config)) {
            std::cout << "Configuration saved to: " << saveConfigFile << std::endl;
        }
        else {
            std::cerr << "Failed to save configuration to: " << saveConfigFile << std::endl;
        }
    }

    Configuration::printConfig(config);

    fs::create_directories(config.outputDirectory);

    if (config.runAllCompositions) {
        std::cout << "\n=== Running all compositions for "
            << config.totalAtoms << " atoms ===" << std::endl;

        auto compositions = Configuration::generateAllCompositions(config.totalAtoms);

        std::string summaryFile = config.outputDirectory + "/all_compositions_summary.txt";
        std::ofstream summary(summaryFile);
        summary << "# Composition\tBest_Energy\tAvg_Energy\tStd_Dev\tRuns\n";

        double globalBestEnergy = 1e10;
        std::string globalBestComposition;

        for (const auto& comp : compositions) {
            int nA = comp.first;
            int nB = comp.second;
            runOptimizationForComposition(config, nA, nB);
        }

        summary.close();
        std::cout << "\n=== All compositions completed ===" << std::endl;
        std::cout << "Summary saved to: " << summaryFile << std::endl;

    }
    else {
        runOptimizationForComposition(config, config.numElementA, config.numElementB);
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "Optimization completed successfully!" << std::endl;
    std::cout << "Results saved to: " << config.outputDirectory << std::endl;
    std::cout << "========================================" << std::endl;

    std::cout << "\nPress Enter to exit..." << std::endl;
    std::cin.get();

    return 0;
}
