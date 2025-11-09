#pragma once
#include "Common.h"
#include "IDE.h"
#include "SaNSDE.h"

enum class AlgorithmType {
    IDE,
    SaNSDE,
    PSO
};

class Configuration {
public:
    struct SystemConfig {
        // General parameters
        int totalAtoms = 38;
        int numElementA = 38;
        int numElementB = 0;
        std::string elementA = "Pt";
        std::string elementB = "Co";
        std::string initialStructuresDir = "";

        // Algorithm selection
        AlgorithmType algorithm = AlgorithmType::IDE;

        // IDE-specific parameters
        IDE::Parameters ideParams;

        // SaNSDE-specific parameters
        SaNSDE::Parameters sansdeParams;

        // File and output parameters
        std::string potentialFile = "data/gupta_PtCo.txt";
        std::string outputDirectory = "results";
        bool saveIntermediates = true;
        int saveFrequency = 10;

        // Run control
        int numRuns = 1;
        bool runAllCompositions = false;  // New feature for running all combinations

        // Advanced options
        bool verbose = true;
        int randomSeed = -1;  // -1 means use time-based seed
        double convergenceTolerance = 1e-6;
        int stallGenerations = 50;  // Stop if no improvement for this many generations
    };

    static bool loadFromFile(const std::string& filename, SystemConfig& config);
    static bool saveToFile(const std::string& filename, const SystemConfig& config);
    static SystemConfig getDefaultConfig();
    static void printConfig(const SystemConfig& config);
    static std::vector<std::pair<int, int>> generateAllCompositions(int totalAtoms);

    static std::vector<BinaryAlloyCluster> loadInitialStructures(
        const std::string& directory,
        const std::string& elementA,
        const std::string& elementB
    );

    static bool loadStructureFromFile(
        const std::string& filename,
        BinaryAlloyCluster& cluster,
        const std::string& elementA,
        const std::string& elementB
    );

};
