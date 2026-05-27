#pragma once
#include "Common.h"
#include "CDE.h"
#include "SaNSDE.h"
#include "pso.h"
#include "BasinHopping.h"
#include "OrderParameterGuided.h"

enum class AlgorithmType {
    CDE,
    SaNSDE,
    PSO,
    BH,
    OPG
};

enum class PotentialType {
    Gupta,
    FinnisSinclair,
    SuttonChen
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

        // Algorithm selection
        AlgorithmType algorithm = AlgorithmType::CDE;

        // Potential type selection
        PotentialType potentialType = PotentialType::Gupta;

        // CDE-specific parameters
        CDE::Parameters cdeParams;

        // SaNSDE-specific parameters
        SaNSDE::Parameters sansdeParams;

        // PSO-specific parameters
        PSO::Parameters psoParams;

        // BH-specific parameters
        BasinHopping::Parameters bhParams;

        // OPG-specific parameters
        OrderParameterGuided::Parameters opgParams;

        // File and output parameters
        std::string potentialFile = "data/gupta_PtCo.txt";
        std::string outputDirectory = "results";

        // Run control
        int numRuns = 1;
        bool runAllCompositions = false;  // New feature for running all combinations

        // Advanced options
        bool verbose = true;
    };

    static bool loadFromFile(const std::string& filename, SystemConfig& config);
    static bool saveToFile(const std::string& filename, const SystemConfig& config);
    static SystemConfig getDefaultConfig();
    static void printConfig(const SystemConfig& config);
    static std::vector<std::pair<int, int>> generateAllCompositions(int totalAtoms);

    // Helper function to convert string to PotentialType
    static PotentialType stringToPotentialType(const std::string& str);
    static std::string potentialTypeToString(PotentialType type);
};
