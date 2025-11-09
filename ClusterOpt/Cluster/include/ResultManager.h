#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"

class ResultManager {
private:
    std::string outputDirectory;
    std::ofstream energyFile;
    std::ofstream bestEnergyFile;

    double historicalBestEnergy;
    BinaryAlloyCluster historicalBestCluster;
    bool hasHistoricalBest;
    std::string historyFilePath;

    int currentRun;  // Added to track current run number

public:
    ResultManager(const std::string& dir);
    ~ResultManager();

    void logGeneration(int gen, const std::vector<BinaryAlloyCluster>& population,
        const BinaryAlloyCluster& best);
    void saveXYZ(const BinaryAlloyCluster& cluster, const std::string& filename);
    void saveDiamond(const BinaryAlloyCluster& cluster, const std::string& filename);

    bool loadHistoricalBest(const std::string& compositionName);
    void saveGenerationBest(int generation, const BinaryAlloyCluster& cluster);
    bool updateHistoricalBest(const BinaryAlloyCluster& cluster);

    // New methods for handling runs
    void startNewRun(int runNumber);
    void closeCurrentRun();

    double getHistoricalBestEnergy() const { return historicalBestEnergy; }
    bool hasHistory() const { return hasHistoricalBest; }
    const BinaryAlloyCluster& getHistoricalBestCluster() const { return historicalBestCluster; }
};
