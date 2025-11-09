#include "../include/ResultManager.h"

ResultManager::ResultManager(const std::string& dir)
    : outputDirectory(dir),
    historicalBestEnergy(1e10),  // Initialize to large positive value
    historicalBestCluster(1, "A", "B"),
    hasHistoricalBest(false),
    currentRun(0) {

    fs::create_directories(dir);

    energyFile.open(dir + "/energy.txt");
    // Don't open bestEnergyFile here - will be opened per run

    historyFilePath = dir + "/historical_best.txt";
}

ResultManager::~ResultManager() {
    if (energyFile.is_open()) energyFile.close();
    if (bestEnergyFile.is_open()) bestEnergyFile.close();
}

void ResultManager::startNewRun(int runNumber) {
    currentRun = runNumber;

    // Close previous file if open
    if (bestEnergyFile.is_open()) {
        bestEnergyFile.close();
    }

    // Create filename with run number prefix
    std::string filename;
    if (runNumber > 0) {
        filename = outputDirectory + "/" + std::to_string(runNumber) + "_best_energy_per_generation.txt";
    }
    else {
        filename = outputDirectory + "/best_energy_per_generation.txt";
    }

    // Open in write mode (not append) to reset the file
    bestEnergyFile.open(filename, std::ios::out);

    if (bestEnergyFile.is_open()) {
        // Optionally write header
        bestEnergyFile << "# Generation\tEnergy\tComposition" << std::endl;
    }
}

void ResultManager::closeCurrentRun() {
    if (bestEnergyFile.is_open()) {
        bestEnergyFile.close();
    }
}

bool ResultManager::loadHistoricalBest(const std::string& compositionName) {
    std::ifstream file(historyFilePath);
    if (!file.is_open()) {
        std::cout << "No historical record found. Starting fresh." << std::endl;
        hasHistoricalBest = false;
        historicalBestEnergy = 1e10;
        return false;
    }

    // Read the file in Diamond format:
    // Line 1: number of atoms
    // Line 2: composition and energy
    // Lines 3+: atom coordinates

    std::string line;

    // Read number of atoms
    if (!std::getline(file, line)) {
        std::cout << "historical_best.txt is empty" << std::endl;
        file.close();
        hasHistoricalBest = false;
        historicalBestEnergy = 1e10;
        return false;
    }

    int numAtoms = std::stoi(line);

    // Read composition and energy
    if (!std::getline(file, line)) {
        std::cout << "historical_best.txt missing composition line" << std::endl;
        file.close();
        hasHistoricalBest = false;
        historicalBestEnergy = 1e10;
        return false;
    }

    std::istringstream iss(line);
    std::string comp;
    double energy;

    if (!(iss >> comp >> energy)) {
        std::cout << "Failed to parse composition and energy" << std::endl;
        file.close();
        hasHistoricalBest = false;
        historicalBestEnergy = 1e10;
        return false;
    }

    std::cout << "Found historical record: " << comp << " with energy " << energy << " eV" << std::endl;

    // Always load the historical best regardless of composition match
    historicalBestEnergy = energy;

    // Parse composition string (e.g., "Pt19Co19")
    std::string elemA = "", elemB = "";
    int numA = 0, numB = 0;
    size_t i = 0;

    // Get first element
    while (i < comp.length() && !std::isdigit(comp[i])) {
        elemA += comp[i++];
    }
    // Get first number
    while (i < comp.length() && std::isdigit(comp[i])) {
        numA = numA * 10 + (comp[i++] - '0');
    }
    // Get second element
    while (i < comp.length() && !std::isdigit(comp[i])) {
        elemB += comp[i++];
    }
    // Get second number
    while (i < comp.length() && std::isdigit(comp[i])) {
        numB = numB * 10 + (comp[i++] - '0');
    }

    historicalBestCluster = BinaryAlloyCluster(numA, numB, elemA, elemB);

    // Read coordinates
    for (int idx = 0; idx < numAtoms; ++idx) {
        if (!std::getline(file, line)) break;

        std::istringstream coordStream(line);
        std::string symbol;
        double x, y, z;
        if (coordStream >> symbol >> x >> y >> z) {
            historicalBestCluster.setAtomPosition(idx, x, y, z);
        }
    }

    historicalBestCluster.setEnergy(energy);
    hasHistoricalBest = true;

    file.close();

    if (comp != compositionName) {
        std::cout << "Warning: Historical best is for " << comp
            << " but current run is for " << compositionName << std::endl;
        std::cout << "Historical best energy (" << energy << " eV) will be used as reference." << std::endl;
    }
    else {
        std::cout << "Successfully loaded historical best for " << comp << std::endl;
    }

    return true;
}

void ResultManager::saveGenerationBest(int generation, const BinaryAlloyCluster& cluster) {
    if (!bestEnergyFile.is_open()) return;

    bestEnergyFile << generation << "\t"
        << std::setprecision(10) << cluster.getEnergy() << "\t"
        << cluster.getCompositionString() << std::endl;
    bestEnergyFile.flush();
}

bool ResultManager::updateHistoricalBest(const BinaryAlloyCluster& cluster) {
    double energy = cluster.getEnergy();

    // Only update if better than historical best
    if (!hasHistoricalBest || energy < historicalBestEnergy) {
        historicalBestEnergy = energy;
        historicalBestCluster = cluster;
        hasHistoricalBest = true;

        std::ofstream file(historyFilePath);
        if (file.is_open()) {
            // Write in Diamond format
            file << cluster.getNumAtoms() << std::endl;
            file << cluster.getCompositionString() << "  "
                << std::setprecision(12) << energy << std::endl;

            for (int i = 0; i < cluster.getNumAtoms(); ++i) {
                auto pos = cluster.getAtomPosition(i);
                file << cluster.getElementSymbol(i) << "\t"
                    << std::setprecision(10)
                    << pos[0] << "\t" << pos[1] << "\t" << pos[2] << std::endl;
            }

            file.close();

            std::cout << "\n*** New historical best found! ***" << std::endl;
            std::cout << "Energy: " << energy << " eV" << std::endl;

            // Save historical best structures
            saveXYZ(cluster, "historical_best.xyz");
            saveDiamond(cluster, "historical_best.txt");

            return true;
        }
    }

    return false;
}

void ResultManager::logGeneration(int gen, const std::vector<BinaryAlloyCluster>& population,
    const BinaryAlloyCluster& best) {
    if (!energyFile.is_open()) return;

    energyFile << gen << "\t" << best.getEnergy();
    for (const auto& ind : population) {
        energyFile << "\t" << ind.getEnergy();
    }
    energyFile << std::endl;
}

void ResultManager::saveXYZ(const BinaryAlloyCluster& cluster, const std::string& filename) {
    std::ofstream file(outputDirectory + "/" + filename);
    if (!file.is_open()) return;

    int n = cluster.getNumAtoms();
    file << n << std::endl;
    file << "Energy = " << std::setprecision(10) << cluster.getEnergy() << " eV, "
        << cluster.getCompositionString() << std::endl;

    for (int i = 0; i < n; ++i) {
        std::string symbol = cluster.getElementSymbol(i);
        auto pos = cluster.getAtomPosition(i);
        file << symbol << " " << std::setprecision(10)
            << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    }

    file.close();
}

void ResultManager::saveDiamond(const BinaryAlloyCluster& cluster, const std::string& filename) {
    std::ofstream file(outputDirectory + "/" + filename);
    if (!file.is_open()) return;

    int n = cluster.getNumAtoms();

    file << n << std::endl;
    file << cluster.getCompositionString() << "  "
        << std::setprecision(10) << cluster.getEnergy() << std::endl;

    const double* x = cluster.data();
    const double* y = x + n;
    const double* z = y + n;

    for (int i = 0; i < n; ++i) {
        std::string symbol = cluster.getElementSymbol(i);
        file << symbol << "\t" << std::setprecision(10)
            << x[i] << "\t" << y[i] << "\t" << z[i] << std::endl;
    }

    file.close();
}