#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "GuptaPotential.h"
#include "LocalOptimizer.h"
#include "ResultManager.h"
#include <thread>
#include <mutex>

struct PSO_Individual {
    BinaryAlloyCluster cluster;
    std::vector<double> velocity;
    BinaryAlloyCluster pbest;
    double energy;
    double pbestEnergy;

    PSO_Individual(int nA, int nB, const std::string& elemA, const std::string& elemB)
        : cluster(nA, nB, elemA, elemB),
        pbest(nA, nB, elemA, elemB),
        energy(1e10),
        pbestEnergy(1e10) {
        int numAtoms = nA + nB;
        velocity.resize(3 * numAtoms, 0.0);
    }
};

class PSO_Population {
private:
    std::vector<PSO_Individual> population;
    PSO_Individual globalBest;
    GuptaPotential* potential;
    LocalOptimizer* localOptimizer;
    int populationSize;
    int numAtoms;
    int numElementA;
    int numElementB;
    std::string elementA;
    std::string elementB;
    double omega;
    double omegaMax;
    double omegaMin;
    double c1;
    double c2;
    double vmax;
    double psoRatio;
    int currentGeneration;
    int maxGenerations;

public:
    PSO_Population(int popSize, const BinaryAlloyCluster& initial,
        GuptaPotential* pot, LocalOptimizer* opt, int maxGen);

    void evolve();
    void updateVelocityAndPosition(int idx);
    void generateRandomStructures(int startIdx, int count);
    double evaluateCluster(BinaryAlloyCluster& cluster);
    void updatePersonalAndGlobalBest();
    void swapAtoms();
    void receiveIndividual(const PSO_Individual& ind);
    const PSO_Individual& getGlobalBest() const { return globalBest; }

private:
    void initializeVelocities();
    void updateOmega();
    void clampVelocity(std::vector<double>& vel);
};

class PSO {
private:
    std::vector<std::unique_ptr<PSO_Population>> populations;
    PSO_Individual globalBest;
    GuptaPotential* potential;
    LocalOptimizer* localOptimizer;
    int generation;
    int populationSize;
    int maxGenerations;
    int exchangeInterval;
    bool useMultiPopulation;
    std::mutex globalBestMutex;

public:
    PSO(GuptaPotential* pot, LocalOptimizer* opt);
    BinaryAlloyCluster optimize(const BinaryAlloyCluster& initial, ResultManager* resultManager = nullptr);
    void initialize(const BinaryAlloyCluster& initial);
    void evolve();
    void exchangeBestIndividuals();
    void updateGlobalBest(const PSO_Individual& candidate);
    const PSO_Individual& getBestIndividual() const { return globalBest; }
    int getGeneration() const { return generation; }
};