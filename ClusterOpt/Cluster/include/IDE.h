#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "GuptaPotential.h"
#include "LocalOptimizer.h"
#include "ResultManager.h"
#include <thread>
#include <mutex>

enum IDE_MutationStrategy {
    IDE_RAND1,
    IDE_RAND2,
    IDE_BEST1,
    IDE_BEST2,
    IDE_RAND_TO_BEST1
};

struct IDE_Individual {
    BinaryAlloyCluster cluster;
    double energy;

    IDE_Individual(int nA, int nB, const std::string& elemA, const std::string& elemB)
        : cluster(nA, nB, elemA, elemB), energy(1e10) {}
};

class IDE_Population {
private:
    IDE_MutationStrategy strategy;
    std::vector<IDE_Individual> population;
    std::vector<IDE_Individual> mutantPopulation;
    std::vector<IDE_Individual> trialPopulation;
    IDE_Individual bestIndividual;

    GuptaPotential* potential;
    LocalOptimizer* localOptimizer;
    int populationSize;
    int numAtoms;
    int numElementA;
    int numElementB;
    std::string elementA;
    std::string elementB;
    int localSearchCount;

public:
    IDE_Population(IDE_MutationStrategy strat, int popSize,
        const BinaryAlloyCluster& initial,
        GuptaPotential* pot, LocalOptimizer* opt);
    void evolve();
    void mutation();
    void crossover();
    void selection();
    void swapAtoms();
    double evaluateCluster(BinaryAlloyCluster& cluster);
    void sphereCutSplice(const BinaryAlloyCluster& parent1,
        const BinaryAlloyCluster& parent2,
        BinaryAlloyCluster& child1,
        BinaryAlloyCluster& child2);

    void receiveIndividual(const IDE_Individual& ind);
    const IDE_Individual& getBestIndividual() const { return bestIndividual; }
};

class IDE {
public:
    struct Parameters {
        int populationSize = 30;
        int maxGenerations = 500;
        int exchangeInterval = 50;
        bool useLocalSearch = true;
        int localSearchFrequency = 1;
        bool useMultiPopulation = true;
        bool useThreading = false;  // 默认启用多线程
    };

private:
    Parameters params;
    std::vector<std::unique_ptr<IDE_Population>> populations;
    IDE_Individual globalBest;
    GuptaPotential* potential;
    LocalOptimizer* localOptimizer;

    int generation;
    int evaluationCount;

    std::mutex globalBestMutex;

public:
    IDE(const Parameters& p, GuptaPotential* pot, LocalOptimizer* opt);

    BinaryAlloyCluster optimize(const BinaryAlloyCluster& initial,
        ResultManager* resultManager = nullptr);

    void initialize(const BinaryAlloyCluster& initial);
    void evolve();
    void exchangeBestIndividuals();
    void updateGlobalBest(const IDE_Individual& candidate);

    const IDE_Individual& getBestIndividual() const { return globalBest; }
    int getGeneration() const { return generation; }
};
