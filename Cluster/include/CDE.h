#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "PotentialBase.h"
#include "LocalOptimizer.h"
#include "ResultManager.h"
#include <thread>
#include <mutex>

enum CDE_MutationStrategy {
    CDE_RAND1,
    CDE_BEST1,
    CDE_RAND_TO_BEST1
};

struct CDE_Individual {
    BinaryAlloyCluster cluster;
    double energy;

    CDE_Individual(int nA, int nB, const std::string& elemA, const std::string& elemB)
        : cluster(nA, nB, elemA, elemB), energy(1e10) {}
};

class CDE_Population {
private:
    CDE_MutationStrategy strategy;
    std::vector<CDE_Individual> population;
    std::vector<CDE_Individual> mutantPopulation;
    std::vector<CDE_Individual> trialPopulation;
    CDE_Individual bestIndividual;

    PotentialBase* potential;
    LocalOptimizer* localOptimizer;
    int populationSize;
    int numAtoms;
    int numElementA;
    int numElementB;
    std::string elementA;
    std::string elementB;
    int localSearchCount;

public:
    CDE_Population(CDE_MutationStrategy strat, int popSize,
        const BinaryAlloyCluster& initial,
        PotentialBase* pot, LocalOptimizer* opt);
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

    void receiveIndividual(const CDE_Individual& ind);
    const CDE_Individual& getBestIndividual() const { return bestIndividual; }
};

class CDE {
public:
    struct Parameters {
        int populationSize = 30;
        int maxGenerations = 500;
        int exchangeInterval = 50;
        bool useLocalSearch = true;
        int localSearchFrequency = 1;
        bool useMultiPopulation = true;
        bool useThreading = false;
    };

private:
    Parameters params;
    std::vector<std::unique_ptr<CDE_Population>> populations;
    CDE_Individual globalBest;
    PotentialBase* potential;
    LocalOptimizer* localOptimizer;

    int generation;
    int evaluationCount;

    std::mutex globalBestMutex;

public:
    CDE(const Parameters& p, PotentialBase* pot, LocalOptimizer* opt);

    BinaryAlloyCluster optimize(const BinaryAlloyCluster& initial,
        ResultManager* resultManager = nullptr);

    void initialize(const BinaryAlloyCluster& initial);
    void evolve();
    void exchangeBestIndividuals();
    void updateGlobalBest(const CDE_Individual& candidate);

    const CDE_Individual& getBestIndividual() const { return globalBest; }
    int getGeneration() const { return generation; }
};
