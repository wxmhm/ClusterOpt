#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "PotentialBase.h"
#include "LocalOptimizer.h"
#include "ResultManager.h"
#include <thread>
#include <mutex>
#include <vector>
#include <memory>

enum class MutationStrategy {
    RAND1,
    RAND2,
    BEST1,
    RAND_TO_BEST1,
    CURRENT_TO_BEST1
};


class SaNSDE_Population {
public:
    struct Parameters {
        int populationSize = 30;
        int maxGenerations = 1000;
        int learningPeriod = 50;
        double F_min = 0.1, F_max = 1.0;
        double CR_min = 0.0, CR_max = 1.0;
        double p_min = 0.1;
        int neighborhoodSizeMin = 2;
        int neighborhoodSizeMax = 15;
        int memorySize = 100;
        bool useLocalSearch = true;
        int localSearchFrequency = 1;
    };

private:
    Parameters params;
    PotentialBase* potential;
    LocalOptimizer* localOptimizer;

    int populationId;
    std::string elementA;
    std::string elementB;
    int numElementA;
    int numElementB;

    std::vector<BinaryAlloyCluster> population;
    std::vector<BinaryAlloyCluster> mutantPopulation;
    std::vector<BinaryAlloyCluster> trialPopulation;
    BinaryAlloyCluster bestCluster;

    std::vector<double> F_values;
    std::vector<double> CR_values;
    std::vector<int> strategies;

    std::vector<MutationStrategy> availableStrategies;
    std::vector<double> strategyProbabilities;
    std::vector<int> strategySuccesses;
    std::vector<int> strategyFailures;

    struct SuccessMemory {
        double F, CR;
        int strategy;
        double improvement;
    };
    std::deque<SuccessMemory> successMemory;

    int generation;
    int evaluationCount;
    int localSearchCount;

public:
    SaNSDE_Population(int id, const Parameters& p, PotentialBase* pot,
        LocalOptimizer* opt, const BinaryAlloyCluster& initial);

    void evolve(int currentGen);
    void initializePopulation(const BinaryAlloyCluster& initial);
    double evaluateCluster(BinaryAlloyCluster& cluster);

    const BinaryAlloyCluster& getBestCluster() const { return bestCluster; }
    double getBestEnergy() const { return bestCluster.getEnergy(); }
    void receiveIndividual(const BinaryAlloyCluster& cluster);

private:
    void mutation(int index);
    void crossover(int index);
    void selection();
    void adaptParameters(int index);
    void updateStrategyProbabilities();
    int selectStrategy();
    double generateF(int index);
    double generateCR(int index);
    void neighborhoodSearch(int index, int currentGen);
    void sphereCutSplice(const BinaryAlloyCluster& parent1, const BinaryAlloyCluster& parent2,
        BinaryAlloyCluster& child1, BinaryAlloyCluster& child2);
};


class SaNSDE {
public:
    using Parameters = SaNSDE_Population::Parameters;

private:
    Parameters params;
    PotentialBase* potential;
    LocalOptimizer* localOptimizer;

    static const int NUM_POPULATIONS = 3;
    std::vector<std::unique_ptr<SaNSDE_Population>> populations;

    BinaryAlloyCluster globalBest;
    std::mutex globalBestMutex;

    int generation;
    bool useThreading;

public:
    SaNSDE(const Parameters& p, PotentialBase* pot, LocalOptimizer* opt,
        bool useThreading = false);

    void initialize(const BinaryAlloyCluster& initial);
    void evolve();
    BinaryAlloyCluster optimize(const BinaryAlloyCluster& initial,
        ResultManager* resultManager = nullptr);

    const BinaryAlloyCluster& getBestCluster() const { return globalBest; }
    int getGeneration() const { return generation; }
    int getEvaluationCount() const;
    int getLocalSearchCount() const;
    double getDiversity() const;

private:
    void updateGlobalBest(const BinaryAlloyCluster& candidate);
    void migrationBetweenPopulations();
};
