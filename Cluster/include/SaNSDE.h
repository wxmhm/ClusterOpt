#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "PotentialBase.h"
#include "ResultManager.h"
#include "NELbfgs.h"
#include <thread>
#include <mutex>
#include <vector>
#include <memory>
#include <deque>

enum class MutationStrategy {
    RAND1,
    RAND2,
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
        bool useThreading = false;
    };

private:
    Parameters params;
    PotentialBase* potential;

    ThreadPool* threadPool;
    std::vector<std::unique_ptr<NELbfgs>> lbfgsPool;

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
    std::atomic<int> evaluationCount{0};
    std::atomic<int> localSearchCount{0};

    // Parallel evaluation batch data
    struct SpliceInfo {
        int child2Idx;
        int popIndex;
    };
    std::vector<BinaryAlloyCluster> child2Pool;
    std::vector<double> child2Energies;
    std::vector<SpliceInfo> spliceInfos;
    std::vector<bool> usedSphereCutSplice;

public:
    SaNSDE_Population(int id, const Parameters& p, PotentialBase* pot,
        const BinaryAlloyCluster& initial, ThreadPool* pool);

    void evolve(int currentGen);
    void initializePopulation(const BinaryAlloyCluster& initial);
    double evaluateCluster(BinaryAlloyCluster& cluster, NELbfgs* lbfgs);

    const BinaryAlloyCluster& getBestCluster() const { return bestCluster; }
    double getBestEnergy() const { return bestCluster.getEnergy(); }
    int getEvaluationCount() const { return evaluationCount.load(); }
    int getLocalSearchCount() const { return localSearchCount.load(); }
    void receiveIndividual(const BinaryAlloyCluster& cluster);

private:
    // Phased parallel evolution
    void generateTrials(int currentGen);
    void evaluateTrials();
    void doNeighborhoodSearch(int currentGen);
    void applySelection();

    // Per-individual operations (used by generateTrials)
    void mutation(int index);
    void adaptParameters(int index);
    int selectStrategy();
    double generateF(int idx);
    double generateCR(int idx);

    // Strategy adaptation
    void updateStrategyProbabilities();

    // Crossover helpers
    void sphereCutSplice(const BinaryAlloyCluster& parent1, const BinaryAlloyCluster& parent2,
        BinaryAlloyCluster& child1, BinaryAlloyCluster& child2);
};


class SaNSDE {
public:
    using Parameters = SaNSDE_Population::Parameters;

private:
    Parameters params;
    PotentialBase* potential;

    static const int NUM_POPULATIONS = 3;
    std::vector<std::unique_ptr<SaNSDE_Population>> populations;

    std::unique_ptr<ThreadPool> threadPool;

    BinaryAlloyCluster globalBest;
    std::mutex globalBestMutex;

    int generation;
    bool useThreading;

public:
    SaNSDE(const Parameters& p, PotentialBase* pot, bool useThreading = false);

    void initialize(const BinaryAlloyCluster& initial);
    void evolve();
    BinaryAlloyCluster optimize(const BinaryAlloyCluster& initial,
        ResultManager* resultManager = nullptr);

    const BinaryAlloyCluster& getBestCluster() const { return globalBest; }
    int getGeneration() const { return generation; }

private:
    void updateGlobalBest(const BinaryAlloyCluster& candidate);
    void migrationBetweenPopulations();
};
