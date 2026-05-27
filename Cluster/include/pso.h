#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "PotentialBase.h"
#include "ResultManager.h"
#include "NELbfgs.h"
#include <thread>
#include <mutex>

struct PSO_Individual {
    BinaryAlloyCluster cluster;
    std::vector<double> velocity;
    BinaryAlloyCluster pbest;
    double energy;
    double pbestEnergy;

    // HCLPSO: per-atom exemplar indices + self-adaptive learning factor
    std::vector<int> exemplarIndices;
    int stagnationSinceRefresh;
    int improvementCount;    // consecutive improvements (for adaptive c)
    double adaptiveC;         // per-particle learning factor

    PSO_Individual(int nA, int nB, const std::string& elemA, const std::string& elemB)
        : cluster(nA, nB, elemA, elemB),
        pbest(nA, nB, elemA, elemB),
        energy(1e10),
        pbestEnergy(1e10),
        stagnationSinceRefresh(0),
        improvementCount(0),
        adaptiveC(1.49445) {
        int numAtoms = nA + nB;
        velocity.resize(3 * numAtoms, 0.0);
        exemplarIndices.resize(numAtoms, 0);
    }
};

class PSO_Population {
public:
    struct Parameters {
        int populationSize = 30;
        int maxGenerations = 200;
        double omegaMax = 0.9;
        double omegaMin = 0.4;
        double c1 = 2.5;
        double c2 = 1.5;
        double vmax = 0.1;
        double psoRatio = 0.6;
        bool useLocalSearch = true;
        int localSearchFrequency = 1;
        bool useThreading = false;

        // Stagnation / perturbation
        int stagnationThreshold = 10;
        int numPerturbations = 20;
        double perturbationScale = 0.05;
        int stagnationLimit = 15;

        // HCLPSO: Comprehensive Learning PSO (no gbest in velocity update)
        double c = 1.49445;          // base learning factor
        double cAdaptRate = 0.1;     // per-particle c adaptation rate
        int exemplarGap = 7;         // generations before refreshing exemplars
        bool useRingTopology = false; // ring neighborhood for exemplar selection
        double explorerRatio = 0.7;  // fraction of particles as CLPSO explorers
    };

private:
    Parameters params;
    PotentialBase* potential;

    ThreadPool* threadPool;
    std::vector<std::unique_ptr<NELbfgs>> lbfgsPool;

    std::vector<PSO_Individual> population;
    PSO_Individual globalBest;

    int numAtoms;
    int numElementA;
    int numElementB;
    std::string elementA;
    std::string elementB;
    double omega;
    int currentGeneration;

    std::atomic<int> evaluationCount{0};
    std::atomic<int> localSearchCount{0};

    // V2 stagnation tracking (retained for perturbGlobalBest)
    std::vector<int> stagnationCounters;
    int globalStagnationCounter;

    // CLPSO: particle fitness ranking for Pc computation
    std::vector<int> fitnessRank;

    // Parallel evaluation batch data
    struct EvalTarget {
        BinaryAlloyCluster* cluster;
        double* energy;
    };
    std::vector<EvalTarget> evalTargets;

public:
    PSO_Population(const Parameters& p, PotentialBase* pot,
        const BinaryAlloyCluster& initial, ThreadPool* pool);

    void evolve();
    void initializePopulation(const BinaryAlloyCluster& initial);
    double evaluateCluster(BinaryAlloyCluster& cluster, NELbfgs* lbfgs);

    const PSO_Individual& getGlobalBest() const { return globalBest; }
    int getEvaluationCount() const { return evaluationCount.load(); }
    int getLocalSearchCount() const { return localSearchCount.load(); }
    void receiveIndividual(const PSO_Individual& ind);

private:
    // Phased parallel evolution
    void generateCandidates();
    void evaluateCandidates();
    void applyResults();
    void perturbGlobalBest();

    // Per-particle operations
    void updateVelocityAndPosition(int idx);
    void generateRandomStructure(int idx);
    void swapAtoms();

    void initializeVelocities();
    void updateOmega();
    void clampVelocity(std::vector<double>& vel);

    // CLPSO methods
    void updateFitnessRanking();
    void selectExemplars(int particleIdx);
    void refreshExemplarsIfNeeded(int particleIdx);
};


class PSO {
public:
    using Parameters = PSO_Population::Parameters;

private:
    Parameters params;
    PotentialBase* potential;

    static const int NUM_POPULATIONS = 3;
    std::vector<std::unique_ptr<PSO_Population>> populations;

    std::unique_ptr<ThreadPool> threadPool;

    PSO_Individual globalBest;
    std::mutex globalBestMutex;

    int generation;
    int exchangeInterval;

public:
    PSO(const Parameters& p, PotentialBase* pot);

    void initialize(const BinaryAlloyCluster& initial);
    void evolve();
    BinaryAlloyCluster optimize(const BinaryAlloyCluster& initial,
        ResultManager* resultManager = nullptr);

    const PSO_Individual& getBestIndividual() const { return globalBest; }
    int getGeneration() const { return generation; }

private:
    void updateGlobalBest(const PSO_Individual& candidate);
    void exchangeBestIndividuals();
};
