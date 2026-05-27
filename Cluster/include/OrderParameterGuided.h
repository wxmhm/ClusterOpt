#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "PotentialBase.h"
#include "NELbfgs.h"
#include "ResultManager.h"
#include <vector>
#include <memory>
#include <map>

// Order Parameter Guided structure search with adaptive operator selection.
// Uses Steinhardt bond order parameters (Q4, Q6) to track structural diversity,
// Thompson sampling to learn which perturbation works best in each OP region,
// and UCB-style archive selection to balance exploration and exploitation.

class OrderParameterGuided {
public:
    enum class OperatorType { SINGLE = 0, COLLECTIVE = 1, SWAP = 2, COUNT = 3 };

    struct Parameters {
        // Execution
        int maxSteps = 5000;
        int archiveMaxSize = 100;
        bool useThreading = true;
        bool useLocalSearch = true;
        int localSearchFrequency = 1;

        // Metropolis acceptance
        double temperature = 0.1;
        double Tmin = 0.01;
        double Tmax = 2.0;

        // Perturbation
        double smallPerturbScale = 0.15;
        double largePerturbScale = 0.6;
        int collectiveGroupSize = 5;
        double swapProbability = 0.15;
        int stepsBeforeLargePerturb = 30;

        // OP-space
        double opDiversityWeight = 0.5;
        double opGridResolution = 0.03;

        // Thompson sampling for adaptive operator selection
        int operatorLearnWindow = 50;     // steps before Thompson sampling kicks in
        int opRegionBins = 5;             // bins per OP dimension for region discretization
        double thompsonStrength = 2.0;    // prior strength (higher = more exploration)

        // Adaptive temperature
        int stuckThreshold = 25;          // (unused with SA schedule)
        double heatFactor = 1.5;          // (unused with SA schedule)
        double coolFactor = 0.9;          // (unused with SA schedule)
    };

    struct OPSignature {
        double Q4 = 0.0;
        double Q6 = 0.0;
        double W4 = 0.0;
        double W6 = 0.0;

        double distance(const OPSignature& other) const {
            double dQ4 = Q4 - other.Q4;
            double dQ6 = Q6 - other.Q6;
            double dW4 = W4 - other.W4;
            double dW6 = W6 - other.W6;
            return std::sqrt(dQ4*dQ4 + dQ6*dQ6 + dW4*dW4 + dW6*dW6);
        }
    };

private:
    Parameters params;
    PotentialBase* potential;
    std::unique_ptr<ThreadPool> threadPool;
    std::vector<std::unique_ptr<NELbfgs>> lbfgsPool;

    struct ArchiveEntry {
        BinaryAlloyCluster structure;
        OPSignature signature;
        double energy;
        int visitCount;
    };
    std::vector<ArchiveEntry> archive;

    // Thompson sampling: per-OP-region operator success/failure counts
    struct OperatorStats {
        int successes[3] = {1, 1, 1};  // prior: Beta(1,1) uniform
        int failures[3] = {1, 1, 1};
    };
    std::map<int, OperatorStats> operatorStats;  // keyed by discretized OP region
    int totalOperatorAttempts;

    BinaryAlloyCluster globalBest;
    double globalBestEnergy;
    double currentTemperature;

    int numAtoms, numElementA, numElementB;
    std::string elementA, elementB;

    std::atomic<int> stepCount{0};
    std::atomic<int> acceptedCount{0};
    int stepsSinceNewBasin;
    int totalVisits;

public:
    OrderParameterGuided(const Parameters& p, PotentialBase* pot);
    BinaryAlloyCluster optimize(const BinaryAlloyCluster& initial,
        ResultManager* resultManager = nullptr);

    double getBestEnergy() const { return globalBestEnergy; }
    int getStepCount() const { return stepCount.load(); }

private:
    void initialize(const BinaryAlloyCluster& initial);
    double evaluateCluster(BinaryAlloyCluster& cluster, NELbfgs* lbfgs);

    OPSignature computeOPSignature(const BinaryAlloyCluster& cluster);

    // Archive management
    int findClosestArchiveEntry(const OPSignature& sig) const;
    void updateArchive(const BinaryAlloyCluster& cluster, const OPSignature& sig);
    int selectStructureFromArchive() const;  // UCB selection

    // Adaptive operator selection
    int discretizeOPRegion(const OPSignature& sig) const;
    OperatorType selectOperator(const OPSignature& sig);
    void updateOperatorStats(const OPSignature& sig, OperatorType op, bool success);

    // Perturbation operators
    BinaryAlloyCluster perturbSingle(BinaryAlloyCluster src, double scale);
    BinaryAlloyCluster perturbCollective(BinaryAlloyCluster src);
    BinaryAlloyCluster perturbSwap(const BinaryAlloyCluster& src);

    // Adaptive temperature
    void updateTemperature(int step);

    bool metropolisAccept(double oldE, double newE, double T);
};
