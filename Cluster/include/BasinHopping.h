#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "PotentialBase.h"
#include "NELbfgs.h"
#include "ResultManager.h"
#include <vector>
#include <memory>

// Basin Hopping with surface-biased moves and structural memory
// BH is the gold standard for cluster global optimization (Wales & Doye, 1997)
// Extended with: coordination-weighted atom selection, OP-space structural memory,
// adaptive step size, collective displacements, type swaps, parallel tempering

class BasinHopping {
public:
    struct Parameters {
        // Temperature schedule
        double T0 = 1.0;              // initial temperature (eV)
        double Tfinal = 0.01;         // final temperature
        double coolingRate = 0.97;    // multiplicative cooling per step (slower)
        int stepsPerT = 50;           // MC steps per temperature

        // Move probabilities
        double probSingle = 0.3;      // single atom displacement
        double probCollective = 0.4;  // collective group displacement
        double probSwap = 0.3;        // atom type swap

        // Move parameters
        double singleStepSize = 0.3;  // single atom displacement sigma
        double collectiveStepSize = 0.5; // collective displacement sigma
        int collectiveSize = 4;       // atoms in collective move

        // Parallel tempering
        int numReplicas = 5;          // replicas at different T (increased)
        int exchangeInterval = 10;    // attempt replica exchange every N steps

        // Execution
        int maxTotalSteps = 5000;     // total MC steps
        bool useThreading = true;
        bool useLocalSearch = true;
        int localSearchFrequency = 1;

        // Surface bias: atoms with fewer neighbors are more likely to be moved
        bool useSurfaceBias = true;
        double surfaceBiasStrength = 0.7;  // 0=uniform, 1=strong surface preference

        // Structural memory: detect and escape revisited OP-space regions
        bool useStructuralMemory = true;
        double opMemoryResolution = 0.03;   // OP-space resolution for memory
        int opMemoryMaxSize = 200;           // max entries in OP memory
        int escapeSteps = 5;                // collective moves when revisiting

        // Adaptive step size
        bool useAdaptiveStep = true;
        double targetAcceptRate = 0.5;       // target Metropolis acceptance rate
        double stepAdaptRate = 0.05;         // rate of step size adjustment
    };

private:
    Parameters params;
    PotentialBase* potential;
    std::unique_ptr<ThreadPool> threadPool;
    std::vector<std::unique_ptr<NELbfgs>> lbfgsPool;

    // Replica for parallel tempering
    struct Replica {
        BinaryAlloyCluster structure;
        double energy;
        double temperature;
        double stepSize;    // adaptive per-replica step size
        int acceptCount;    // for adaptive step tracking
        int rejectCount;
    };
    std::vector<Replica> replicas;

    // OP-space structural memory
    struct OPMemoryEntry {
        double Q4, Q6;       // Steinhardt order parameters
        int visitCount;
    };
    std::vector<OPMemoryEntry> opMemory;

    // Best found
    BinaryAlloyCluster globalBest;
    double globalBestEnergy;

    int numAtoms, numElementA, numElementB;
    std::string elementA, elementB;

    // Statistics
    std::atomic<int> stepCount{0};
    std::atomic<int> acceptedCount{0};
    std::atomic<int> swapAccepted{0};

public:
    BasinHopping(const Parameters& p, PotentialBase* pot);
    BinaryAlloyCluster optimize(const BinaryAlloyCluster& initial,
        ResultManager* resultManager = nullptr);

    double getBestEnergy() const { return globalBestEnergy; }
    int getStepCount() const { return stepCount.load(); }
    int getAcceptedCount() const { return acceptedCount.load(); }

private:
    void initialize(const BinaryAlloyCluster& initial);
    double evaluateCluster(BinaryAlloyCluster& cluster, NELbfgs* lbfgs);

    // Move operators
    BinaryAlloyCluster singleAtomMove(const BinaryAlloyCluster& src, int replicaIdx);
    BinaryAlloyCluster collectiveMove(const BinaryAlloyCluster& src);
    BinaryAlloyCluster swapMove(const BinaryAlloyCluster& src);

    // Surface bias
    std::vector<int> computeCoordinationNumbers(const BinaryAlloyCluster& cluster) const;
    int selectSurfaceBiasedAtom(const std::vector<int>& coordNumbers) const;

    // Structural memory
    void updateOPMemory(const BinaryAlloyCluster& cluster);
    bool isOPRegionVisited(double Q4, double Q6) const;
    void computeQuickOP(const BinaryAlloyCluster& cluster, double& Q4, double& Q6) const;

    // Adaptive step size
    void updateAdaptiveStep(int replicaIdx);

    // Metropolis
    bool metropolisAccept(double oldE, double newE, double T);

    // Parallel tempering
    void attemptReplicaExchange();
};
