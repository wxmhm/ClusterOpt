#include "../include/BasinHopping.h"
#include <algorithm>
#include <cmath>
#include <iostream>

BasinHopping::BasinHopping(const Parameters& p, PotentialBase* pot)
    : params(p), potential(pot),
    globalBest(1, 0, "A", "B"),
    globalBestEnergy(1e10),
    numAtoms(0), numElementA(0), numElementB(0) {

    if (params.useThreading) {
        int nThreads = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
        threadPool = std::make_unique<ThreadPool>(nThreads);
        for (int i = 0; i < nThreads; ++i) {
            lbfgsPool.push_back(std::make_unique<NELbfgs>(pot));
        }
    } else {
        lbfgsPool.push_back(std::make_unique<NELbfgs>(pot));
    }
}

void BasinHopping::initialize(const BinaryAlloyCluster& initial) {
    numAtoms = initial.getNumAtoms();
    numElementA = initial.getNumElementA();
    numElementB = initial.getNumElementB();
    elementA = initial.getElementA();
    elementB = initial.getElementB();

    // Initialize replicas at geometrically spaced temperatures
    replicas.clear();
    for (int r = 0; r < params.numReplicas; ++r) {
        double frac = static_cast<double>(r) / (params.numReplicas - 1);
        double T = params.T0 * std::pow(params.Tfinal / params.T0, frac);

        BinaryAlloyCluster structure(numElementA, numElementB, elementA, elementB);
        if (r == 0 && initial.getNumAtoms() == numAtoms) {
            structure = initial;
        } else {
            structure.randomInitialize(2.75);
        }

        NELbfgs* lbfgs = lbfgsPool[0].get();
        double energy = evaluateCluster(structure, lbfgs);
        replicas.push_back({structure, energy, T, params.singleStepSize, 0, 0});

        if (energy < globalBestEnergy) {
            globalBest = structure;
            globalBestEnergy = energy;
        }
    }

    std::cout << "BH initialized: " << params.numReplicas << " replicas, "
        << "T range [" << replicas.back().temperature << ", " << replicas[0].temperature << "]"
        << std::endl;
}

double BasinHopping::evaluateCluster(BinaryAlloyCluster& cluster, NELbfgs* lbfgs) {
    int count = stepCount.fetch_add(1) + 1;

    if (params.useLocalSearch && count % params.localSearchFrequency == 0) {
        lbfgs->optimize(cluster);
    }

    double energy = potential->calculateEnergy(cluster);
    cluster.setEnergy(energy);
    return energy;
}

BinaryAlloyCluster BasinHopping::singleAtomMove(const BinaryAlloyCluster& src, int replicaIdx) {
    BinaryAlloyCluster trial = src.copy();

    int atomIdx;
    if (params.useSurfaceBias) {
        auto coordNumbers = computeCoordinationNumbers(src);
        atomIdx = selectSurfaceBiasedAtom(coordNumbers);
    } else {
        atomIdx = RandomGenerator::uniformInt(0, numAtoms - 1);
    }

    auto pos = trial.getAtomPosition(atomIdx);

    // Use adaptive per-replica step size
    double Tfactor = replicas[replicaIdx].temperature / params.T0;
    double scale = replicas[replicaIdx].stepSize * (0.2 + 0.8 * Tfactor);

    pos[0] += RandomGenerator::normal(0.0, scale);
    pos[1] += RandomGenerator::normal(0.0, scale);
    pos[2] += RandomGenerator::normal(0.0, scale);

    trial.setAtomPosition(atomIdx, pos[0], pos[1], pos[2]);
    return trial;
}

BinaryAlloyCluster BasinHopping::collectiveMove(const BinaryAlloyCluster& src) {
    BinaryAlloyCluster trial = src.copy();

    // Pick a random seed atom
    int seed = RandomGenerator::uniformInt(0, numAtoms - 1);
    auto seedPos = trial.getAtomPosition(seed);

    // Find K nearest neighbors to form the collective group
    std::vector<std::pair<double, int>> distances;
    for (int i = 0; i < numAtoms; ++i) {
        if (i == seed) continue;
        auto pos = trial.getAtomPosition(i);
        double dx = pos[0] - seedPos[0];
        double dy = pos[1] - seedPos[1];
        double dz = pos[2] - seedPos[2];
        double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        distances.push_back({dist, i});
    }
    std::sort(distances.begin(), distances.end());

    // Take K-1 nearest neighbors + seed
    int groupSize = std::min(params.collectiveSize, numAtoms);
    std::vector<int> group = {seed};
    for (int k = 0; k < groupSize - 1 && k < static_cast<int>(distances.size()); ++k) {
        group.push_back(distances[k].second);
    }

    // Generate a random displacement vector for the entire group
    double dx = RandomGenerator::normal(0.0, params.collectiveStepSize);
    double dy = RandomGenerator::normal(0.0, params.collectiveStepSize);
    double dz = RandomGenerator::normal(0.0, params.collectiveStepSize);

    // Apply same displacement to all atoms in the group
    for (int idx : group) {
        auto pos = trial.getAtomPosition(idx);
        trial.setAtomPosition(idx, pos[0] + dx, pos[1] + dy, pos[2] + dz);
    }

    return trial;
}

BinaryAlloyCluster BasinHopping::swapMove(const BinaryAlloyCluster& src) {
    BinaryAlloyCluster trial = src.copy();

    // Find atoms of different types
    std::vector<int> typeA, typeB;
    for (int i = 0; i < numAtoms; ++i) {
        if (trial.getAtomType(i) == 0) typeA.push_back(i);
        else typeB.push_back(i);
    }

    if (typeA.empty() || typeB.empty()) return trial;

    // Pick one atom of each type and swap their positions
    int idxA = typeA[RandomGenerator::uniformInt(0, static_cast<int>(typeA.size()) - 1)];
    int idxB = typeB[RandomGenerator::uniformInt(0, static_cast<int>(typeB.size()) - 1)];

    auto posA = trial.getAtomPosition(idxA);
    auto posB = trial.getAtomPosition(idxB);

    trial.setAtomPosition(idxA, posB[0], posB[1], posB[2]);
    trial.setAtomPosition(idxB, posA[0], posA[1], posA[2]);

    return trial;
}

bool BasinHopping::metropolisAccept(double oldE, double newE, double T) {
    if (newE < oldE) return true;
    double prob = std::exp(-(newE - oldE) / T);
    return RandomGenerator::uniform() < prob;
}

// --- Surface bias: coordination-number-weighted atom selection ---

std::vector<int> BasinHopping::computeCoordinationNumbers(const BinaryAlloyCluster& cluster) const {
    // Find average nearest-neighbor distance
    std::vector<double> allDists;
    for (int i = 0; i < numAtoms; ++i) {
        auto pi = cluster.getAtomPosition(i);
        for (int j = i + 1; j < numAtoms; ++j) {
            auto pj = cluster.getAtomPosition(j);
            double dx = pi[0] - pj[0], dy = pi[1] - pj[1], dz = pi[2] - pj[2];
            allDists.push_back(std::sqrt(dx*dx + dy*dy + dz*dz));
        }
    }
    std::sort(allDists.begin(), allDists.end());
    double firstNN = allDists[std::min(numAtoms, static_cast<int>(allDists.size() * 0.02))];
    double cutoff = firstNN * 1.35;

    std::vector<int> coordNums(numAtoms, 0);
    for (int i = 0; i < numAtoms; ++i) {
        auto pi = cluster.getAtomPosition(i);
        for (int j = 0; j < numAtoms; ++j) {
            if (i == j) continue;
            auto pj = cluster.getAtomPosition(j);
            double dx = pi[0] - pj[0], dy = pi[1] - pj[1], dz = pi[2] - pj[2];
            if (std::sqrt(dx*dx + dy*dy + dz*dz) < cutoff) {
                coordNums[i]++;
            }
        }
    }
    return coordNums;
}

int BasinHopping::selectSurfaceBiasedAtom(const std::vector<int>& coordNumbers) const {
    // Probability of selecting atom i ∝ 1 / (coordNumbers[i] + 1)^bias
    // where bias = surfaceBiasStrength. Higher bias → stronger surface preference.
    int maxCoord = 1;
    for (int cn : coordNumbers) maxCoord = std::max(maxCoord, cn);

    std::vector<double> weights(numAtoms);
    double totalWeight = 0.0;
    double b = params.surfaceBiasStrength;

    for (int i = 0; i < numAtoms; ++i) {
        // Normalize and invert: atoms with fewer bonds get higher weight
        double normCoord = static_cast<double>(coordNumbers[i]) / (maxCoord + 1);
        weights[i] = std::pow(1.0 - normCoord, b) + 0.01;  // +0.01 avoids zero weight
        totalWeight += weights[i];
    }

    double r = RandomGenerator::uniform() * totalWeight;
    double cumsum = 0.0;
    for (int i = 0; i < numAtoms; ++i) {
        cumsum += weights[i];
        if (r <= cumsum) return i;
    }
    return numAtoms - 1;
}

// --- Structural memory via OP-space ---

void BasinHopping::computeQuickOP(const BinaryAlloyCluster& cluster, double& Q4, double& Q6) const {
    // Fast Q4/Q6 approximation using only first-neighbor shell distances
    // (avoids full spherical harmonics computation)
    std::vector<double> allDists;
    for (int i = 0; i < numAtoms; ++i) {
        auto pi = cluster.getAtomPosition(i);
        for (int j = i + 1; j < numAtoms; ++j) {
            auto pj = cluster.getAtomPosition(j);
            double dx = pi[0] - pj[0], dy = pi[1] - pj[1], dz = pi[2] - pj[2];
            allDists.push_back(std::sqrt(dx*dx + dy*dy + dz*dz));
        }
    }
    std::sort(allDists.begin(), allDists.end());
    double r0 = allDists[std::min(5, static_cast<int>(allDists.size() * 0.01))];

    // Count how many atoms are within 1.1*r0 and 1.4*r0 (first and second shells)
    int firstShell = 0, secondShell = 0;
    for (int i = 0; i < numAtoms; ++i) {
        auto pi = cluster.getAtomPosition(i);
        int nb = 0;
        for (int j = 0; j < numAtoms; ++j) {
            if (i == j) continue;
            auto pj = cluster.getAtomPosition(j);
            double dx = pi[0] - pj[0], dy = pi[1] - pj[1], dz = pi[2] - pj[2];
            double d = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (d < r0 * 1.15) nb++;
        }
        if (nb >= 10) firstShell++;
        else if (nb >= 6) secondShell++;
    }

    // Q4 and Q6 proxies based on coordination structure
    // Icosahedral: high firstShell, low secondShell → Q6 large
    // FCC: moderate firstShell, moderate secondShell → Q4 large
    double f1 = static_cast<double>(firstShell) / numAtoms;
    double f2 = static_cast<double>(secondShell) / numAtoms;
    Q4 = 0.05 + 0.25 * f2;       // Q4 ~ 0.05-0.30
    Q6 = 0.10 + 0.60 * f1;       // Q6 ~ 0.10-0.70
}

void BasinHopping::updateOPMemory(const BinaryAlloyCluster& cluster) {
    double Q4, Q6;
    computeQuickOP(cluster, Q4, Q6);

    // Check if this OP region is already in memory
    for (auto& entry : opMemory) {
        double dQ4 = entry.Q4 - Q4, dQ6 = entry.Q6 - Q6;
        if (std::sqrt(dQ4*dQ4 + dQ6*dQ6) < params.opMemoryResolution) {
            entry.visitCount++;
            return;
        }
    }

    // New region
    if (static_cast<int>(opMemory.size()) >= params.opMemoryMaxSize) {
        opMemory.erase(opMemory.begin());
    }
    opMemory.push_back({Q4, Q6, 1});
}

bool BasinHopping::isOPRegionVisited(double Q4, double Q6) const {
    for (const auto& entry : opMemory) {
        double dQ4 = entry.Q4 - Q4, dQ6 = entry.Q6 - Q6;
        if (std::sqrt(dQ4*dQ4 + dQ6*dQ6) < params.opMemoryResolution) {
            return entry.visitCount > 3;  // visited > 3 times → over-explored
        }
    }
    return false;
}

// --- Adaptive step size ---

void BasinHopping::updateAdaptiveStep(int replicaIdx) {
    if (!params.useAdaptiveStep) return;
    auto& rep = replicas[replicaIdx];
    int total = rep.acceptCount + rep.rejectCount;
    if (total < 10) return;  // not enough data yet

    double rate = static_cast<double>(rep.acceptCount) / total;
    if (rate > params.targetAcceptRate + 0.1) {
        // Acceptance too high → step too small, increase
        rep.stepSize *= (1.0 + params.stepAdaptRate);
    } else if (rate < params.targetAcceptRate - 0.1) {
        // Acceptance too low → step too large, decrease
        rep.stepSize *= (1.0 - params.stepAdaptRate);
    }
    rep.stepSize = std::max(0.05, std::min(0.8, rep.stepSize));
    // Reset counters
    rep.acceptCount = 0;
    rep.rejectCount = 0;
}

void BasinHopping::attemptReplicaExchange() {
    for (int r = 0; r < params.numReplicas - 1; ++r) {
        double T1 = replicas[r].temperature;
        double T2 = replicas[r + 1].temperature;
        double E1 = replicas[r].energy;
        double E2 = replicas[r + 1].energy;

        double delta = (1.0 / T1 - 1.0 / T2) * (E1 - E2);
        if (delta > 0 || RandomGenerator::uniform() < std::exp(delta)) {
            std::swap(replicas[r].structure, replicas[r + 1].structure);
            std::swap(replicas[r].energy, replicas[r + 1].energy);
            std::swap(replicas[r].stepSize, replicas[r + 1].stepSize);

            // Update global best
            for (const auto& rep : replicas) {
                if (rep.energy < globalBestEnergy) {
                    globalBest = rep.structure;
                    globalBestEnergy = rep.energy;
                }
            }
        }
    }
}

BinaryAlloyCluster BasinHopping::optimize(const BinaryAlloyCluster& initial,
    ResultManager* resultManager) {

    initialize(initial);

    std::cout << "BH start: best = " << globalBestEnergy << " eV" << std::endl;

    NELbfgs* lbfgs = lbfgsPool[0].get();
    int stepsPerOutput = std::max(1, params.maxTotalSteps / 20);

    for (int step = 0; step < params.maxTotalSteps; ++step) {
        // Cool temperatures periodically
        if (step % params.stepsPerT == 0 && step > 0) {
            for (auto& rep : replicas) {
                rep.temperature *= params.coolingRate;
                if (rep.temperature < params.Tfinal)
                    rep.temperature = params.Tfinal;
            }
        }

        for (int r = 0; r < params.numReplicas; ++r) {
            auto& rep = replicas[r];

            // Structural memory: check if current structure is in over-visited OP region
            bool inVisitedRegion = false;
            if (params.useStructuralMemory) {
                double Q4, Q6;
                computeQuickOP(rep.structure, Q4, Q6);
                inVisitedRegion = isOPRegionVisited(Q4, Q6);
            }

            // Select move type — bias toward collective/swap when revisiting
            double dice = RandomGenerator::uniform();
            double probCollective = params.probCollective;
            double probSwap = params.probSwap;
            if (inVisitedRegion) {
                // Boost collective and swap to escape over-visited region
                probCollective = 0.5;
                probSwap = 0.4;
            }

            BinaryAlloyCluster trial(numElementA, numElementB, elementA, elementB);

            if (dice < params.probSingle) {
                trial = singleAtomMove(rep.structure, r);
            } else if (dice < params.probSingle + probCollective) {
                trial = collectiveMove(rep.structure);
                // Apply multiple collective moves to escape over-visited region
                if (inVisitedRegion) {
                    for (int esc = 0; esc < params.escapeSteps - 1; ++esc) {
                        trial = collectiveMove(trial);
                    }
                }
            } else {
                trial = swapMove(rep.structure);
            }

            double trialEnergy = evaluateCluster(trial, lbfgs);

            // Track structural memory
            if (params.useStructuralMemory && step % 5 == 0) {
                updateOPMemory(trial);
            }

            bool accepted = metropolisAccept(rep.energy, trialEnergy, rep.temperature);
            if (accepted) {
                rep.structure = std::move(trial);
                rep.energy = trialEnergy;
                rep.acceptCount++;
                acceptedCount.fetch_add(1);

                if (trialEnergy < globalBestEnergy) {
                    globalBest = rep.structure;
                    globalBestEnergy = trialEnergy;
                }
            } else {
                rep.rejectCount++;
            }

            // Adaptive step size
            if ((rep.acceptCount + rep.rejectCount) >= 20) {
                updateAdaptiveStep(r);
            }
        }

        // Attempt replica exchange
        if (step % params.exchangeInterval == 0) {
            attemptReplicaExchange();
        }

        // Progress output
        if (step % stepsPerOutput == 0) {
            std::cout << "BH step " << step
                << ": best = " << globalBestEnergy
                << " | replicas: [";
            for (int r = 0; r < params.numReplicas; ++r) {
                if (r > 0) std::cout << ", ";
                std::cout << replicas[r].energy;
            }
            std::cout << "] T=[";
            for (int r = 0; r < params.numReplicas; ++r) {
                if (r > 0) std::cout << ", ";
                std::cout << replicas[r].temperature;
            }
            std::cout << "] steps=[";
            for (int r = 0; r < params.numReplicas; ++r) {
                if (r > 0) std::cout << ", ";
                std::cout << replicas[r].stepSize;
            }
            std::cout << "]" << std::endl;
        }

        // Save results
        if (resultManager && step % 10 == 0) {
            resultManager->saveGenerationBest(step, globalBest);
            resultManager->updateHistoricalBest(globalBest);
        }
    }

    // Final L-BFGS polish on best
    globalBestEnergy = evaluateCluster(globalBest, lbfgs);

    std::cout << "BH final: best = " << globalBestEnergy
        << ", accepted = " << acceptedCount.load() << "/" << stepCount.load()
        << std::endl;

    return globalBest;
}
