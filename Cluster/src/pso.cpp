#include "../include/pso.h"

// ==================== PSO_Population Implementation ====================

PSO_Population::PSO_Population(const Parameters& p, PotentialBase* pot,
    const BinaryAlloyCluster& initial, ThreadPool* pool)
    : params(p), potential(pot), threadPool(pool),
    currentGeneration(0), omega(params.omegaMax),
    globalBest(initial.getNumElementA(), initial.getNumElementB(),
        initial.getElementA(), initial.getElementB()) {

    numElementA = initial.getNumElementA();
    numElementB = initial.getNumElementB();
    elementA = initial.getElementA();
    elementB = initial.getElementB();
    numAtoms = numElementA + numElementB;

    // Create NELbfgs instances for parallel evaluation
    int numLbfgs = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
    for (int i = 0; i < numLbfgs; ++i) {
        lbfgsPool.push_back(std::make_unique<NELbfgs>(pot));
    }

    initializePopulation(initial);
}

void PSO_Population::initializePopulation(const BinaryAlloyCluster& initial) {
    population.reserve(params.populationSize);
    for (int i = 0; i < params.populationSize; ++i) {
        population.emplace_back(numElementA, numElementB, elementA, elementB);
    }

    population[0].cluster = initial;
    population[0].energy = evaluateCluster(population[0].cluster, lbfgsPool[0].get());
    population[0].pbest = population[0].cluster;
    population[0].pbestEnergy = population[0].energy;

    int countStructureFiles = population[0].cluster.countStructureFiles();
    int minLoadAttempts = (std::min)(10, params.populationSize - 1);

    for (int i = 1; i < params.populationSize; ++i) {
        bool loadSuccess = false;
        if (countStructureFiles > 0 && i <= minLoadAttempts) {
            int fileIndex = ((i - 1) % countStructureFiles) + 1;
            loadSuccess = population[i].cluster.loadStructureInitialize(fileIndex, numAtoms);
        }
        if (!loadSuccess) {
            population[i].cluster.randomInitialize(2.75);
        }
        population[i].energy = evaluateCluster(population[i].cluster, lbfgsPool[0].get());
        population[i].pbest = population[i].cluster;
        population[i].pbestEnergy = population[i].energy;
    }

    initializeVelocities();

    stagnationCounters.assign(params.populationSize, 0);
    globalStagnationCounter = 0;

    // HCLPSO: initialize fitness ranking and exemplars for all particles
    updateFitnessRanking();
    for (int i = 0; i < params.populationSize; ++i) {
        selectExemplars(i);
    }

    globalBest = population[0];
    for (const auto& ind : population) {
        if (ind.energy < globalBest.energy) {
            globalBest = ind;
        }
    }
}

void PSO_Population::initializeVelocities() {
    for (auto& ind : population) {
        for (size_t i = 0; i < ind.velocity.size(); ++i) {
            ind.velocity[i] = RandomGenerator::uniform(-params.vmax, params.vmax);
        }
    }
}

void PSO_Population::updateOmega() {
    omega = params.omegaMax - (params.omegaMax - params.omegaMin) *
        static_cast<double>(currentGeneration) / params.maxGenerations;
}

void PSO_Population::clampVelocity(std::vector<double>& vel) {
    for (auto& v : vel) {
        if (v > params.vmax) v = params.vmax;
        if (v < -params.vmax) v = -params.vmax;
    }
}

double PSO_Population::evaluateCluster(BinaryAlloyCluster& cluster, NELbfgs* lbfgs) {
    int count = evaluationCount.fetch_add(1) + 1;

    if (params.useLocalSearch && count % params.localSearchFrequency == 0) {
        lbfgs->optimize(cluster);
        localSearchCount.fetch_add(1);
    }

    double energy = potential->calculateEnergy(cluster);
    cluster.setEnergy(energy);
    return energy;
}

// ==================== Evolution ====================

void PSO_Population::evolve() {
    if (threadPool) {
        // === Parallel path ===
        generateCandidates();
        evaluateCandidates();
        applyResults();
    } else {
        // === Serial fallback ===
        NELbfgs* lbfgs = lbfgsPool[0].get();
        int numPSO = static_cast<int>(params.populationSize * params.psoRatio);
        int numRandom = params.populationSize - numPSO;

        for (int i = 0; i < numPSO; ++i) {
            auto& ind = population[i];
            auto& coords = ind.cluster.getCoordinates();

            // HCLPSO: all particles use CLPSO exemplar-based velocity update
            // No gbest term — prevents population collapse
            refreshExemplarsIfNeeded(i);
            for (int a = 0; a < numAtoms; ++a) {
                int exemplarIdx = ind.exemplarIndices[a];
                const auto& exemplarCoords = population[exemplarIdx].pbest.getCoordinates();
                for (int d = 0; d < 3; ++d) {
                    int j = 3 * a + d;
                    double r = RandomGenerator::uniform();
                    ind.velocity[j] = omega * ind.velocity[j] +
                        ind.adaptiveC * r * (exemplarCoords[j] - coords[j]);
                }
            }
            clampVelocity(ind.velocity);
            for (size_t j = 0; j < coords.size(); ++j) {
                coords[j] += ind.velocity[j];
            }
            ind.energy = evaluateCluster(ind.cluster, lbfgs);
        }

        // Velocity reset for stagnant particles
        for (int i = 0; i < numPSO; ++i) {
            if (population[i].energy < population[i].pbestEnergy) {
                stagnationCounters[i] = 0;
            } else {
                stagnationCounters[i]++;
                if (stagnationCounters[i] > params.stagnationLimit) {
                    for (auto& v : population[i].velocity) {
                        v = RandomGenerator::uniform(-params.vmax, params.vmax);
                    }
                    stagnationCounters[i] = 0;
                }
            }
        }

        for (int i = 0; i < numRandom; ++i) {
            int idx = numPSO + i;
            if (idx >= params.populationSize) break;
            stagnationCounters[idx] = 0;
            population[idx].cluster.randomInitialize(2.75);
            population[idx].energy = evaluateCluster(population[idx].cluster, lbfgs);
            for (auto& v : population[idx].velocity) {
                v = RandomGenerator::uniform(-params.vmax, params.vmax);
            }
        }

        // Update pbest, adaptive c, and gbest
        double prevBest = globalBest.energy;
        for (auto& ind : population) {
            if (ind.energy < ind.pbestEnergy) {
                ind.pbest = ind.cluster;
                ind.pbestEnergy = ind.energy;
                ind.improvementCount++;
                // Reduce c to exploit around this good region
                ind.adaptiveC = std::max(0.5, ind.adaptiveC * (1.0 - params.cAdaptRate));
            } else {
                ind.improvementCount = 0;
                // Increase c to explore more broadly
                ind.adaptiveC = std::min(2.5, ind.adaptiveC * (1.0 + params.cAdaptRate));
            }
            if (ind.energy < globalBest.energy) {
                globalBest.cluster = ind.cluster;
                globalBest.energy = ind.energy;
                globalBest.pbestEnergy = ind.energy;
            }
        }

        // Swap atoms for binary alloy
        if (numElementA > 0 && numElementA < numAtoms) {
            for (int idx = 0; idx < params.populationSize; ++idx) {
                if (fabs(population[idx].energy - globalBest.energy) < 0.2 &&
                    RandomGenerator::uniform() < 0.9) {
                    BinaryAlloyCluster newCluster = population[idx].cluster;
                    std::vector<int> typeA, typeB;
                    for (int i = 0; i < numAtoms; ++i) {
                        if (newCluster.getAtomType(i) == 0) typeA.push_back(i);
                        else typeB.push_back(i);
                    }
                    if (!typeA.empty() && !typeB.empty()) {
                        int idxA = typeA[RandomGenerator::uniformInt(0, typeA.size() - 1)];
                        int idxB = typeB[RandomGenerator::uniformInt(0, typeB.size() - 1)];
                        auto posA = newCluster.getAtomPosition(idxA);
                        auto posB = newCluster.getAtomPosition(idxB);
                        newCluster.setAtomPosition(idxA, posB[0], posB[1], posB[2]);
                        newCluster.setAtomPosition(idxB, posA[0], posA[1], posA[2]);
                        double newEnergy = evaluateCluster(newCluster, lbfgs);
                        if (newEnergy < population[idx].energy) {
                            population[idx].cluster = newCluster;
                            population[idx].energy = newEnergy;
                            if (newEnergy < population[idx].pbestEnergy) {
                                population[idx].pbest = newCluster;
                                population[idx].pbestEnergy = newEnergy;
                            }
                            if (newEnergy < globalBest.energy) {
                                globalBest.cluster = newCluster;
                                globalBest.energy = newEnergy;
                            }
                        }
                    }
                }
            }
        }

        // Global best perturbation
        if (globalBest.energy < prevBest) {
            globalStagnationCounter = 0;
        } else {
            globalStagnationCounter++;
            if (globalStagnationCounter > params.stagnationThreshold) {
                perturbGlobalBest();
                globalStagnationCounter = 0;
            }
        }
    }

    currentGeneration++;
    updateOmega();
}

// ==================== Parallel evolution phases ====================

void PSO_Population::generateCandidates() {
    int numPSO = static_cast<int>(params.populationSize * params.psoRatio);
    int numRandom = params.populationSize - numPSO;

    evalTargets.clear();
    evalTargets.reserve(params.populationSize);

    // HCLPSO updates: velocity + position (no evaluation), all particles use CLPSO
    for (int i = 0; i < numPSO; ++i) {
        // Reset velocity for stagnant particles
        if (stagnationCounters[i] > params.stagnationLimit) {
            for (auto& v : population[i].velocity) {
                v = RandomGenerator::uniform(-params.vmax, params.vmax);
            }
            stagnationCounters[i] = 0;
        }

        auto& ind = population[i];
        auto& coords = ind.cluster.getCoordinates();

        // HCLPSO: per-atom exemplar velocity update with self-adaptive c
        refreshExemplarsIfNeeded(i);
        for (int a = 0; a < numAtoms; ++a) {
            int exemplarIdx = ind.exemplarIndices[a];
            const auto& exemplarCoords = population[exemplarIdx].pbest.getCoordinates();
            for (int d = 0; d < 3; ++d) {
                int j = 3 * a + d;
                double r = RandomGenerator::uniform();
                ind.velocity[j] = omega * ind.velocity[j] +
                    ind.adaptiveC * r * (exemplarCoords[j] - coords[j]);
            }
        }
        clampVelocity(ind.velocity);
        for (size_t j = 0; j < coords.size(); ++j) {
            coords[j] += ind.velocity[j];
        }

        evalTargets.push_back({&ind.cluster, &ind.energy});
    }

    // Random structures: generate (no evaluation)
    for (int i = 0; i < numRandom; ++i) {
        int idx = numPSO + i;
        if (idx >= params.populationSize) break;
        population[idx].cluster.randomInitialize(2.75);
        for (auto& v : population[idx].velocity) {
            v = RandomGenerator::uniform(-params.vmax, params.vmax);
        }
        evalTargets.push_back({&population[idx].cluster, &population[idx].energy});
    }
}

void PSO_Population::evaluateCandidates() {
    if (evalTargets.empty()) return;

    int poolSize = static_cast<int>(lbfgsPool.size());
    std::atomic<int> lbfgsIdx{0};

    for (auto& target : evalTargets) {
        threadPool->submit([this, &target, &lbfgsIdx, poolSize]() {
            int idx = lbfgsIdx.fetch_add(1) % poolSize;
            *target.energy = evaluateCluster(*target.cluster, lbfgsPool[idx].get());
        });
    }
    threadPool->waitAll();
}

void PSO_Population::applyResults() {
    double prevBest = globalBest.energy;

    // Update pbest, adaptive c, and gbest
    int numPSO = static_cast<int>(params.populationSize * params.psoRatio);
    for (int i = 0; i < params.populationSize; ++i) {
        if (population[i].energy < population[i].pbestEnergy) {
            population[i].pbest = population[i].cluster;
            population[i].pbestEnergy = population[i].energy;
            if (i < numPSO) stagnationCounters[i] = 0;
            // Adaptive c: reduce to exploit good region
            population[i].improvementCount++;
            population[i].adaptiveC = std::max(0.5,
                population[i].adaptiveC * (1.0 - params.cAdaptRate));
            // CLPSO: pbest improved, reset exemplar refresh counter
            population[i].stagnationSinceRefresh = 0;
        } else if (i < numPSO) {
            stagnationCounters[i]++;
            population[i].improvementCount = 0;
            // Adaptive c: increase to explore more
            population[i].adaptiveC = std::min(2.5,
                population[i].adaptiveC * (1.0 + params.cAdaptRate));
            // CLPSO: track stagnation for exemplar refresh
            population[i].stagnationSinceRefresh++;
        }
        if (population[i].energy < globalBest.energy) {
            globalBest.cluster = population[i].cluster;
            globalBest.energy = population[i].energy;
            globalBest.pbestEnergy = population[i].energy;
        }
    }

    // Update fitness ranking after pbest changes
    updateFitnessRanking();

    // Random structure particles always reset stagnation
    for (int i = numPSO; i < params.populationSize; ++i) {
        stagnationCounters[i] = 0;
    }

    // Swap atoms for binary alloy
    if (numElementA > 0 && numElementA < numAtoms) {
        swapAtoms();
    }

    // Global best perturbation
    if (globalBest.energy < prevBest) {
        globalStagnationCounter = 0;
    } else {
        globalStagnationCounter++;
        if (globalStagnationCounter > params.stagnationThreshold) {
            perturbGlobalBest();
            globalStagnationCounter = 0;
        }
    }
}

void PSO_Population::perturbGlobalBest() {
    int poolSize = static_cast<int>(lbfgsPool.size());
    std::atomic<int> lbfgsIdx{0};

    // Variable neighborhood: try multiple scales from small to large
    // Small perturbations exploit current basin, large ones explore new basins
    std::vector<double> scales = {
        params.perturbationScale * 0.2,   // small: fine refinement
        params.perturbationScale * 0.5,   // medium: local exploration
        params.perturbationScale,          // normal: basin exploration
        params.perturbationScale * 2.0    // large: jump to new basin
    };

    int perScale = params.numPerturbations / static_cast<int>(scales.size());
    if (perScale < 1) perScale = 1;

    std::vector<BinaryAlloyCluster> perturbed;
    std::vector<double> energies;
    perturbed.reserve(params.numPerturbations);
    energies.reserve(params.numPerturbations);

    for (double scale : scales) {
        for (int p = 0; p < perScale; ++p) {
            BinaryAlloyCluster copy = globalBest.cluster;
            auto& coords = copy.getCoordinates();
            for (size_t j = 0; j < coords.size(); ++j) {
                coords[j] += RandomGenerator::normal(0.0, scale);
            }
            perturbed.push_back(std::move(copy));
            energies.push_back(globalBest.energy);
        }
    }

    int total = static_cast<int>(perturbed.size());

    // Evaluate all in parallel
    if (threadPool && total > 0) {
        for (int p = 0; p < total; ++p) {
            threadPool->submit([this, &perturbed, &energies, &lbfgsIdx, poolSize, p]() {
                int idx = lbfgsIdx.fetch_add(1) % poolSize;
                energies[p] = evaluateCluster(perturbed[p], lbfgsPool[idx].get());
            });
        }
        threadPool->waitAll();
    } else {
        NELbfgs* lbfgs = lbfgsPool[0].get();
        for (int p = 0; p < total; ++p) {
            energies[p] = evaluateCluster(perturbed[p], lbfgs);
        }
    }

    // Find best perturbed structure
    for (int p = 0; p < total; ++p) {
        if (energies[p] < globalBest.energy) {
            globalBest.cluster = std::move(perturbed[p]);
            globalBest.energy = energies[p];
            globalBest.pbestEnergy = energies[p];
        }
    }
}

void PSO_Population::swapAtoms() {
    if (threadPool) {
        // Parallel path
        int poolSize = static_cast<int>(lbfgsPool.size());
        std::atomic<int> lbfgsIdx{0};

        struct SwapJob {
            BinaryAlloyCluster newCluster;
            int popIndex;
            double energy;
            bool valid;
        };

        std::vector<SwapJob> swapJobs;
        swapJobs.reserve(params.populationSize);

        for (int idx = 0; idx < params.populationSize; ++idx) {
            if (fabs(population[idx].energy - globalBest.energy) < 0.2 &&
                RandomGenerator::uniform() < 0.9) {
                BinaryAlloyCluster newCluster = population[idx].cluster;
                std::vector<int> typeA, typeB;
                for (int i = 0; i < numAtoms; ++i) {
                    if (newCluster.getAtomType(i) == 0) typeA.push_back(i);
                    else typeB.push_back(i);
                }
                if (!typeA.empty() && !typeB.empty()) {
                    int idxA = typeA[RandomGenerator::uniformInt(0, typeA.size() - 1)];
                    int idxB = typeB[RandomGenerator::uniformInt(0, typeB.size() - 1)];
                    auto posA = newCluster.getAtomPosition(idxA);
                    auto posB = newCluster.getAtomPosition(idxB);
                    newCluster.setAtomPosition(idxA, posB[0], posB[1], posB[2]);
                    newCluster.setAtomPosition(idxB, posA[0], posA[1], posA[2]);
                    swapJobs.push_back({newCluster, idx, 0.0, true});
                }
            }
        }

        if (!swapJobs.empty()) {
            for (auto& job : swapJobs) {
                threadPool->submit([this, &job, &lbfgsIdx, poolSize]() {
                    int id = lbfgsIdx.fetch_add(1) % poolSize;
                    job.energy = evaluateCluster(job.newCluster, lbfgsPool[id].get());
                });
            }
            threadPool->waitAll();
        }

        for (auto& job : swapJobs) {
            if (job.valid && job.energy < population[job.popIndex].energy) {
                population[job.popIndex].cluster = job.newCluster;
                population[job.popIndex].energy = job.energy;
                if (job.energy < population[job.popIndex].pbestEnergy) {
                    population[job.popIndex].pbest = job.newCluster;
                    population[job.popIndex].pbestEnergy = job.energy;
                }
                if (job.energy < globalBest.energy) {
                    globalBest.cluster = job.newCluster;
                    globalBest.energy = job.energy;
                }
            }
        }
    } else {
        // Serial fallback
        NELbfgs* lbfgs = lbfgsPool[0].get();
        for (int idx = 0; idx < params.populationSize; ++idx) {
            if (fabs(population[idx].energy - globalBest.energy) < 0.2 &&
                RandomGenerator::uniform() < 0.9) {
                BinaryAlloyCluster newCluster = population[idx].cluster;
                std::vector<int> typeA, typeB;
                for (int i = 0; i < numAtoms; ++i) {
                    if (newCluster.getAtomType(i) == 0) typeA.push_back(i);
                    else typeB.push_back(i);
                }
                if (!typeA.empty() && !typeB.empty()) {
                    int idxA = typeA[RandomGenerator::uniformInt(0, typeA.size() - 1)];
                    int idxB = typeB[RandomGenerator::uniformInt(0, typeB.size() - 1)];
                    auto posA = newCluster.getAtomPosition(idxA);
                    auto posB = newCluster.getAtomPosition(idxB);
                    newCluster.setAtomPosition(idxA, posB[0], posB[1], posB[2]);
                    newCluster.setAtomPosition(idxB, posA[0], posA[1], posA[2]);
                    double newEnergy = evaluateCluster(newCluster, lbfgs);
                    if (newEnergy < population[idx].energy) {
                        population[idx].cluster = newCluster;
                        population[idx].energy = newEnergy;
                        if (newEnergy < population[idx].pbestEnergy) {
                            population[idx].pbest = newCluster;
                            population[idx].pbestEnergy = newEnergy;
                        }
                        if (newEnergy < globalBest.energy) {
                            globalBest.cluster = newCluster;
                            globalBest.energy = newEnergy;
                        }
                    }
                }
            }
        }
    }
}

// ==================== CLPSO Methods ====================

void PSO_Population::updateFitnessRanking() {
    // Sort particle indices by pbest energy (ascending = better fitness)
    std::vector<std::pair<double, int>> rankData;
    rankData.reserve(population.size());
    for (int i = 0; i < static_cast<int>(population.size()); ++i) {
        rankData.push_back({population[i].pbestEnergy, i});
    }
    std::sort(rankData.begin(), rankData.end());

    fitnessRank.resize(population.size());
    for (int r = 0; r < static_cast<int>(rankData.size()); ++r) {
        fitnessRank[rankData[r].second] = r; // rank 0 = best
    }
}

void PSO_Population::selectExemplars(int particleIdx) {
    int ps = static_cast<int>(population.size());

    // Learning probability: Pc = 0.05 + 0.45 * (exp(10*(rank)/(ps-1)) - 1) / (exp(10) - 1)
    // Better-ranked particles (low rank) have lower Pc → more self-reliance
    double rank = static_cast<double>(fitnessRank[particleIdx]);
    double pc = 0.05 + 0.45 * (std::exp(10.0 * rank / (ps - 1)) - 1.0) / (std::exp(10.0) - 1.0);

    auto& exemplars = population[particleIdx].exemplarIndices;

    for (int a = 0; a < numAtoms; ++a) {
        if (RandomGenerator::uniform() < pc) {
            if (params.useRingTopology) {
                // Ring topology: pick exemplar from neighbors (i-1, i, i+1)
                int cands[3] = {
                    (particleIdx - 1 + ps) % ps,
                    particleIdx,
                    (particleIdx + 1) % ps
                };
                // Pick the one with best pbest energy
                int bestCand = cands[0];
                double bestE = population[cands[0]].pbestEnergy;
                for (int c = 1; c < 3; ++c) {
                    if (population[cands[c]].pbestEnergy < bestE) {
                        bestE = population[cands[c]].pbestEnergy;
                        bestCand = cands[c];
                    }
                }
                exemplars[a] = bestCand;
            } else {
                // Tournament: randomly pick 2, choose the one with better pbest energy
                int r1 = RandomGenerator::uniformInt(0, ps - 1);
                int r2 = RandomGenerator::uniformInt(0, ps - 1);
                exemplars[a] = (population[r1].pbestEnergy < population[r2].pbestEnergy) ? r1 : r2;
            }
        } else {
            // Learn from own pbest
            exemplars[a] = particleIdx;
        }
    }
}

void PSO_Population::refreshExemplarsIfNeeded(int particleIdx) {
    auto& ind = population[particleIdx];

    if (ind.stagnationSinceRefresh >= params.exemplarGap) {
        selectExemplars(particleIdx);
        ind.stagnationSinceRefresh = 0;
    }
}

void PSO_Population::receiveIndividual(const PSO_Individual& ind) {
    int worstIdx = 0;
    double worstEnergy = population[0].energy;
    for (int i = 1; i < params.populationSize; ++i) {
        if (population[i].energy > worstEnergy) {
            worstEnergy = population[i].energy;
            worstIdx = i;
        }
    }
    if (ind.energy < worstEnergy) {
        population[worstIdx] = ind;
        if (ind.energy < globalBest.energy) {
            globalBest = ind;
        }
    }
}

// ==================== PSO Main Class Implementation ====================

PSO::PSO(const Parameters& p, PotentialBase* pot)
    : params(p), potential(pot), generation(0),
    globalBest(1, 0, "A", "B") {

    exchangeInterval = 50;

    if (params.useThreading) {
        int nThreads = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
        threadPool = std::make_unique<ThreadPool>(nThreads);
    }
}

void PSO::initialize(const BinaryAlloyCluster& initial) {
    populations.clear();

    ThreadPool* pool = threadPool.get();

    for (int i = 0; i < NUM_POPULATIONS; ++i) {
        populations.emplace_back(std::make_unique<PSO_Population>(
            params, potential, initial, pool));
    }

    globalBest = populations[0]->getGlobalBest();
    for (const auto& pop : populations) {
        if (pop->getGlobalBest().energy < globalBest.energy) {
            globalBest = pop->getGlobalBest();
        }
    }
}

void PSO::updateGlobalBest(const PSO_Individual& candidate) {
    std::lock_guard<std::mutex> lock(globalBestMutex);
    if (candidate.energy < globalBest.energy) {
        globalBest = candidate;
    }
}

void PSO::evolve() {
    generation++;

    for (auto& pop : populations) {
        pop->evolve();
        updateGlobalBest(pop->getGlobalBest());
    }

    if (generation % exchangeInterval == 0) {
        exchangeBestIndividuals();
    }
}

void PSO::exchangeBestIndividuals() {
    for (size_t i = 0; i < populations.size(); ++i) {
        const auto& bestInd = populations[i]->getGlobalBest();
        for (size_t j = 0; j < populations.size(); ++j) {
            if (i != j) {
                populations[j]->receiveIndividual(bestInd);
            }
        }
    }
}

BinaryAlloyCluster PSO::optimize(const BinaryAlloyCluster& initial, ResultManager* resultManager) {
    initialize(initial);

    for (int gen = 0; gen < params.maxGenerations; ++gen) {
        evolve();

        if (resultManager) {
            resultManager->saveGenerationBest(generation, globalBest.cluster);
            resultManager->updateHistoricalBest(globalBest.cluster);

            std::vector<BinaryAlloyCluster> dummyPop;
            for (const auto& pop : populations) {
                dummyPop.push_back(pop->getGlobalBest().cluster);
            }
            resultManager->logGeneration(generation, dummyPop, globalBest.cluster);
        }

        if (gen % 1 == 0) {
            std::cout << "Gen " << gen << ": Best = " << globalBest.energy;
            if (populations.size() > 1) {
                std::cout << " [";
                for (size_t i = 0; i < populations.size(); ++i) {
                    if (i > 0) std::cout << ", ";
                    std::cout << populations[i]->getGlobalBest().energy;
                }
                std::cout << "]";
            }
            if (resultManager && resultManager->hasHistory()) {
                std::cout << ", Historical = " << resultManager->getHistoricalBestEnergy();
            }
            std::cout << std::endl;
        }
    }

    return globalBest.cluster;
}
