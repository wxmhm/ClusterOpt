#include "../include/CDE.h"

// ==================== CDE_Population Implementation ====================

CDE_Population::CDE_Population(CDE_MutationStrategy strat, int popSize,
    const BinaryAlloyCluster& initial,
    PotentialBase* pot, ThreadPool* pool)
    : strategy(strat), populationSize(popSize), potential(pot),
    threadPool(pool), localSearchCount(0),
    bestIndividual(initial.getNumElementA(), initial.getNumElementB(),
        initial.getElementA(), initial.getElementB()) {

    numElementA = initial.getNumElementA();
    numElementB = initial.getNumElementB();
    elementA = initial.getElementA();
    elementB = initial.getElementB();
    numAtoms = numElementA + numElementB;

    // Create NELbfgs instances for parallel evaluation
    // Use at least 1, up to the number of hardware threads
    int numLbfgs = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
    for (int i = 0; i < numLbfgs; ++i) {
        lbfgsPool.push_back(std::make_unique<NELbfgs>(pot));
    }

    population.reserve(populationSize);
    mutantPopulation.reserve(populationSize);
    trialPopulation.reserve(populationSize);

    for (int i = 0; i < populationSize; ++i) {
        population.emplace_back(numElementA, numElementB, elementA, elementB);
        mutantPopulation.emplace_back(numElementA, numElementB, elementA, elementB);
        trialPopulation.emplace_back(numElementA, numElementB, elementA, elementB);
    }

    population[0].cluster = initial;
    population[0].energy = evaluateCluster(population[0].cluster, lbfgsPool[0].get());

    int countStructureFiles = population[0].cluster.countStructureFiles();
    int minLoadAttempts = 10;

    for (int i = 1; i < populationSize; ++i) {
        bool loadSuccess = false;

        if (countStructureFiles > 0 && i <= minLoadAttempts) {
            int fileIndex = ((i - 1) % countStructureFiles) + 1;
            loadSuccess = population[i].cluster.loadStructureInitialize(fileIndex, numAtoms);
        }

        if (!loadSuccess)
            population[i].cluster.randomInitialize(2.75);

        population[i].energy = evaluateCluster(population[i].cluster, lbfgsPool[0].get());
    }

    bestIndividual = population[0];
    for (const auto& ind : population) {
        if (ind.energy < bestIndividual.energy) {
            bestIndividual = ind;
        }
    }
}

double CDE_Population::evaluateCluster(BinaryAlloyCluster& cluster, NELbfgs* lbfgs) {
    if (lbfgs) {
        lbfgs->optimize(cluster);
    }

    double energy = potential->calculateEnergy(cluster);
    cluster.setEnergy(energy);
    return energy;
}

void CDE_Population::evolve() {
    mutation();
    crossover();
    selection();
    if (numElementA > 0 && numElementA < numAtoms) {
        swapAtoms();
    }
}

void CDE_Population::mutation() {
    for (int i = 0; i < populationSize; ++i) {
        double F = RandomGenerator::uniform(0.1, 0.9);

        switch (strategy) {
        case CDE_RAND1: {
            auto indices = RandomGenerator::permutation(populationSize, 3);
            auto& coords = mutantPopulation[i].cluster.getCoordinates();
            const auto& coords0 = population[indices[0]].cluster.getCoordinates();
            const auto& coords1 = population[indices[1]].cluster.getCoordinates();
            const auto& coords2 = population[indices[2]].cluster.getCoordinates();

            for (size_t j = 0; j < coords.size(); ++j) {
                coords[j] = coords2[j] + F * (coords0[j] - coords1[j]);
            }
            break;
        }

        case CDE_BEST1: {
            auto indices = RandomGenerator::permutation(populationSize, 2);
            auto& coords = mutantPopulation[i].cluster.getCoordinates();
            const auto& bestCoords = bestIndividual.cluster.getCoordinates();
            const auto& coords0 = population[indices[0]].cluster.getCoordinates();
            const auto& coords1 = population[indices[1]].cluster.getCoordinates();

            for (size_t j = 0; j < coords.size(); ++j) {
                coords[j] = bestCoords[j] + F * (coords0[j] - coords1[j]);
            }
            break;
        }

        case CDE_RAND_TO_BEST1: {
            auto indices = RandomGenerator::permutation(populationSize, 2);
            auto& coords = mutantPopulation[i].cluster.getCoordinates();
            const auto& currentCoords = population[i].cluster.getCoordinates();
            const auto& bestCoords = bestIndividual.cluster.getCoordinates();
            const auto& coords0 = population[indices[0]].cluster.getCoordinates();
            const auto& coords1 = population[indices[1]].cluster.getCoordinates();

            for (size_t j = 0; j < coords.size(); ++j) {
                coords[j] = currentCoords[j] + F * (bestCoords[j] - currentCoords[j]) +
                    F * (coords0[j] - coords1[j]);
            }
            break;
        }
        }
    }
}

void CDE_Population::crossover() {
    if (threadPool) {
        // === Parallel path ===
        int poolSize = static_cast<int>(lbfgsPool.size());
        std::atomic<int> lbfgsIdx{0};

        struct EvalTask {
            BinaryAlloyCluster* cluster;
            double* energy;
        };

        struct SpliceInfo {
            int child2Idx;   // index into child2Pool
            int popIndex;    // index into trialPopulation
        };

        std::vector<EvalTask> evalTasks;
        std::vector<BinaryAlloyCluster> child2Pool;
        std::vector<double> child2Energies;
        std::vector<SpliceInfo> spliceInfos;
        evalTasks.reserve(populationSize * 2);
        child2Pool.reserve(populationSize);
        child2Energies.reserve(populationSize);
        spliceInfos.reserve(populationSize);

        for (int i = 0; i < populationSize; ++i) {
            double CR = RandomGenerator::uniform();

            if (CR < 0.3) {
                int jrand = RandomGenerator::uniformInt(0, 3 * numAtoms - 1);
                auto& trialCoords = trialPopulation[i].cluster.getCoordinates();
                const auto& currentCoords = population[i].cluster.getCoordinates();
                const auto& mutantCoords = mutantPopulation[i].cluster.getCoordinates();

                for (int j = 0; j < 3 * numAtoms; ++j) {
                    if (RandomGenerator::uniform() < CR || j == jrand) {
                        trialCoords[j] = mutantCoords[j];
                    }
                    else {
                        trialCoords[j] = currentCoords[j];
                    }
                }
                evalTasks.push_back({ &trialPopulation[i].cluster, &trialPopulation[i].energy });
            }
            else {
                child2Pool.emplace_back(numElementA, numElementB, elementA, elementB);
                BinaryAlloyCluster& child2 = child2Pool.back();
                sphereCutSplice(population[i].cluster, mutantPopulation[i].cluster,
                    trialPopulation[i].cluster, child2);

                child2Energies.push_back(0.0);
                spliceInfos.push_back({ static_cast<int>(child2Pool.size()) - 1, i });
                evalTasks.push_back({ &trialPopulation[i].cluster, &trialPopulation[i].energy });
                evalTasks.push_back({ &child2, &child2Energies.back() });
            }
        }

        if (!evalTasks.empty()) {
            for (auto& task : evalTasks) {
                threadPool->submit([this, &task, &lbfgsIdx, poolSize]() {
                    int idx = lbfgsIdx.fetch_add(1) % poolSize;
                    *task.energy = evaluateCluster(*task.cluster, lbfgsPool[idx].get());
                });
            }
            threadPool->waitAll();
        }

        for (auto& info : spliceInfos) {
            double child2Energy = child2Energies[info.child2Idx];
            if (child2Energy < trialPopulation[info.popIndex].energy) {
                trialPopulation[info.popIndex].energy = child2Energy;
                trialPopulation[info.popIndex].cluster = std::move(child2Pool[info.child2Idx]);
            }
        }
    }
    else {
        // === Serial fallback ===
        NELbfgs* lbfgs = lbfgsPool[0].get();
        for (int i = 0; i < populationSize; ++i) {
            double CR = RandomGenerator::uniform();

            if (CR < 0.3) {
                int jrand = RandomGenerator::uniformInt(0, 3 * numAtoms - 1);
                auto& trialCoords = trialPopulation[i].cluster.getCoordinates();
                const auto& currentCoords = population[i].cluster.getCoordinates();
                const auto& mutantCoords = mutantPopulation[i].cluster.getCoordinates();

                for (int j = 0; j < 3 * numAtoms; ++j) {
                    if (RandomGenerator::uniform() < CR || j == jrand) {
                        trialCoords[j] = mutantCoords[j];
                    }
                    else {
                        trialCoords[j] = currentCoords[j];
                    }
                }
                trialPopulation[i].energy = evaluateCluster(trialPopulation[i].cluster, lbfgs);
            }
            else {
                BinaryAlloyCluster child2(numElementA, numElementB, elementA, elementB);
                sphereCutSplice(population[i].cluster, mutantPopulation[i].cluster,
                    trialPopulation[i].cluster, child2);

                trialPopulation[i].energy = evaluateCluster(trialPopulation[i].cluster, lbfgs);
                double energy2 = evaluateCluster(child2, lbfgs);

                if (energy2 < trialPopulation[i].energy) {
                    trialPopulation[i].cluster = std::move(child2);
                    trialPopulation[i].energy = energy2;
                }
            }
        }
    }
}

void CDE_Population::selection() {
    for (int i = 0; i < populationSize; ++i) {
        if (trialPopulation[i].energy < population[i].energy) {
            population[i] = trialPopulation[i];

            if (population[i].energy < bestIndividual.energy) {
                bestIndividual = population[i];
            }
        }
    }
}

void CDE_Population::swapAtoms() {
    if (threadPool) {
        // === Parallel path ===
        int poolSize = static_cast<int>(lbfgsPool.size());
        std::atomic<int> lbfgsIdx{0};

        struct SwapJob {
            BinaryAlloyCluster newCluster;
            int populationIdx;
            double energy;
            bool valid;
        };

        std::vector<SwapJob> swapJobs;
        swapJobs.reserve(populationSize);

        for (int idx = 0; idx < populationSize; ++idx) {
            if (fabs(population[idx].energy - bestIndividual.energy) < 0.2 && RandomGenerator::uniform() < 0.9) {
                BinaryAlloyCluster newCluster = population[idx].cluster;

                std::vector<int> typeA_indices, typeB_indices;
                for (int i = 0; i < numAtoms; ++i) {
                    if (newCluster.getAtomType(i) == 0) {
                        typeA_indices.push_back(i);
                    }
                    else {
                        typeB_indices.push_back(i);
                    }
                }

                if (!typeA_indices.empty() && !typeB_indices.empty()) {
                    int idxA = typeA_indices[RandomGenerator::uniformInt(0, static_cast<int>(typeA_indices.size()) - 1)];
                    int idxB = typeB_indices[RandomGenerator::uniformInt(0, static_cast<int>(typeB_indices.size()) - 1)];

                    auto posA = newCluster.getAtomPosition(idxA);
                    auto posB = newCluster.getAtomPosition(idxB);
                    newCluster.setAtomPosition(idxA, posB[0], posB[1], posB[2]);
                    newCluster.setAtomPosition(idxB, posA[0], posA[1], posA[2]);

                    swapJobs.push_back({ newCluster, idx, 0.0, true });
                }
            }
        }

        if (!swapJobs.empty()) {
            for (auto& job : swapJobs) {
                threadPool->submit([this, &job, &lbfgsIdx, poolSize]() {
                    int lbfgsId = lbfgsIdx.fetch_add(1) % poolSize;
                    job.energy = evaluateCluster(job.newCluster, lbfgsPool[lbfgsId].get());
                });
            }
            threadPool->waitAll();
        }

        for (auto& job : swapJobs) {
            if (job.valid && job.energy < population[job.populationIdx].energy) {
                population[job.populationIdx].cluster = job.newCluster;
                population[job.populationIdx].energy = job.energy;

                if (job.energy < bestIndividual.energy) {
                    bestIndividual.cluster = job.newCluster;
                    bestIndividual.energy = job.energy;
                }
            }
        }
    }
    else {
        // === Serial fallback ===
        NELbfgs* lbfgs = lbfgsPool[0].get();
        for (int idx = 0; idx < populationSize; ++idx) {
            if (fabs(population[idx].energy - bestIndividual.energy) < 0.2 && RandomGenerator::uniform() < 0.9) {
                BinaryAlloyCluster newCluster = population[idx].cluster;

                std::vector<int> typeA_indices, typeB_indices;
                for (int i = 0; i < numAtoms; ++i) {
                    if (newCluster.getAtomType(i) == 0) {
                        typeA_indices.push_back(i);
                    }
                    else {
                        typeB_indices.push_back(i);
                    }
                }

                if (!typeA_indices.empty() && !typeB_indices.empty()) {
                    int idxA = typeA_indices[RandomGenerator::uniformInt(0, static_cast<int>(typeA_indices.size()) - 1)];
                    int idxB = typeB_indices[RandomGenerator::uniformInt(0, static_cast<int>(typeB_indices.size()) - 1)];

                    auto posA = newCluster.getAtomPosition(idxA);
                    auto posB = newCluster.getAtomPosition(idxB);
                    newCluster.setAtomPosition(idxA, posB[0], posB[1], posB[2]);
                    newCluster.setAtomPosition(idxB, posA[0], posA[1], posA[2]);

                    double newEnergy = evaluateCluster(newCluster, lbfgs);

                    if (newEnergy < population[idx].energy) {
                        population[idx].cluster = newCluster;
                        population[idx].energy = newEnergy;

                        if (newEnergy < bestIndividual.energy) {
                            bestIndividual.cluster = newCluster;
                            bestIndividual.energy = newEnergy;
                        }
                    }
                }
            }
        }
    }
}

void CDE_Population::sphereCutSplice(const BinaryAlloyCluster& parent1,
    const BinaryAlloyCluster& parent2,
    BinaryAlloyCluster& child1,
    BinaryAlloyCluster& child2) {
    int n = parent1.getNumAtoms();

    BinaryAlloyCluster p1 = parent1.copy();
    BinaryAlloyCluster p2 = parent2.copy();
    p1.centerAtOrigin();
    p2.centerAtOrigin();

    std::vector<double> dist1(n), dist2(n);
    for (int i = 0; i < n; ++i) {
        auto pos1 = p1.getAtomPosition(i);
        auto pos2 = p2.getAtomPosition(i);
        dist1[i] = std::sqrt(pos1[0] * pos1[0] + pos1[1] * pos1[1] + pos1[2] * pos1[2]);
        dist2[i] = std::sqrt(pos2[0] * pos2[0] + pos2[1] * pos2[1] + pos2[2] * pos2[2]);
    }

    auto sorted1 = dist1;
    auto sorted2 = dist2;
    std::sort(sorted1.begin(), sorted1.end());
    std::sort(sorted2.begin(), sorted2.end());

    std::vector<int> validCuts;
    for (int k = 1; k < n - 1; ++k) {
        if (sorted1[k - 1] < sorted2[k] && sorted2[k - 1] < sorted1[k]) {
            validCuts.push_back(k);
        }
    }

    if (!validCuts.empty()) {
        int cutIdx = validCuts[RandomGenerator::uniformInt(0, validCuts.size() - 1)];
        double cutRadius = (sorted1[cutIdx - 1] + sorted2[cutIdx]) / 2.0;

        int idx1 = 0, idx2 = 0;
        for (int i = 0; i < n && (idx1 < n || idx2 < n); ++i) {
            if (dist1[i] < cutRadius && idx1 < n) {
                auto pos = p1.getAtomPosition(i);
                child1.setAtomPosition(idx1++, pos[0], pos[1], pos[2]);
            }
            if (dist2[i] >= cutRadius && idx1 < n) {
                auto pos = p2.getAtomPosition(i);
                child1.setAtomPosition(idx1++, pos[0], pos[1], pos[2]);
            }

            if (dist1[i] >= cutRadius && idx2 < n) {
                auto pos = p1.getAtomPosition(i);
                child2.setAtomPosition(idx2++, pos[0], pos[1], pos[2]);
            }
            if (dist2[i] < cutRadius && idx2 < n) {
                auto pos = p2.getAtomPosition(i);
                child2.setAtomPosition(idx2++, pos[0], pos[1], pos[2]);
            }
        }
    }
    else {
        child1 = parent1;
        child2 = parent2;
    }
}


void CDE_Population::receiveIndividual(const CDE_Individual& ind) {
    int worstIdx = 0;
    double worstEnergy = population[0].energy;

    for (int i = 1; i < populationSize; ++i) {
        if (population[i].energy > worstEnergy) {
            worstEnergy = population[i].energy;
            worstIdx = i;
        }
    }

    if (ind.energy < worstEnergy) {
        population[worstIdx] = ind;

        if (ind.energy < bestIndividual.energy) {
            bestIndividual = ind;
        }
    }
}


// ==================== CDE Main Class Implementation ====================

CDE::CDE(const Parameters& p, PotentialBase* pot)
    : params(p), potential(pot),
    generation(0), evaluationCount(0),
    globalBest(1, 0, "A", "B") {

    if (params.useThreading) {
        int nThreads = std::max(1, params.numThreads);
        threadPool = std::make_unique<ThreadPool>(nThreads);
    }
}

void CDE::initialize(const BinaryAlloyCluster& initial) {
    populations.clear();

    ThreadPool* pool = threadPool.get();

    if (params.useMultiPopulation) {
        populations.emplace_back(std::make_unique<CDE_Population>(
            CDE_RAND1, params.populationSize, initial, potential, pool));

        populations.emplace_back(std::make_unique<CDE_Population>(
            CDE_BEST1, params.populationSize, initial, potential, pool));

        populations.emplace_back(std::make_unique<CDE_Population>(
            CDE_RAND_TO_BEST1, params.populationSize, initial, potential, pool));
    }
    else {
        populations.emplace_back(std::make_unique<CDE_Population>(
            CDE_RAND1, params.populationSize, initial, potential, pool));
    }

    globalBest = populations[0]->getBestIndividual();
    for (const auto& pop : populations) {
        if (pop->getBestIndividual().energy < globalBest.energy) {
            globalBest = pop->getBestIndividual();
        }
    }
}

void CDE::updateGlobalBest(const CDE_Individual& candidate) {
    std::lock_guard<std::mutex> lock(globalBestMutex);
    if (candidate.energy < globalBest.energy) {
        globalBest = candidate;
    }
}

void CDE::evolve() {
    generation++;

    if (params.useThreading && threadPool) {
        // Mutation and crossover for all populations
        for (auto& pop : populations) {
            pop->mutation();
            pop->crossover();
        }

        // Selection (fast, serial)
        for (auto& pop : populations) {
            pop->selection();
            updateGlobalBest(pop->getBestIndividual());
        }

        // Swap atoms (uses thread pool internally)
        for (auto& pop : populations) {
            pop->swapAtoms();
            updateGlobalBest(pop->getBestIndividual());
        }
    }
    else {
        for (auto& pop : populations) {
            pop->evolve();
        }

        for (const auto& pop : populations) {
            if (pop->getBestIndividual().energy < globalBest.energy) {
                globalBest = pop->getBestIndividual();
            }
        }
    }

    if (params.useMultiPopulation && generation % params.exchangeInterval == 0) {
        exchangeBestIndividuals();
    }
}

void CDE::exchangeBestIndividuals() {
    for (size_t i = 0; i < populations.size(); ++i) {
        const auto& bestInd = populations[i]->getBestIndividual();
        for (size_t j = 0; j < populations.size(); ++j) {
            if (i != j) {
                populations[j]->receiveIndividual(bestInd);
            }
        }
    }
}

BinaryAlloyCluster CDE::optimize(const BinaryAlloyCluster& initial,
    ResultManager* resultManager) {
    initialize(initial);

    for (int gen = 0; gen < params.maxGenerations; ++gen) {
        evolve();

        if (resultManager) {
            resultManager->saveGenerationBest(generation, globalBest.cluster);
            resultManager->updateHistoricalBest(globalBest.cluster);

            std::vector<BinaryAlloyCluster> dummyPop;
            for (const auto& pop : populations) {
                dummyPop.push_back(pop->getBestIndividual().cluster);
            }
            resultManager->logGeneration(generation, dummyPop, globalBest.cluster);
        }

        if (gen % 1 == 0) {
            std::cout << "Gen " << gen << ": Best = " << globalBest.energy;

            if (params.useMultiPopulation) {
                std::cout << " [";
                for (size_t i = 0; i < populations.size(); ++i) {
                    if (i > 0) std::cout << ", ";
                    std::cout << populations[i]->getBestIndividual().energy;
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
