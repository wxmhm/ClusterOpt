#include "../include/SaNSDE.h"

// ==================== SaNSDE_Population Implementation ====================

SaNSDE_Population::SaNSDE_Population(int id, const Parameters& p, PotentialBase* pot,
    const BinaryAlloyCluster& initial, ThreadPool* pool)
    : populationId(id), params(p), potential(pot), threadPool(pool),
    generation(0),
    elementA(initial.getElementA()), elementB(initial.getElementB()),
    numElementA(initial.getNumElementA()), numElementB(initial.getNumElementB()),
    bestCluster(initial.getNumElementA(), initial.getNumElementB(),
        initial.getElementA(), initial.getElementB()) {

    availableStrategies = {
        MutationStrategy::RAND1,
        MutationStrategy::RAND_TO_BEST1,
        MutationStrategy::RAND2,
        MutationStrategy::CURRENT_TO_BEST1
    };

    int numStrategies = static_cast<int>(availableStrategies.size());
    strategyProbabilities.resize(numStrategies, 1.0 / numStrategies);
    strategySuccesses.resize(numStrategies, 0);
    strategyFailures.resize(numStrategies, 0);

    population.reserve(params.populationSize);
    mutantPopulation.reserve(params.populationSize);
    trialPopulation.reserve(params.populationSize);

    F_values.resize(params.populationSize, 0.5);
    CR_values.resize(params.populationSize, 0.5);
    strategies.resize(params.populationSize, 0);

    // Create NELbfgs instances for parallel evaluation
    int numLbfgs = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
    for (int i = 0; i < numLbfgs; ++i) {
        lbfgsPool.push_back(std::make_unique<NELbfgs>(pot));
    }

    initializePopulation(initial);
}

void SaNSDE_Population::initializePopulation(const BinaryAlloyCluster& initial) {
    population.clear();
    mutantPopulation.clear();
    trialPopulation.clear();

    for (int i = 0; i < params.populationSize; ++i) {
        population.emplace_back(numElementA, numElementB, elementA, elementB);
        mutantPopulation.emplace_back(numElementA, numElementB, elementA, elementB);
        trialPopulation.emplace_back(numElementA, numElementB, elementA, elementB);
    }

    population[0] = initial;
    population[0].setEnergy(evaluateCluster(population[0], lbfgsPool[0].get()));

    int countStructureFiles = population[0].countStructureFiles();
    int minLoadAttempts = 10;

    for (int i = 1; i < params.populationSize; ++i) {
        bool loadSuccess = false;
        if (countStructureFiles > 0 && i <= minLoadAttempts) {
            int fileIndex = ((i - 1) % countStructureFiles) + 1;
            loadSuccess = population[i].loadStructureInitialize(fileIndex, numElementA+ numElementB);
        }

        if (!loadSuccess)
            population[i].randomInitialize(2.75);

        population[i].setEnergy(evaluateCluster(population[i], lbfgsPool[0].get()));

        F_values[i] = 0.5 + 0.3 * RandomGenerator::uniform();
        CR_values[i] = RandomGenerator::uniform();
        strategies[i] = RandomGenerator::uniformInt(0, static_cast<int>(availableStrategies.size()) - 1);
    }

    bestCluster = population[0];
    for (const auto& ind : population) {
        if (ind.getEnergy() < bestCluster.getEnergy()) {
            bestCluster = ind;
        }
    }
}

double SaNSDE_Population::evaluateCluster(BinaryAlloyCluster& cluster, NELbfgs* lbfgs) {
    int count = evaluationCount.fetch_add(1) + 1;

    if (params.useLocalSearch && count % params.localSearchFrequency == 0) {
        lbfgs->optimize(cluster);
        localSearchCount.fetch_add(1);
    }

    double energy = potential->calculateEnergy(cluster);
    cluster.setEnergy(energy);
    return energy;
}

void SaNSDE_Population::evolve(int currentGen) {
    generation = currentGen;

    if (generation % params.learningPeriod == 0) {
        updateStrategyProbabilities();
    }

    if (threadPool) {
        // === Parallel path ===
        generateTrials(currentGen);
        evaluateTrials();
        doNeighborhoodSearch(currentGen);
        applySelection();
    } else {
        // === Serial fallback ===
        NELbfgs* lbfgs = lbfgsPool[0].get();

        for (int i = 0; i < params.populationSize; ++i) {
            strategies[i] = selectStrategy();
            adaptParameters(i);
            mutation(i);

            // Crossover
            int n = population[0].getNumAtoms();
            if (RandomGenerator::uniform() < CR_values[i]) {
                int jrand = RandomGenerator::uniformInt(0, 3 * n - 1);
                auto& trialCoords = trialPopulation[i].getCoordinates();
                const auto& currentCoords = population[i].getCoordinates();
                const auto& mutantCoords = mutantPopulation[i].getCoordinates();

                for (int j = 0; j < 3 * n; ++j) {
                    if (RandomGenerator::uniform() < CR_values[i] || j == jrand) {
                        trialCoords[j] = mutantCoords[j];
                    } else {
                        trialCoords[j] = currentCoords[j];
                    }
                }
            } else {
                BinaryAlloyCluster child2(numElementA, numElementB, elementA, elementB);
                sphereCutSplice(population[i], mutantPopulation[i],
                    trialPopulation[i], child2);

                double e1 = evaluateCluster(trialPopulation[i], lbfgs);
                double e2 = evaluateCluster(child2, lbfgs);

                if (e2 < e1) {
                    trialPopulation[i] = child2;
                }
            }

            // Neighborhood search
            if (generation % 10 == 0 && RandomGenerator::uniform() < 0.2) {
                int ns = params.neighborhoodSizeMax -
                    (generation * (params.neighborhoodSizeMax - params.neighborhoodSizeMin)) /
                    params.maxGenerations;

                BinaryAlloyCluster best = population[i];
                double bestEnergy = best.getEnergy();

                for (int j = 0; j < ns; ++j) {
                    BinaryAlloyCluster neighbor = population[i];
                    auto& coords = neighbor.getCoordinates();

                    for (size_t k = 0; k < coords.size(); ++k) {
                        if (RandomGenerator::uniform() < 0.1) {
                            coords[k] += RandomGenerator::normal(0, 0.01);
                        }
                    }

                    double energy = evaluateCluster(neighbor, lbfgs);
                    if (energy < bestEnergy) {
                        best = neighbor;
                        bestEnergy = energy;
                    }
                }

                if (bestEnergy < population[i].getEnergy()) {
                    population[i] = best;
                }
            }
        }

        // Selection
        for (int i = 0; i < params.populationSize; ++i) {
            double trialEnergy = evaluateCluster(trialPopulation[i], lbfgs);

            if (trialEnergy < population[i].getEnergy()) {
                SuccessMemory sm;
                sm.F = F_values[i];
                sm.CR = CR_values[i];
                sm.strategy = strategies[i];
                sm.improvement = population[i].getEnergy() - trialEnergy;

                successMemory.push_back(sm);
                if (successMemory.size() > static_cast<size_t>(params.memorySize)) {
                    successMemory.pop_front();
                }

                strategySuccesses[strategies[i]]++;
                population[i] = trialPopulation[i];

                if (trialEnergy < bestCluster.getEnergy()) {
                    bestCluster = trialPopulation[i];
                }
            } else {
                strategyFailures[strategies[i]]++;
            }
        }
    }
}

// ==================== Parallel evolution phases ====================

void SaNSDE_Population::generateTrials(int currentGen) {
    child2Pool.clear();
    child2Energies.clear();
    spliceInfos.clear();
    usedSphereCutSplice.assign(params.populationSize, false);

    for (int i = 0; i < params.populationSize; ++i) {
        strategies[i] = selectStrategy();
        adaptParameters(i);
        mutation(i);

        int n = population[0].getNumAtoms();
        if (RandomGenerator::uniform() < CR_values[i]) {
            // Binomial crossover - generate trial without evaluation
            int jrand = RandomGenerator::uniformInt(0, 3 * n - 1);
            auto& trialCoords = trialPopulation[i].getCoordinates();
            const auto& currentCoords = population[i].getCoordinates();
            const auto& mutantCoords = mutantPopulation[i].getCoordinates();

            for (int j = 0; j < 3 * n; ++j) {
                if (RandomGenerator::uniform() < CR_values[i] || j == jrand) {
                    trialCoords[j] = mutantCoords[j];
                } else {
                    trialCoords[j] = currentCoords[j];
                }
            }
        } else {
            // SphereCutSplice - generate both children, evaluate later
            usedSphereCutSplice[i] = true;
            child2Pool.emplace_back(numElementA, numElementB, elementA, elementB);
            sphereCutSplice(population[i], mutantPopulation[i],
                trialPopulation[i], child2Pool.back());
            child2Energies.push_back(0.0);
            spliceInfos.push_back({static_cast<int>(child2Pool.size()) - 1, i});
        }
    }
}

void SaNSDE_Population::evaluateTrials() {
    int poolSize = static_cast<int>(lbfgsPool.size());
    std::atomic<int> lbfgsIdx{0};

    struct EvalTask {
        BinaryAlloyCluster* cluster;
        double* energy;
    };

    std::vector<EvalTask> evalTasks;
    evalTasks.reserve(params.populationSize + spliceInfos.size());

    // Collect all trial evaluations
    for (int i = 0; i < params.populationSize; ++i) {
        evalTasks.push_back({&trialPopulation[i], nullptr});
    }

    // Use trialPopulation[i].energy as storage for pre-computed energy
    // We'll store energies in a temp vector and assign later
    std::vector<double> trialEnergies(params.populationSize);

    for (int i = 0; i < params.populationSize; ++i) {
        evalTasks[i].energy = &trialEnergies[i];
    }

    // Collect child2 evaluations
    for (size_t s = 0; s < spliceInfos.size(); ++s) {
        evalTasks.push_back({&child2Pool[spliceInfos[s].child2Idx],
                             &child2Energies[spliceInfos[s].child2Idx]});
    }

    // Submit all to thread pool
    for (auto& task : evalTasks) {
        threadPool->submit([this, &task, &lbfgsIdx, poolSize]() {
            int idx = lbfgsIdx.fetch_add(1) % poolSize;
            *task.energy = evaluateCluster(*task.cluster, lbfgsPool[idx].get());
        });
    }
    threadPool->waitAll();

    // Assign trial energies
    for (int i = 0; i < params.populationSize; ++i) {
        trialPopulation[i].setEnergy(trialEnergies[i]);
    }

    // Resolve sphereCutSplice: pick better of trial vs child2
    for (auto& info : spliceInfos) {
        double child2Energy = child2Energies[info.child2Idx];
        if (child2Energy < trialPopulation[info.popIndex].getEnergy()) {
            trialPopulation[info.popIndex] = std::move(child2Pool[info.child2Idx]);
            trialPopulation[info.popIndex].setEnergy(child2Energy);
        }
    }
}

void SaNSDE_Population::doNeighborhoodSearch(int currentGen) {
    if (generation % 10 != 0) return;

    int poolSize = static_cast<int>(lbfgsPool.size());
    std::atomic<int> lbfgsIdx{0};

    struct NSData {
        int popIndex;
        int ns;
        std::vector<BinaryAlloyCluster> neighbors;
        std::vector<double> energies;
    };

    std::vector<NSData> allNS;
    allNS.reserve(params.populationSize);

    // Generate all neighbors (serial, fast)
    for (int i = 0; i < params.populationSize; ++i) {
        if (RandomGenerator::uniform() >= 0.2) continue;

        int ns = params.neighborhoodSizeMax -
            (generation * (params.neighborhoodSizeMax - params.neighborhoodSizeMin)) /
            params.maxGenerations;

        NSData data;
        data.popIndex = i;
        data.ns = ns;
        data.neighbors.reserve(ns);
        data.energies.resize(ns);

        for (int j = 0; j < ns; ++j) {
            data.neighbors.emplace_back(numElementA, numElementB, elementA, elementB);
            data.neighbors[j] = population[i];
            auto& coords = data.neighbors[j].getCoordinates();

            for (size_t k = 0; k < coords.size(); ++k) {
                if (RandomGenerator::uniform() < 0.1) {
                    coords[k] += RandomGenerator::normal(0, 0.01);
                }
            }
        }

        allNS.push_back(std::move(data));
    }

    if (allNS.empty()) return;

    // Submit all neighbor evaluations to thread pool
    for (auto& nsData : allNS) {
        for (int j = 0; j < nsData.ns; ++j) {
            threadPool->submit([this, &nsData, j, &lbfgsIdx, poolSize]() {
                int idx = lbfgsIdx.fetch_add(1) % poolSize;
                nsData.energies[j] = evaluateCluster(nsData.neighbors[j], lbfgsPool[idx].get());
            });
        }
    }
    threadPool->waitAll();

    // Apply best neighbor for each individual
    for (auto& nsData : allNS) {
        double bestEnergy = population[nsData.popIndex].getEnergy();
        int bestIdx = -1;

        for (int j = 0; j < nsData.ns; ++j) {
            if (nsData.energies[j] < bestEnergy) {
                bestEnergy = nsData.energies[j];
                bestIdx = j;
            }
        }

        if (bestIdx >= 0) {
            population[nsData.popIndex] = nsData.neighbors[bestIdx];
            population[nsData.popIndex].setEnergy(bestEnergy);
        }
    }
}

void SaNSDE_Population::applySelection() {
    for (int i = 0; i < params.populationSize; ++i) {
        double trialEnergy = trialPopulation[i].getEnergy();

        if (trialEnergy < population[i].getEnergy()) {
            SuccessMemory sm;
            sm.F = F_values[i];
            sm.CR = CR_values[i];
            sm.strategy = strategies[i];
            sm.improvement = population[i].getEnergy() - trialEnergy;

            successMemory.push_back(sm);
            if (successMemory.size() > static_cast<size_t>(params.memorySize)) {
                successMemory.pop_front();
            }

            strategySuccesses[strategies[i]]++;
            population[i] = trialPopulation[i];

            if (trialEnergy < bestCluster.getEnergy()) {
                bestCluster = trialPopulation[i];
            }
        } else {
            strategyFailures[strategies[i]]++;
        }
    }
}

// ==================== Per-individual operations ====================

void SaNSDE_Population::mutation(int idx) {
    int n = population[0].getNumAtoms();
    auto strategy = availableStrategies[strategies[idx]];

    switch (strategy) {
    case MutationStrategy::RAND1: {
        auto indices = RandomGenerator::permutation(params.populationSize, 3);
        auto& coords = mutantPopulation[idx].getCoordinates();
        for (int j = 0; j < 3 * n; ++j) {
            coords[j] = population[indices[2]].getCoordinates()[j] +
                F_values[idx] * (population[indices[0]].getCoordinates()[j] -
                    population[indices[1]].getCoordinates()[j]);
        }
        break;
    }

    case MutationStrategy::RAND_TO_BEST1: {
        auto indices = RandomGenerator::permutation(params.populationSize, 2);
        auto& coords = mutantPopulation[idx].getCoordinates();
        for (int j = 0; j < 3 * n; ++j) {
            coords[j] = population[idx].getCoordinates()[j] +
                F_values[idx] * (bestCluster.getCoordinates()[j] -
                    population[idx].getCoordinates()[j]) +
                F_values[idx] * (population[indices[0]].getCoordinates()[j] -
                    population[indices[1]].getCoordinates()[j]);
        }
        break;
    }

    case MutationStrategy::RAND2: {
        auto indices = RandomGenerator::permutation(params.populationSize, 5);
        auto& coords = mutantPopulation[idx].getCoordinates();
        for (int j = 0; j < 3 * n; ++j) {
            coords[j] = population[indices[4]].getCoordinates()[j] +
                F_values[idx] * (population[indices[0]].getCoordinates()[j] +
                    population[indices[1]].getCoordinates()[j] -
                    population[indices[2]].getCoordinates()[j] -
                    population[indices[3]].getCoordinates()[j]);
        }
        break;
    }

    case MutationStrategy::CURRENT_TO_BEST1: {
        auto indices = RandomGenerator::permutation(params.populationSize, 2);
        double K = 0.5 * (F_values[idx] + 1);
        auto& coords = mutantPopulation[idx].getCoordinates();
        for (int j = 0; j < 3 * n; ++j) {
            coords[j] = population[idx].getCoordinates()[j] +
                K * (bestCluster.getCoordinates()[j] -
                    population[idx].getCoordinates()[j]) +
                F_values[idx] * (population[indices[0]].getCoordinates()[j] -
                    population[indices[1]].getCoordinates()[j]);
        }
        break;
    }

    default:
        mutantPopulation[idx] = population[idx];
        break;
    }
}

void SaNSDE_Population::adaptParameters(int idx) {
    F_values[idx] = generateF(idx);
    CR_values[idx] = generateCR(idx);
}

double SaNSDE_Population::generateF(int idx) {
    double F;
    if (!successMemory.empty() && RandomGenerator::uniform() < 0.5) {
        int memIdx = RandomGenerator::uniformInt(0, static_cast<int>(successMemory.size()) - 1);
        F = RandomGenerator::cauchy(successMemory[memIdx].F, 0.1);
    } else {
        F = RandomGenerator::cauchy(0.5, 0.3);
    }

    return (std::max)(params.F_min, (std::min)(params.F_max, F));
}

double SaNSDE_Population::generateCR(int idx) {
    double CR;
    if (!successMemory.empty() && RandomGenerator::uniform() < 0.5) {
        int memIdx = RandomGenerator::uniformInt(0, static_cast<int>(successMemory.size()) - 1);
        CR = RandomGenerator::normal(successMemory[memIdx].CR, 0.1);
    } else {
        CR = RandomGenerator::normal(0.5, 0.1);
    }

    return (std::max)(params.CR_min, (std::min)(params.CR_max, CR));
}

int SaNSDE_Population::selectStrategy() {
    double r = RandomGenerator::uniform();
    double cumProb = 0;

    for (size_t i = 0; i < strategyProbabilities.size(); ++i) {
        cumProb += strategyProbabilities[i];
        if (r <= cumProb) {
            return static_cast<int>(i);
        }
    }

    return static_cast<int>(strategyProbabilities.size()) - 1;
}

void SaNSDE_Population::updateStrategyProbabilities() {
    int numStrategies = static_cast<int>(availableStrategies.size());
    double sumSuccess = 0;

    for (int i = 0; i < numStrategies; ++i) {
        double total = static_cast<double>(strategySuccesses[i] + strategyFailures[i]);
        if (total > 0) {
            strategyProbabilities[i] = static_cast<double>(strategySuccesses[i]) / total;
            sumSuccess += strategyProbabilities[i];
        } else {
            strategyProbabilities[i] = params.p_min;
        }

        strategySuccesses[i] = 0;
        strategyFailures[i] = 0;
    }

    if (sumSuccess > 0) {
        for (int i = 0; i < numStrategies; ++i) {
            strategyProbabilities[i] = strategyProbabilities[i] / sumSuccess *
                (1 - numStrategies * params.p_min) + params.p_min;
        }
    } else {
        std::fill(strategyProbabilities.begin(), strategyProbabilities.end(),
            1.0 / numStrategies);
    }
}

void SaNSDE_Population::sphereCutSplice(const BinaryAlloyCluster& parent1,
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
        int cutIdx = validCuts[RandomGenerator::uniformInt(0, static_cast<int>(validCuts.size()) - 1)];
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
    } else {
        child1 = parent1;
        child2 = parent2;
    }
}

void SaNSDE_Population::receiveIndividual(const BinaryAlloyCluster& cluster) {
    int worstIdx = 0;
    double worstEnergy = population[0].getEnergy();

    for (int i = 1; i < params.populationSize; ++i) {
        if (population[i].getEnergy() > worstEnergy) {
            worstEnergy = population[i].getEnergy();
            worstIdx = i;
        }
    }

    if (cluster.getEnergy() < worstEnergy) {
        population[worstIdx] = cluster;

        if (cluster.getEnergy() < bestCluster.getEnergy()) {
            bestCluster = cluster;
        }
    }
}

// ==================== SaNSDE Main Class Implementation ====================

SaNSDE::SaNSDE(const Parameters& p, PotentialBase* pot, bool useThreading)
    : params(p), potential(pot), generation(0), useThreading(useThreading),
    globalBest(1, 0, "A", "B") {

    if (useThreading) {
        int nThreads = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
        threadPool = std::make_unique<ThreadPool>(nThreads);
    }
}

void SaNSDE::initialize(const BinaryAlloyCluster& initial) {
    populations.clear();

    ThreadPool* pool = threadPool.get();

    for (int i = 0; i < NUM_POPULATIONS; ++i) {
        populations.emplace_back(
            std::make_unique<SaNSDE_Population>(i, params, potential, initial, pool)
        );
    }

    globalBest = populations[0]->getBestCluster();
    for (const auto& pop : populations) {
        if (pop->getBestEnergy() < globalBest.getEnergy()) {
            globalBest = pop->getBestCluster();
        }
    }

    generation = 0;
}

void SaNSDE::updateGlobalBest(const BinaryAlloyCluster& candidate) {
    std::lock_guard<std::mutex> lock(globalBestMutex);
    if (candidate.getEnergy() < globalBest.getEnergy()) {
        globalBest = candidate;
    }
}

void SaNSDE::evolve() {
    generation++;

    for (auto& pop : populations) {
        pop->evolve(generation);
        updateGlobalBest(pop->getBestCluster());
    }

    if (generation % 50 == 0) {
        migrationBetweenPopulations();
    }
}

void SaNSDE::migrationBetweenPopulations() {
    for (int i = 0; i < NUM_POPULATIONS; ++i) {
        int nextPop = (i + 1) % NUM_POPULATIONS;
        populations[nextPop]->receiveIndividual(populations[i]->getBestCluster());
    }
}

BinaryAlloyCluster SaNSDE::optimize(const BinaryAlloyCluster& initial, ResultManager* resultManager) {
    initialize(initial);

    for (int gen = 0; gen < params.maxGenerations; ++gen) {
        evolve();

        if (resultManager) {
            resultManager->saveGenerationBest(generation, globalBest);
            resultManager->updateHistoricalBest(globalBest);

            std::vector<BinaryAlloyCluster> dummyPop;
            for (const auto& pop : populations) {
                dummyPop.push_back(pop->getBestCluster());
            }
            resultManager->logGeneration(generation, dummyPop, globalBest);
        }

        if (gen % 1 == 0) {
            std::cout << "Gen " << gen << ": Best = " << globalBest.getEnergy();
            std::cout << " [";
            for (size_t i = 0; i < populations.size(); ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << populations[i]->getBestEnergy();
            }
            std::cout << "]";

            if (resultManager && resultManager->hasHistory()) {
                std::cout << ", Historical = " << resultManager->getHistoricalBestEnergy();
            }

            std::cout << std::endl;
        }
    }

    return globalBest;
}
