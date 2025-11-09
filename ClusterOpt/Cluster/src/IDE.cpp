#include "../include/IDE.h"

// ==================== IDE_Population Implementation ====================

IDE_Population::IDE_Population(IDE_MutationStrategy strat, int popSize,
    const BinaryAlloyCluster& initial,
    GuptaPotential* pot, LocalOptimizer* opt)
    : strategy(strat), populationSize(popSize), potential(pot),
    localOptimizer(opt), localSearchCount(0),
    bestIndividual(initial.getNumElementA(), initial.getNumElementB(),
        initial.getElementA(), initial.getElementB()) {

    numElementA = initial.getNumElementA();
    numElementB = initial.getNumElementB();
    elementA = initial.getElementA();
    elementB = initial.getElementB();
    numAtoms = numElementA + numElementB;

    population.reserve(populationSize);
    mutantPopulation.reserve(populationSize);
    trialPopulation.reserve(populationSize);

    for (int i = 0; i < populationSize; ++i) {
        population.emplace_back(numElementA, numElementB, elementA, elementB);
        mutantPopulation.emplace_back(numElementA, numElementB, elementA, elementB);
        trialPopulation.emplace_back(numElementA, numElementB, elementA, elementB);
    }

    population[0].cluster = initial;
    population[0].energy = evaluateCluster(population[0].cluster);


    int countStructureFiles = population[0].cluster.countStructureFiles();
    int minLoadAttempts = 10;

    for (int i = 1; i < populationSize; ++i) {
        bool loadSuccess = false;

        if (countStructureFiles > 0 && i <= minLoadAttempts) {
            int fileIndex = ((i - 1) % countStructureFiles) + 1;
            loadSuccess = population[i].cluster.loadStructureInitialize(fileIndex, numAtoms);
            //std::cout << loadSuccess << "\t" << fileIndex << std::endl;
        }

        if (!loadSuccess)
            population[i].cluster.randomInitialize(2.75);

        population[i].energy = evaluateCluster(population[i].cluster);
    }

    bestIndividual = population[0];
    for (const auto& ind : population) {
        if (ind.energy < bestIndividual.energy) {
            bestIndividual = ind;
        }
    }
}

double IDE_Population::evaluateCluster(BinaryAlloyCluster& cluster) {
    if (localOptimizer) {
        localOptimizer->optimize(cluster);
    }

    double energy = potential->calculateEnergy(cluster);
    cluster.setEnergy(energy);
    return energy;
}

void IDE_Population::evolve() {
    mutation();
    crossover();
    selection();
    if (numElementA > 0 && numElementA < numAtoms) {
        swapAtoms();
    }
}

void IDE_Population::mutation() {
    for (int i = 0; i < populationSize; ++i) {
        double F = RandomGenerator::uniform(0.1, 0.9);

        switch (strategy) {
        case IDE_RAND1: {
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

        case IDE_BEST1: {
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

        case IDE_RAND2: {
            auto indices = RandomGenerator::permutation(populationSize, 5);
            auto& coords = mutantPopulation[i].cluster.getCoordinates();
            const auto& coords0 = population[indices[0]].cluster.getCoordinates();
            const auto& coords1 = population[indices[1]].cluster.getCoordinates();
            const auto& coords2 = population[indices[2]].cluster.getCoordinates();
            const auto& coords3 = population[indices[3]].cluster.getCoordinates();
            const auto& coords4 = population[indices[4]].cluster.getCoordinates();

            for (size_t j = 0; j < coords.size(); ++j) {
                coords[j] = coords4[j] + F * (coords0[j] + coords1[j] - coords2[j] - coords3[j]);
            }
            break;
        }

        case IDE_BEST2: {
            auto indices = RandomGenerator::permutation(populationSize, 4);
            auto& coords = mutantPopulation[i].cluster.getCoordinates();
            const auto& bestCoords = bestIndividual.cluster.getCoordinates();
            const auto& coords0 = population[indices[0]].cluster.getCoordinates();
            const auto& coords1 = population[indices[1]].cluster.getCoordinates();
            const auto& coords2 = population[indices[2]].cluster.getCoordinates();
            const auto& coords3 = population[indices[3]].cluster.getCoordinates();

            for (size_t j = 0; j < coords.size(); ++j) {
                coords[j] = bestCoords[j] + F * (coords0[j] + coords1[j] - coords2[j] - coords3[j]);
            }
            break;
        }

        case IDE_RAND_TO_BEST1: {
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

void IDE_Population::crossover() {
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
            trialPopulation[i].energy = evaluateCluster(trialPopulation[i].cluster);
        }
        else {
            BinaryAlloyCluster child2(numElementA, numElementB, elementA, elementB);
            sphereCutSplice(population[i].cluster, mutantPopulation[i].cluster,
                trialPopulation[i].cluster, child2);

            trialPopulation[i].energy = evaluateCluster(trialPopulation[i].cluster);
            double energy2 = evaluateCluster(child2);

            if (energy2 < trialPopulation[i].energy) {
                trialPopulation[i].cluster = child2;
                trialPopulation[i].energy = energy2;
            }
        }
    }
}

void IDE_Population::selection() {
    for (int i = 0; i < populationSize; ++i) {
        if (trialPopulation[i].energy < population[i].energy) {
            population[i] = trialPopulation[i];

            if (population[i].energy < bestIndividual.energy) {
                bestIndividual = population[i];
            }
        }
    }
}

void IDE_Population::swapAtoms() {
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

                double newEnergy = evaluateCluster(newCluster);

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

void IDE_Population::sphereCutSplice(const BinaryAlloyCluster& parent1,
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


void IDE_Population::receiveIndividual(const IDE_Individual& ind) {
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


// ==================== IDE Main Class Implementation ====================

IDE::IDE(const Parameters& p, GuptaPotential* pot, LocalOptimizer* opt)
    : params(p), potential(pot), localOptimizer(opt),
    generation(0), evaluationCount(0),
    globalBest(1, 0, "A", "B") {
}

void IDE::initialize(const BinaryAlloyCluster& initial) {
    populations.clear();

    if (params.useMultiPopulation) {
        populations.emplace_back(std::make_unique<IDE_Population>(
            IDE_RAND1, params.populationSize, initial, potential, localOptimizer));

        populations.emplace_back(std::make_unique<IDE_Population>(
            IDE_BEST1, params.populationSize, initial, potential, localOptimizer));

        populations.emplace_back(std::make_unique<IDE_Population>(
            IDE_RAND_TO_BEST1, params.populationSize, initial, potential, localOptimizer));
    }
    else {
        populations.emplace_back(std::make_unique<IDE_Population>(
            IDE_RAND1, params.populationSize, initial, potential, localOptimizer));
    }

    globalBest = populations[0]->getBestIndividual();
    for (const auto& pop : populations) {
        if (pop->getBestIndividual().energy < globalBest.energy) {
            globalBest = pop->getBestIndividual();
        }
    }
}

void IDE::updateGlobalBest(const IDE_Individual& candidate) {
    std::lock_guard<std::mutex> lock(globalBestMutex);
    if (candidate.energy < globalBest.energy) {
        globalBest = candidate;
    }
}

void IDE::evolve() {
    generation++;

    if (params.useThreading && params.useMultiPopulation) {
        std::vector<std::thread> threads;

        for (auto& pop : populations) {
            threads.emplace_back([this, &pop]() {
                pop->evolve();
                updateGlobalBest(pop->getBestIndividual());
                });
        }

        for (auto& thread : threads) {
            thread.join();
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

void IDE::exchangeBestIndividuals() {
    for (size_t i = 0; i < populations.size(); ++i) {
        const auto& bestInd = populations[i]->getBestIndividual();
        for (size_t j = 0; j < populations.size(); ++j) {
            if (i != j) {
                populations[j]->receiveIndividual(bestInd);
            }
        }
    }
}

BinaryAlloyCluster IDE::optimize(const BinaryAlloyCluster& initial,
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
