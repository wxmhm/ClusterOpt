#include "../include/PSO.h"

PSO_Population::PSO_Population(int popSize, const BinaryAlloyCluster& initial,
    GuptaPotential* pot, LocalOptimizer* opt, int maxGen)
    : populationSize(popSize), potential(pot), localOptimizer(opt),
    currentGeneration(0), maxGenerations(maxGen),
    globalBest(initial.getNumElementA(), initial.getNumElementB(),
        initial.getElementA(), initial.getElementB()) {

    numElementA = initial.getNumElementA();
    numElementB = initial.getNumElementB();
    elementA = initial.getElementA();
    elementB = initial.getElementB();
    numAtoms = numElementA + numElementB;

    omegaMax = 0.9;
    omegaMin = 0.4;
    omega = omegaMax;
    c1 = 2.0;
    c2 = 2.0;
    vmax = 0.1;
    psoRatio = 0.6;

    population.reserve(populationSize);
    for (int i = 0; i < populationSize; ++i) {
        population.emplace_back(numElementA, numElementB, elementA, elementB);
    }

    population[0].cluster = initial;
    population[0].energy = evaluateCluster(population[0].cluster);
    population[0].pbest = population[0].cluster;
    population[0].pbestEnergy = population[0].energy;

    int countStructureFiles = population[0].cluster.countStructureFiles();
    int minLoadAttempts = (std::min)(10, populationSize - 1);

    for (int i = 1; i < populationSize; ++i) {
        bool loadSuccess = false;
        if (countStructureFiles > 0 && i <= minLoadAttempts) {
            int fileIndex = ((i - 1) % countStructureFiles) + 1;
            loadSuccess = population[i].cluster.loadStructureInitialize(fileIndex, numAtoms);
        }
        if (!loadSuccess) {
            population[i].cluster.randomInitialize(2.75);
        }
        population[i].energy = evaluateCluster(population[i].cluster);
        population[i].pbest = population[i].cluster;
        population[i].pbestEnergy = population[i].energy;
    }

    initializeVelocities();

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
            ind.velocity[i] = RandomGenerator::uniform(-vmax, vmax);
        }
    }
}

void PSO_Population::updateOmega() {
    omega = omegaMax - (omegaMax - omegaMin) *
        static_cast<double>(currentGeneration) / maxGenerations;
}

void PSO_Population::clampVelocity(std::vector<double>& vel) {
    for (auto& v : vel) {
        if (v > vmax) v = vmax;
        if (v < -vmax) v = -vmax;
    }
}

double PSO_Population::evaluateCluster(BinaryAlloyCluster& cluster) {
    if (localOptimizer) {
        localOptimizer->optimize(cluster);
    }
    double energy = potential->calculateEnergy(cluster);
    cluster.setEnergy(energy);
    return energy;
}

void PSO_Population::evolve() {
    int numPSO = static_cast<int>(populationSize * psoRatio);
    int numRandom = populationSize - numPSO;

    for (int i = 0; i < numPSO; ++i) {
        updateVelocityAndPosition(i);
    }

    if (numRandom > 0) {
        generateRandomStructures(numPSO, numRandom);
    }

    updatePersonalAndGlobalBest();

    if (numElementA > 0 && numElementA < numAtoms) {
        swapAtoms();
    }

    currentGeneration++;
    updateOmega();
}

void PSO_Population::updateVelocityAndPosition(int idx) {
    auto& ind = population[idx];
    auto& coords = ind.cluster.getCoordinates();
    const auto& pbestCoords = ind.pbest.getCoordinates();
    const auto& gbestCoords = globalBest.cluster.getCoordinates();

    for (size_t j = 0; j < ind.velocity.size(); ++j) {
        double r1 = RandomGenerator::uniform();
        double r2 = RandomGenerator::uniform();
        ind.velocity[j] = omega * ind.velocity[j] +
            c1 * r1 * (pbestCoords[j] - coords[j]) +
            c2 * r2 * (gbestCoords[j] - coords[j]);
    }

    clampVelocity(ind.velocity);

    for (size_t j = 0; j < coords.size(); ++j) {
        coords[j] += ind.velocity[j];
    }

    ind.energy = evaluateCluster(ind.cluster);
}

void PSO_Population::generateRandomStructures(int startIdx, int count) {
    for (int i = 0; i < count; ++i) {
        int idx = startIdx + i;
        if (idx >= populationSize) break;
        population[idx].cluster.randomInitialize(2.75);
        population[idx].energy = evaluateCluster(population[idx].cluster);
        for (auto& v : population[idx].velocity) {
            v = RandomGenerator::uniform(-vmax, vmax);
        }
    }
}

void PSO_Population::updatePersonalAndGlobalBest() {
    for (auto& ind : population) {
        if (ind.energy < ind.pbestEnergy) {
            ind.pbest = ind.cluster;
            ind.pbestEnergy = ind.energy;
        }
        if (ind.energy < globalBest.energy) {
            globalBest.cluster = ind.cluster;
            globalBest.energy = ind.energy;
            globalBest.pbestEnergy = ind.energy;
        }
    }
}

void PSO_Population::swapAtoms() {
    for (int idx = 0; idx < populationSize; ++idx) {
        if (fabs(population[idx].energy - globalBest.energy) < 0.2 &&
            RandomGenerator::uniform() < 0.9) {
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
                int idxA = typeA_indices[RandomGenerator::uniformInt(0, typeA_indices.size() - 1)];
                int idxB = typeB_indices[RandomGenerator::uniformInt(0, typeB_indices.size() - 1)];
                auto posA = newCluster.getAtomPosition(idxA);
                auto posB = newCluster.getAtomPosition(idxB);
                newCluster.setAtomPosition(idxA, posB[0], posB[1], posB[2]);
                newCluster.setAtomPosition(idxB, posA[0], posA[1], posA[2]);
                double newEnergy = evaluateCluster(newCluster);
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

void PSO_Population::receiveIndividual(const PSO_Individual& ind) {
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
        if (ind.energy < globalBest.energy) {
            globalBest = ind;
        }
    }
}

PSO::PSO(GuptaPotential* pot, LocalOptimizer* opt)
    : potential(pot), localOptimizer(opt), generation(0),
    globalBest(1, 0, "A", "B") {

    populationSize = 30;
    maxGenerations = 200;
    exchangeInterval = 50;
    useMultiPopulation = true;
}

void PSO::initialize(const BinaryAlloyCluster& initial) {
    populations.clear();
    if (useMultiPopulation) {
        for (int i = 0; i < 3; ++i) {
            populations.emplace_back(std::make_unique<PSO_Population>(
                populationSize, initial, potential, localOptimizer, maxGenerations));
        }
    }
    else {
        populations.emplace_back(std::make_unique<PSO_Population>(
            populationSize, initial, potential, localOptimizer, maxGenerations));
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
    }
    for (const auto& pop : populations) {
        if (pop->getGlobalBest().energy < globalBest.energy) {
            globalBest = pop->getGlobalBest();
        }
    }
    if (useMultiPopulation && generation % exchangeInterval == 0) {
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
    for (int gen = 0; gen < maxGenerations; ++gen) {
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
            if (useMultiPopulation) {
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