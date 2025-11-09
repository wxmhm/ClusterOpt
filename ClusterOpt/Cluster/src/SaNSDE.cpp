#include "../include/SaNSDE.h"

// ==================== SaNSDE_Population Implementation ====================

SaNSDE_Population::SaNSDE_Population(int id, const Parameters& p, GuptaPotential* pot,
    LocalOptimizer* opt, const BinaryAlloyCluster& initial)
    : populationId(id), params(p), potential(pot), localOptimizer(opt),
    generation(0), evaluationCount(0), localSearchCount(0),
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
    population[0].setEnergy(evaluateCluster(population[0]));

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

        population[i].setEnergy(evaluateCluster(population[i]));

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

double SaNSDE_Population::evaluateCluster(BinaryAlloyCluster& cluster) {
    evaluationCount++;

    if (params.useLocalSearch && evaluationCount % params.localSearchFrequency == 0) {
        localOptimizer->optimize(cluster);
        localSearchCount++;
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

    for (int i = 0; i < params.populationSize; ++i) {
        strategies[i] = selectStrategy();
        adaptParameters(i);
        mutation(i);
        crossover(i);

        if (generation % 10 == 0 && RandomGenerator::uniform() < 0.2) {
            neighborhoodSearch(i, currentGen);
        }
    }

    selection();
}

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

void SaNSDE_Population::crossover(int idx) {
    int n = population[0].getNumAtoms();

    if (RandomGenerator::uniform() < CR_values[idx]) {
        int jrand = RandomGenerator::uniformInt(0, 3 * n - 1);
        auto& trialCoords = trialPopulation[idx].getCoordinates();
        const auto& currentCoords = population[idx].getCoordinates();
        const auto& mutantCoords = mutantPopulation[idx].getCoordinates();

        for (int j = 0; j < 3 * n; ++j) {
            if (RandomGenerator::uniform() < CR_values[idx] || j == jrand) {
                trialCoords[j] = mutantCoords[j];
            }
            else {
                trialCoords[j] = currentCoords[j];
            }
        }
    }
    else {
        BinaryAlloyCluster child2(population[0].getNumElementA(),
            population[0].getNumElementB(),
            population[0].getElementA(),
            population[0].getElementB());
        sphereCutSplice(population[idx], mutantPopulation[idx],
            trialPopulation[idx], child2);

        double energy1 = evaluateCluster(trialPopulation[idx]);
        double energy2 = evaluateCluster(child2);

        if (energy2 < energy1) {
            trialPopulation[idx] = child2;
        }
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
        for (int i = 0; i < n; ++i) {
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

void SaNSDE_Population::selection() {
    for (int i = 0; i < params.populationSize; ++i) {
        double trialEnergy = evaluateCluster(trialPopulation[i]);

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
        }
        else {
            strategyFailures[strategies[i]]++;
        }
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
    }
    else {
        F = RandomGenerator::cauchy(0.5, 0.3);
    }

    return (std::max)(params.F_min, (std::min)(params.F_max, F));
}

double SaNSDE_Population::generateCR(int idx) {
    double CR;
    if (!successMemory.empty() && RandomGenerator::uniform() < 0.5) {
        int memIdx = RandomGenerator::uniformInt(0, static_cast<int>(successMemory.size()) - 1);
        CR = RandomGenerator::normal(successMemory[memIdx].CR, 0.1);
    }
    else {
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
        }
        else {
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
    }
    else {
        std::fill(strategyProbabilities.begin(), strategyProbabilities.end(),
            1.0 / numStrategies);
    }
}

void SaNSDE_Population::neighborhoodSearch(int idx, int currentGen) {
    int ns = params.neighborhoodSizeMax -
        (currentGen * (params.neighborhoodSizeMax - params.neighborhoodSizeMin)) /
        params.maxGenerations;

    BinaryAlloyCluster best = population[idx];
    double bestEnergy = best.getEnergy();

    for (int i = 0; i < ns; ++i) {
        BinaryAlloyCluster neighbor = population[idx];
        auto& coords = neighbor.getCoordinates();

        for (size_t j = 0; j < coords.size(); ++j) {
            if (RandomGenerator::uniform() < 0.1) {
                coords[j] += RandomGenerator::normal(0, 0.01);
            }
        }

        double energy = evaluateCluster(neighbor);
        if (energy < bestEnergy) {
            best = neighbor;
            bestEnergy = energy;
        }
    }

    if (bestEnergy < population[idx].getEnergy()) {
        population[idx] = best;
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

SaNSDE::SaNSDE(const Parameters& p, GuptaPotential* pot, LocalOptimizer* opt, bool useThreading)
    : params(p), potential(pot), localOptimizer(opt), generation(0), useThreading(useThreading),
    globalBest(1, 0, "A", "B") {
}

void SaNSDE::initialize(const BinaryAlloyCluster& initial) {
    populations.clear();

    for (int i = 0; i < NUM_POPULATIONS; ++i) {
        populations.emplace_back(
            std::make_unique<SaNSDE_Population>(i, params, potential, localOptimizer, initial)
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

    if (useThreading) {
        std::vector<std::thread> threads;

        for (auto& pop : populations) {
            threads.emplace_back([this, &pop]() {
                pop->evolve(generation);
                updateGlobalBest(pop->getBestCluster());
                });
        }

        for (auto& thread : threads) {
            thread.join();
        }
    }
    else {
        for (auto& pop : populations) {
            pop->evolve(generation);
            updateGlobalBest(pop->getBestCluster());
        }
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

int SaNSDE::getEvaluationCount() const {
    int total = 0;
    for (const auto& pop : populations) {
        // 注：需要在Population类中添加getter
    }
    return total;
}

int SaNSDE::getLocalSearchCount() const {
    int total = 0;
    for (const auto& pop : populations) {
        // 注：需要在Population类中添加getter
    }
    return total;
}

double SaNSDE::getDiversity() const {
    if (populations.empty()) return 0.0;

    double minEnergy = globalBest.getEnergy();
    double maxEnergy = globalBest.getEnergy();

    for (const auto& pop : populations) {
        minEnergy = (std::min)(minEnergy, pop->getBestEnergy());
        maxEnergy = (std::max)(maxEnergy, pop->getBestEnergy());
    }

    if (std::abs(maxEnergy - minEnergy) < Constants::EPSILON) {
        return 1.0;
    }

    return (maxEnergy - minEnergy) / maxEnergy;
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
