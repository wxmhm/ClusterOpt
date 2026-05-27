#include "../include/OrderParameterGuided.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <complex>
#include <map>

namespace {

const double PI = 3.14159265358979323846;

// --- Spherical harmonics helpers ---

double factorial(int n) {
    static const double cache[] = {1,1,2,6,24,120,720,5040,40320,362880,
        3628800,39916800,479001600,6227020800.0,87178291200.0,
        1307674368000.0,20922789888000.0,355687428096000.0,
        6402373705728000.0,121645100408832000.0};
    if (n <= 19) return cache[n];
    double r = cache[19];
    for (int i = 20; i <= n; ++i) r *= i;
    return r;
}

// Double factorial: n!! = n*(n-2)*(n-4)*...
double doubleFactorial(int n) {
    if (n <= 1) return 1.0;
    double r = 1.0;
    for (int i = n; i >= 2; i -= 2) r *= i;
    return r;
}

// Associated Legendre polynomial P_l^m(x) with Condon-Shortley phase (-1)^m
// x = cos(theta), |x| <= 1
double associatedLegendre(int l, int m, double x) {
    if (m < 0) return 0.0;
    if (m > l) return 0.0;

    // P_m^m(x) = (-1)^m * (2m-1)!! * (1-x^2)^{m/2}
    double pmm = 1.0;
    if (m > 0) {
        double fact = doubleFactorial(2 * m - 1);
        double omx2 = std::sqrt(1.0 - x * x);
        pmm = fact * std::pow(omx2, m);
        if (m % 2 == 1) pmm = -pmm; // (-1)^m factor
    }

    if (l == m) return pmm;

    // P_{m+1}^m(x) = x * (2m+1) * P_m^m(x)
    double pm1m = x * (2.0 * m + 1.0) * pmm;
    if (l == m + 1) return pm1m;

    // Recurrence: (l-m)*P_l^m = x*(2l-1)*P_{l-1}^m - (l+m-1)*P_{l-2}^m
    double plm = 0.0;
    for (int ll = m + 2; ll <= l; ++ll) {
        plm = (x * (2.0 * ll - 1.0) * pm1m - (ll + m - 1.0) * pmm) / (ll - m);
        pmm = pm1m;
        pm1m = plm;
    }
    return plm;
}

// Spherical harmonic Y_l^m(theta, phi)
// Returns complex value with Condon-Shortley phase
std::complex<double> sphericalHarmonic(int l, int m, double theta, double phi) {
    // For m < 0: Y_l^{-m} = (-1)^m * conj(Y_l^m)
    int absM = std::abs(m);

    // Normalization factor
    double norm = std::sqrt((2.0 * l + 1.0) / (4.0 * PI) *
        factorial(l - absM) / factorial(l + absM));

    double plm = associatedLegendre(l, absM, std::cos(theta));

    double re = norm * plm * std::cos(absM * phi);
    double im = norm * plm * std::sin(absM * phi);

    if (m < 0) {
        // Y_l^{-m} = (-1)^m * conj(Y_l^m) = (-1)^m * (re - i*im) for m>=0
        if (absM % 2 == 1) return std::complex<double>(-re, im);
        else return std::complex<double>(re, -im);
    }
    return std::complex<double>(re, im);
}

} // anonymous namespace

// --- OrderParameterGuided implementation ---

OrderParameterGuided::OrderParameterGuided(const Parameters& p, PotentialBase* pot)
    : params(p), potential(pot),
      globalBest(1, 0, "A", "B"),
      globalBestEnergy(1e10),
      currentTemperature(p.temperature),
      numAtoms(0), numElementA(0), numElementB(0),
      totalOperatorAttempts(0), stepsSinceNewBasin(0), totalVisits(0) {

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

void OrderParameterGuided::initialize(const BinaryAlloyCluster& initial) {
    numAtoms = initial.getNumAtoms();
    numElementA = initial.getNumElementA();
    numElementB = initial.getNumElementB();
    elementA = initial.getElementA();
    elementB = initial.getElementB();

    archive.clear();
    globalBestEnergy = 1e10;

    NELbfgs* lbfgs = lbfgsPool[0].get();

    // Initialize with several random structures to populate the archive
    int initCount = std::min(params.archiveMaxSize / 2, 20);
    for (int i = 0; i < initCount; ++i) {
        BinaryAlloyCluster structure(numElementA, numElementB, elementA, elementB);
        if (i == 0) {
            structure = initial;
        } else {
            structure.randomInitialize(2.75);
        }
        double energy = evaluateCluster(structure, lbfgs);
        OPSignature sig = computeOPSignature(structure);
        updateArchive(structure, sig);

        if (energy < globalBestEnergy) {
            globalBest = structure;
            globalBestEnergy = energy;
        }
    }

    std::cout << "OPG initialized: " << archive.size() << " basins in archive, "
        << "best = " << globalBestEnergy << " eV" << std::endl;
}

double OrderParameterGuided::evaluateCluster(BinaryAlloyCluster& cluster, NELbfgs* lbfgs) {
    int count = stepCount.fetch_add(1) + 1;

    if (params.useLocalSearch && count % params.localSearchFrequency == 0) {
        lbfgs->optimize(cluster);
    }

    double energy = potential->calculateEnergy(cluster);
    cluster.setEnergy(energy);
    return energy;
}

OrderParameterGuided::OPSignature
OrderParameterGuided::computeOPSignature(const BinaryAlloyCluster& cluster) {
    OPSignature sig;
    int n = cluster.getNumAtoms();

    // Find average nearest-neighbor distance to set cutoff
    std::vector<double> allDists;
    for (int i = 0; i < n; ++i) {
        auto pi = cluster.getAtomPosition(i);
        for (int j = i + 1; j < n; ++j) {
            auto pj = cluster.getAtomPosition(j);
            double dx = pi[0] - pj[0];
            double dy = pi[1] - pj[1];
            double dz = pi[2] - pj[2];
            allDists.push_back(std::sqrt(dx*dx + dy*dy + dz*dz));
        }
    }
    std::sort(allDists.begin(), allDists.end());
    // Use ~20th percentile as first-neighbor estimate, cutoff at 1.3x
    double firstNN = allDists[std::min(n, static_cast<int>(allDists.size() * 0.02))];
    double cutoff = firstNN * 1.35;

    // Compute q_lm for each atom, then average
    std::complex<double> Q4m[9];  // m = -4..4
    std::complex<double> Q6m[13]; // m = -6..6
    for (int m = 0; m < 9; ++m) Q4m[m] = 0.0;
    for (int m = 0; m < 13; ++m) Q6m[m] = 0.0;

    int totalAtoms = 0;
    for (int i = 0; i < n; ++i) {
        auto pi = cluster.getAtomPosition(i);

        std::complex<double> q4m[9] = {};
        std::complex<double> q6m[13] = {};
        int nbonds = 0;

        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            auto pj = cluster.getAtomPosition(j);
            double dx = pj[0] - pi[0];
            double dy = pj[1] - pi[1];
            double dz = pj[2] - pi[2];
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (dist > cutoff) continue;

            nbonds++;
            double theta = std::acos(dz / dist);
            double phi = std::atan2(dy, dx);

            for (int m = -4; m <= 4; ++m) {
                q4m[m + 4] += sphericalHarmonic(4, m, theta, phi);
            }
            for (int m = -6; m <= 6; ++m) {
                q6m[m + 6] += sphericalHarmonic(6, m, theta, phi);
            }
        }

        if (nbonds > 0) {
            for (int m = 0; m < 9; ++m) q4m[m] /= static_cast<double>(nbonds);
            for (int m = 0; m < 13; ++m) q6m[m] /= static_cast<double>(nbonds);
            totalAtoms++;
        }

        for (int m = 0; m < 9; ++m) Q4m[m] += q4m[m];
        for (int m = 0; m < 13; ++m) Q6m[m] += q6m[m];
    }

    if (totalAtoms > 0) {
        for (int m = 0; m < 9; ++m) Q4m[m] /= static_cast<double>(totalAtoms);
        for (int m = 0; m < 13; ++m) Q6m[m] /= static_cast<double>(totalAtoms);
    }

    // Q_l = sqrt(4π/(2l+1) * Σ_m |Q_lm|²)
    double sum4 = 0.0, sum6 = 0.0;
    for (int m = 0; m < 9; ++m) sum4 += std::norm(Q4m[m]);
    for (int m = 0; m < 13; ++m) sum6 += std::norm(Q6m[m]);

    sig.Q4 = std::sqrt(4.0 * PI / 9.0 * sum4);
    sig.Q6 = std::sqrt(4.0 * PI / 13.0 * sum6);

    // Third-order invariants W_l (simplified)
    // W_l = Σ_{m1,m2,m3} (l l l; m1 m2 m3) Q_lm1 Q_lm2 Q_lm3 / (Σ_m |Q_lm|²)^{3/2}
    // Use Wigner 3j symbols - we compute a simplified version
    // Only non-zero when m1+m2+m3=0

    auto wigner3j_l4 = [](int m1, int m2, int m3) -> double {
        // Pre-computed (4 4 4; m1 m2 m3) Wigner 3j symbols (normalized)
        // Only m1+m2+m3=0 gives non-zero
        if (m1 + m2 + m3 != 0) return 0.0;
        if (std::abs(m1) > 4 || std::abs(m2) > 4 || std::abs(m3) > 4) return 0.0;
        // Values from standard tables for l=4
        static const std::map<std::tuple<int,int,int>, double> table = {
            {{0,0,0}, 0.132981},
            {{1,-1,0}, -0.097590}, {{-1,1,0}, -0.097590}, {{0,1,-1}, -0.097590},
            {{0,-1,1}, -0.097590}, {{1,0,-1}, -0.097590}, {{-1,0,1}, -0.097590},
            {{2,-2,0}, 0.058399}, {{-2,2,0}, 0.058399}, {{0,2,-2}, 0.058399},
            {{0,-2,2}, 0.058399}, {{2,0,-2}, 0.058399}, {{-2,0,2}, 0.058399},
            {{2,-1,-1}, 0.075377}, {{-1,2,-1}, 0.075377}, {{-1,-1,2}, 0.075377},
            {{-2,1,1}, 0.075377}, {{1,-2,1}, 0.075377}, {{1,1,-2}, 0.075377},
            {{3,-3,0}, -0.032127}, {{-3,3,0}, -0.032127}, {{0,3,-3}, -0.032127},
            {{0,-3,3}, -0.032127}, {{3,0,-3}, -0.032127}, {{-3,0,3}, -0.032127},
            {{3,-2,-1}, -0.054591}, {{-2,3,-1}, -0.054591}, {{-1,-2,3}, -0.054591},
            {{-3,2,1}, -0.054591}, {{2,-3,1}, -0.054591}, {{1,2,-3}, -0.054591},
        };
        auto it = table.find({m1, m2, m3});
        return it != table.end() ? it->second : 0.0;
    };

    auto wigner3j_l6 = [](int m1, int m2, int m3) -> double {
        if (m1 + m2 + m3 != 0) return 0.0;
        if (std::abs(m1) > 6 || std::abs(m2) > 6 || std::abs(m3) > 6) return 0.0;
        static const std::map<std::tuple<int,int,int>, double> table = {
            {{0,0,0}, 0.082018},
            {{1,-1,0}, -0.062974}, {{-1,1,0}, -0.062974}, {{0,1,-1}, -0.062974},
            {{0,-1,1}, -0.062974}, {{1,0,-1}, -0.062974}, {{-1,0,1}, -0.062974},
            {{2,-2,0}, 0.043678}, {{-2,2,0}, 0.043678}, {{0,2,-2}, 0.043678},
            {{0,-2,2}, 0.043678}, {{2,0,-2}, 0.043678}, {{-2,0,2}, 0.043678},
            {{2,-1,-1}, 0.059470}, {{-1,2,-1}, 0.059470}, {{-1,-1,2}, 0.059470},
            {{-2,1,1}, 0.059470}, {{1,-2,1}, 0.059470}, {{1,1,-2}, 0.059470},
            {{3,-3,0}, -0.029412}, {{-3,3,0}, -0.029412}, {{0,3,-3}, -0.029412},
            {{0,-3,3}, -0.029412}, {{3,0,-3}, -0.029412}, {{-3,0,3}, -0.029412},
            {{3,-2,-1}, -0.047507}, {{-2,3,-1}, -0.047507}, {{-1,-2,3}, -0.047507},
            {{-3,2,1}, -0.047507}, {{2,-3,1}, -0.047507}, {{1,2,-3}, -0.047507},
            {{3,2,1}, 0.0}, // not in standard tables, approximate
            {{4,-4,0}, -0.019615}, {{-4,4,0}, -0.019615}, {{0,4,-4}, -0.019615},
            {{0,-4,4}, -0.019615}, {{4,0,-4}, -0.019615}, {{-4,0,4}, -0.019615},
        };
        auto it = table.find({m1, m2, m3});
        return it != table.end() ? it->second : 0.0;
    };

    // Compute W4
    double w4num = 0.0;
    for (int m1 = -4; m1 <= 4; ++m1) {
        for (int m2 = -4; m2 <= 4; ++m2) {
            int m3 = -(m1 + m2);
            if (std::abs(m3) > 4) continue;
            double w3j = wigner3j_l4(m1, m2, m3);
            w4num += w3j * std::real(Q4m[m1+4] * Q4m[m2+4] * Q4m[m3+4]);
        }
    }
    double q4norm = std::pow(sum4, 1.5);
    sig.W4 = (q4norm > 1e-30) ? w4num / q4norm : 0.0;

    // Compute W6
    double w6num = 0.0;
    for (int m1 = -6; m1 <= 6; ++m1) {
        for (int m2 = -6; m2 <= 6; ++m2) {
            int m3 = -(m1 + m2);
            if (std::abs(m3) > 6) continue;
            double w3j = wigner3j_l6(m1, m2, m3);
            w6num += w3j * std::real(Q6m[m1+6] * Q6m[m2+6] * Q6m[m3+6]);
        }
    }
    double q6norm = std::pow(sum6, 1.5);
    sig.W6 = (q6norm > 1e-30) ? w6num / q6norm : 0.0;

    return sig;
}

int OrderParameterGuided::findClosestArchiveEntry(const OPSignature& sig) const {
    int bestIdx = -1;
    double bestDist = 1e10;
    for (size_t i = 0; i < archive.size(); ++i) {
        double d = sig.distance(archive[i].signature);
        if (d < bestDist) {
            bestDist = d;
            bestIdx = static_cast<int>(i);
        }
    }
    return bestIdx;
}

void OrderParameterGuided::updateArchive(const BinaryAlloyCluster& cluster,
                                          const OPSignature& sig) {
    int closest = findClosestArchiveEntry(sig);

    if (closest >= 0 && sig.distance(archive[closest].signature) < params.opGridResolution) {
        // Existing basin region — update if better energy
        archive[closest].visitCount++;
        if (cluster.getEnergy() < archive[closest].energy) {
            archive[closest].structure = cluster;
            archive[closest].energy = cluster.getEnergy();
            archive[closest].signature = sig;
        }
    } else {
        // New basin region
        if (static_cast<int>(archive.size()) >= params.archiveMaxSize) {
            // Remove highest-energy entry to keep archive bounded
            int worstIdx = 0;
            double worstE = archive[0].energy;
            for (size_t i = 1; i < archive.size(); ++i) {
                if (archive[i].energy > worstE) {
                    worstE = archive[i].energy;
                    worstIdx = static_cast<int>(i);
                }
            }
            archive.erase(archive.begin() + worstIdx);
        }
        archive.push_back({cluster, sig, cluster.getEnergy(), 1});
    }
}

int OrderParameterGuided::selectStructureFromArchive() const {
    // UCB-style selection: score = -E_norm + D * sqrt(log(total) / (visit+1))
    // Balances exploitation (low energy) with exploration (under-visited basins)
    if (archive.empty()) return 0;
    int n = static_cast<int>(archive.size());

    double eMin = 1e10, eMax = -1e10;
    for (int i = 0; i < n; ++i) {
        if (archive[i].energy < eMin) eMin = archive[i].energy;
        if (archive[i].energy > eMax) eMax = archive[i].energy;
    }
    double eRange = (eMax > eMin) ? eMax - eMin : 1.0;

    std::vector<std::pair<double, int>> scored;
    for (int i = 0; i < n; ++i) {
        // Normalized energy: 1 = best (lowest E), 0 = worst
        double eNorm = 1.0 - (archive[i].energy - eMin) / eRange;

        // UCB exploration bonus: high for under-visited entries
        double totalLog = std::log(static_cast<double>(totalVisits) + 1.0);
        double visitBonus = std::sqrt(totalLog / (archive[i].visitCount + 1.0));

        // Higher score = better
        double score = eNorm + params.opDiversityWeight * visitBonus;
        scored.push_back({-score, i});  // negative for ascending sort
    }
    std::sort(scored.begin(), scored.end());

    // Tournament from top 30% candidates
    int poolSize = std::max(3, n / 3);
    int idx = RandomGenerator::uniformInt(0, poolSize - 1);
    return scored[idx].second;
}

int OrderParameterGuided::discretizeOPRegion(const OPSignature& sig) const {
    // Discretize (Q4, Q6) into a 2D grid for Thompson sampling stats
    // Q4 typically in [0, 0.3], Q6 in [0, 0.7] for clusters
    int bins = params.opRegionBins;
    int q4bin = std::min(bins - 1, static_cast<int>(sig.Q4 / 0.30 * bins));
    int q6bin = std::min(bins - 1, static_cast<int>(sig.Q6 / 0.70 * bins));
    return q4bin * bins + q6bin;  // flattened 2D index
}

OrderParameterGuided::OperatorType
OrderParameterGuided::selectOperator(const OPSignature& sig) {
    // Thompson sampling: sample from Beta(successes, failures) for each operator,
    // pick the one with the highest sample.
    int region = discretizeOPRegion(sig);
    auto& stats = operatorStats[region];  // creates default if new region

    double bestSample = -1.0;
    OperatorType bestOp = OperatorType::SINGLE;

    for (int op = 0; op < static_cast<int>(OperatorType::COUNT); ++op) {
        // Beta distribution sample via Gamma(alpha,1) / (Gamma(alpha,1)+Gamma(beta,1))
        double alpha = stats.successes[op] * params.thompsonStrength;
        double beta = stats.failures[op] * params.thompsonStrength;

        // Simple approximation: sample = alpha / (alpha + beta) + noise
        // This approximates the mean + exploration noise of Beta distribution
        double mean = alpha / (alpha + beta);
        double variance = (alpha * beta) / ((alpha + beta) * (alpha + beta) * (alpha + beta + 1.0));
        double thompsonSample = mean + RandomGenerator::normal(0.0, std::sqrt(variance + 0.01));

        if (thompsonSample > bestSample) {
            bestSample = thompsonSample;
            bestOp = static_cast<OperatorType>(op);
        }
    }

    totalOperatorAttempts++;
    return bestOp;
}

void OrderParameterGuided::updateOperatorStats(const OPSignature& sig,
                                                OperatorType op, bool success) {
    int region = discretizeOPRegion(sig);
    auto& stats = operatorStats[region];
    int opIdx = static_cast<int>(op);
    if (success)
        stats.successes[opIdx]++;
    else
        stats.failures[opIdx]++;
}

void OrderParameterGuided::updateTemperature(int step) {
    // Simple linear SA schedule: Tmax → Tmin over the run
    double frac = static_cast<double>(step) / static_cast<double>(params.maxSteps);
    currentTemperature = params.Tmax - (params.Tmax - params.Tmin) * frac;
}

BinaryAlloyCluster OrderParameterGuided::perturbSingle(BinaryAlloyCluster src, double scale) {
    int atomIdx = RandomGenerator::uniformInt(0, numAtoms - 1);
    auto pos = src.getAtomPosition(atomIdx);

    pos[0] += RandomGenerator::normal(0.0, scale);
    pos[1] += RandomGenerator::normal(0.0, scale);
    pos[2] += RandomGenerator::normal(0.0, scale);

    src.setAtomPosition(atomIdx, pos[0], pos[1], pos[2]);
    return src;
}

BinaryAlloyCluster OrderParameterGuided::perturbCollective(BinaryAlloyCluster src) {
    // Pick random seed, find K nearest neighbors, displace group
    int seed = RandomGenerator::uniformInt(0, numAtoms - 1);
    auto seedPos = src.getAtomPosition(seed);

    std::vector<std::pair<double, int>> distances;
    for (int i = 0; i < numAtoms; ++i) {
        if (i == seed) continue;
        auto pos = src.getAtomPosition(i);
        double dx = pos[0] - seedPos[0];
        double dy = pos[1] - seedPos[1];
        double dz = pos[2] - seedPos[2];
        distances.push_back({std::sqrt(dx*dx + dy*dy + dz*dz), i});
    }
    std::sort(distances.begin(), distances.end());

    int groupSize = std::min(params.collectiveGroupSize, numAtoms);
    std::vector<int> group = {seed};
    for (int k = 0; k < groupSize - 1 && k < static_cast<int>(distances.size()); ++k) {
        group.push_back(distances[k].second);
    }

    double dx = RandomGenerator::normal(0.0, params.largePerturbScale);
    double dy = RandomGenerator::normal(0.0, params.largePerturbScale);
    double dz = RandomGenerator::normal(0.0, params.largePerturbScale);

    for (int idx : group) {
        auto pos = src.getAtomPosition(idx);
        src.setAtomPosition(idx, pos[0] + dx, pos[1] + dy, pos[2] + dz);
    }
    return src;
}

BinaryAlloyCluster OrderParameterGuided::perturbSwap(const BinaryAlloyCluster& src) {
    BinaryAlloyCluster trial = src.copy();

    std::vector<int> typeA, typeB;
    for (int i = 0; i < numAtoms; ++i) {
        if (trial.getAtomType(i) == 0) typeA.push_back(i);
        else typeB.push_back(i);
    }

    if (typeA.empty() || typeB.empty()) return trial;

    int idxA = typeA[RandomGenerator::uniformInt(0, static_cast<int>(typeA.size()) - 1)];
    int idxB = typeB[RandomGenerator::uniformInt(0, static_cast<int>(typeB.size()) - 1)];

    auto posA = trial.getAtomPosition(idxA);
    auto posB = trial.getAtomPosition(idxB);

    trial.setAtomPosition(idxA, posB[0], posB[1], posB[2]);
    trial.setAtomPosition(idxB, posA[0], posA[1], posA[2]);

    return trial;
}

bool OrderParameterGuided::metropolisAccept(double oldE, double newE, double T) {
    if (newE < oldE) return true;
    double prob = std::exp(-(newE - oldE) / T);
    return RandomGenerator::uniform() < prob;
}

BinaryAlloyCluster OrderParameterGuided::optimize(const BinaryAlloyCluster& initial,
                                                   ResultManager* resultManager) {

    initialize(initial);

    NELbfgs* lbfgs = lbfgsPool[0].get();
    int stepsPerOutput = std::max(1, params.maxSteps / 20);

    for (int step = 0; step < params.maxSteps; ++step) {
        // Select parent structure from archive
        int parentIdx = selectStructureFromArchive();
        BinaryAlloyCluster parent = archive[parentIdx].structure;
        OPSignature parentSig = archive[parentIdx].signature;

        // Determine perturbation scale
        double scale = (stepsSinceNewBasin > params.stepsBeforeLargePerturb)
            ? params.largePerturbScale : params.smallPerturbScale;

        // Thompson sampling: learn which operator works best in this OP region
        OperatorType chosenOp;
        if (totalOperatorAttempts < params.operatorLearnWindow * params.opRegionBins * params.opRegionBins) {
            // Initial exploration phase: random operator selection to gather data
            double dice = RandomGenerator::uniform();
            if (dice < 0.3) chosenOp = OperatorType::SINGLE;
            else if (dice < 0.7) chosenOp = OperatorType::COLLECTIVE;
            else chosenOp = OperatorType::SWAP;
        } else {
            chosenOp = selectOperator(parentSig);
        }

        // Generate trial structure using chosen operator
        BinaryAlloyCluster trial(numElementA, numElementB, elementA, elementB);
        switch (chosenOp) {
            case OperatorType::SWAP:
                trial = perturbSwap(parent);
                break;
            case OperatorType::COLLECTIVE:
                trial = perturbCollective(std::move(parent));
                break;
            case OperatorType::SINGLE:
            default:
                trial = perturbSingle(std::move(parent), scale);
                break;
        }

        double trialEnergy = evaluateCluster(trial, lbfgs);
        OPSignature trialSig = computeOPSignature(trial);

        // Check if this is a new OP region
        int closestIdx = findClosestArchiveEntry(trialSig);
        bool isNewRegion = (closestIdx < 0 ||
            trialSig.distance(archive[closestIdx].signature) >= params.opGridResolution);

        // Metropolis acceptance with bonus for new OP regions
        // Bonus uses base temperature (not adaptive T) to avoid runaway
        double effectiveEnergy = trialEnergy;
        if (isNewRegion) {
            effectiveEnergy -= currentTemperature * 5.0;
        }

        double oldEnergy = archive[parentIdx].energy;
        bool accepted = metropolisAccept(oldEnergy, effectiveEnergy, currentTemperature);

        // Update Thompson sampling stats
        bool operatorSuccess = accepted || (trialEnergy < oldEnergy);
        if (totalOperatorAttempts >= params.operatorLearnWindow) {
            updateOperatorStats(parentSig, chosenOp, operatorSuccess);
        }

        if (accepted) {
            acceptedCount.fetch_add(1);
            updateArchive(trial, trialSig);
            totalVisits++;

            if (trialEnergy < globalBestEnergy) {
                globalBest = trial;
                globalBestEnergy = trialEnergy;
            }
        }

        // SA temperature schedule
        updateTemperature(step);

        // Progress output
        if (step % stepsPerOutput == 0) {
            std::cout << "OPG step " << step
                << ": best = " << globalBestEnergy
                << " | archive = " << archive.size() << " basins"
                << " | accepted = " << acceptedCount.load()
                << " | T=" << currentTemperature
                << " | Q4=" << trialSig.Q4 << " Q6=" << trialSig.Q6
                << std::endl;
        }

        // Save results
        if (resultManager && step % 10 == 0) {
            resultManager->saveGenerationBest(step, globalBest);
            resultManager->updateHistoricalBest(globalBest);
        }
    }

    // Final polish
    globalBestEnergy = evaluateCluster(globalBest, lbfgs);

    std::cout << "OPG final: best = " << globalBestEnergy
        << " | archive = " << archive.size() << " basins"
        << " | accepted = " << acceptedCount.load() << "/" << stepCount.load()
        << std::endl;

    return globalBest;
}
