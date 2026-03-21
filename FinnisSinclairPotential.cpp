#include "../include/FinnisSinclairPotential.h"
#include <cmath>
#include <fstream>
#include <sstream>

FinnisSinclairPotential::FinnisSinclairPotential() {
    elementA = "A";
    elementB = "B";
    paramsAA = FinnisSinclairParameters(5.0, 1.0, -0.5, 0.1, 4.5, 1.0, 0.0);
    paramsBB = FinnisSinclairParameters(5.0, 1.0, -0.5, 0.1, 4.5, 1.0, 0.0);
    paramsAB = FinnisSinclairParameters(5.0, 1.0, -0.5, 0.1, 4.5, 1.0, 0.0);
}

FinnisSinclairPotential::FinnisSinclairPotential(const std::string& elemA,
    const std::string& elemB) {
    elementA = elemA;
    elementB = elemB;
    paramsAA = FinnisSinclairParameters(5.0, 1.0, -0.5, 0.1, 4.5, 1.0, 0.0);
    paramsBB = FinnisSinclairParameters(5.0, 1.0, -0.5, 0.1, 4.5, 1.0, 0.0);
    paramsAB = FinnisSinclairParameters(5.0, 1.0, -0.5, 0.1, 4.5, 1.0, 0.0);
}

bool FinnisSinclairPotential::loadParameters(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line;
    std::vector<std::vector<double>> params;

    // Read first line - may contain element symbols
    if (std::getline(file, line)) {
        if (line.find(',') == std::string::npos && line.find('#') != 0) {
            std::istringstream iss(line);
            iss >> elementA >> elementB;
        }
        else if (line.find('#') != 0) {
            file.clear();
            file.seekg(0);
        }
    }

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::vector<double> row;
        std::string value;

        while (std::getline(iss, value, ',')) {
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);

            try {
                row.push_back(std::stod(value));
            }
            catch (...) {
                continue;
            }
        }

        if (row.size() == 7) {  // cutoff, c0, c1, c2, d, c, beta
            params.push_back(row);
        }
    }

    if (params.size() >= 3) {
        paramsAA = FinnisSinclairParameters(params[0][0], params[0][1], params[0][2],
            params[0][3], params[0][4], params[0][5], params[0][6]);
        paramsBB = FinnisSinclairParameters(params[1][0], params[1][1], params[1][2],
            params[1][3], params[1][4], params[1][5], params[1][6]);
        paramsAB = FinnisSinclairParameters(params[2][0], params[2][1], params[2][2],
            params[2][3], params[2][4], params[2][5], params[2][6]);
        return true;
    }

    return false;
}

void FinnisSinclairPotential::setParameters(const FinnisSinclairParameters& aa,
    const FinnisSinclairParameters& bb,
    const FinnisSinclairParameters& ab) {
    paramsAA = aa;
    paramsBB = bb;
    paramsAB = ab;
}

void FinnisSinclairPotential::setElements(const std::string& elemA,
    const std::string& elemB) {
    elementA = elemA;
    elementB = elemB;
}

void FinnisSinclairPotential::computeDistanceMatrix(const BinaryAlloyCluster& cluster) const {
    int n = cluster.getNumAtoms();
    if (distanceMatrix.size() != static_cast<size_t>(n * n)) {
        distanceMatrix.resize(n * n);
    }

    const double* x = cluster.data();
    const double* y = x + n;
    const double* z = y + n;

    for (int i = 0; i < n - 1; ++i) {
        distanceMatrix[i * n + i] = 0.0;
        for (int j = i + 1; j < n; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            distanceMatrix[i * n + j] = r;
            distanceMatrix[j * n + i] = r;
        }
    }
    distanceMatrix[n * n - 1] = 0.0;
}

FinnisSinclairParameters FinnisSinclairPotential::getParameters(
    const BinaryAlloyCluster& cluster, int i, int j) const {
    int typeI = cluster.getAtomType(i);
    int typeJ = cluster.getAtomType(j);

    if (typeI == 0 && typeJ == 0) {
        return paramsAA;
    }
    else if (typeI == 1 && typeJ == 1) {
        return paramsBB;
    }
    else {
        return paramsAB;
    }
}

// Pair potential: V(r) = (r - cutoff)^2 * (c0 + c1*r + c2*r^2) for r < cutoff
double FinnisSinclairPotential::pairPotential(double r,
    const FinnisSinclairParameters& params) const {
    if (r >= params.cutoff) return 0.0;

    double dr = r - params.cutoff;
    double poly = params.c0 + params.c1 * r + params.c2 * r * r;
    return dr * dr * poly;
}

// Pair potential derivative: dV/dr
double FinnisSinclairPotential::pairPotentialDerivative(double r,
    const FinnisSinclairParameters& params) const {
    if (r >= params.cutoff) return 0.0;

    double dr = r - params.cutoff;
    double poly = params.c0 + params.c1 * r + params.c2 * r * r;
    double dPoly = params.c1 + 2.0 * params.c2 * r;

    return 2.0 * dr * poly + dr * dr * dPoly;
}

// Density function: phi(r) = (r - d)^2 + beta*(r-d)^3/d for r < d
double FinnisSinclairPotential::densityFunction(double r,
    const FinnisSinclairParameters& params) const {
    if (r >= params.d) return 0.0;

    double dr = r - params.d;
    return dr * dr + params.beta * dr * dr * dr / params.d;
}

// Density function derivative: dphi/dr = 2(r-d) + 3*beta*(r-d)^2/d
double FinnisSinclairPotential::densityFunctionDerivative(double r,
    const FinnisSinclairParameters& params) const {
    if (r >= params.d) return 0.0;

    double dr = r - params.d;
    return 2.0 * dr + 3.0 * params.beta * dr * dr / params.d;
}

double FinnisSinclairPotential::calculateEnergy(const BinaryAlloyCluster& cluster) {
    int n = cluster.getNumAtoms();
    computeDistanceMatrix(cluster);

    double totalEnergy = 0.0;

    // Calculate density for each atom
    std::vector<double> rho(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;

            double r = distanceMatrix[i * n + j];
            FinnisSinclairParameters params = getParameters(cluster, i, j);
            rho[i] += densityFunction(r, params);
        }
    }

    // Calculate total energy
    for (int i = 0; i < n; ++i) {
        double pairEnergy = 0.0;

        for (int j = 0; j < n; ++j) {
            if (i == j) continue;

            double r = distanceMatrix[i * n + j];
            FinnisSinclairParameters params = getParameters(cluster, i, j);
            pairEnergy += pairPotential(r, params);
        }

        FinnisSinclairParameters params = getParameters(cluster, i, i);
        double embeddingEnergy = (rho[i] > Constants::EPSILON) ?
            params.c * std::sqrt(rho[i]) : 0.0;

        totalEnergy += 0.5 * pairEnergy - embeddingEnergy;
    }

    return totalEnergy;
}

void FinnisSinclairPotential::calculateForces(const BinaryAlloyCluster& cluster,
    std::vector<double>& f) {
    int n = cluster.getNumAtoms();
    computeDistanceMatrix(cluster);

    if (f.size() != static_cast<size_t>(3 * n)) {
        f.resize(3 * n);
    }
    std::fill(f.begin(), f.end(), 0.0);

    // Calculate density
    std::vector<double> rho(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;

            double r = distanceMatrix[i * n + j];
            FinnisSinclairParameters params = getParameters(cluster, i, j);
            rho[i] += densityFunction(r, params);
        }
    }

    // Calculate embedding energy derivative coefficients
    std::vector<double> dRho(n);
    for (int i = 0; i < n; ++i) {
        FinnisSinclairParameters params = getParameters(cluster, i, i);
        if (rho[i] > Constants::EPSILON) {
            dRho[i] = -params.c / (2.0 * std::sqrt(rho[i]));
        }
        else {
            dRho[i] = 0.0;
        }
    }

    const double* x = cluster.data();
    const double* y = x + n;
    const double* z = y + n;

    double* fx = f.data();
    double* fy = fx + n;
    double* fz = fy + n;

    // Calculate forces
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double r = distanceMatrix[i * n + j];
            if (r < Constants::EPSILON) continue;

            FinnisSinclairParameters params = getParameters(cluster, i, j);

            // Pair potential contribution
            double dV = pairPotentialDerivative(r, params);

            // Density term contribution
            double dPhi = densityFunctionDerivative(r, params);
            double embeddingForce = (dRho[i] + dRho[j]) * dPhi;

            // Total radial force component
            double forceRadial = -dV - embeddingForce;

            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];

            double forceScale = forceRadial / r;

            fx[i] += forceScale * dx;
            fx[j] -= forceScale * dx;
            fy[i] += forceScale * dy;
            fy[j] -= forceScale * dy;
            fz[i] += forceScale * dz;
            fz[j] -= forceScale * dz;
        }
    }
}

double FinnisSinclairPotential::calculateEnergyWithForces(const BinaryAlloyCluster& cluster,
    std::vector<double>& f) {
    calculateForces(cluster, f);
    return calculateEnergy(cluster);
}