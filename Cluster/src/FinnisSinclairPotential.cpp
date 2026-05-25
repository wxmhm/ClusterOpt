#include "../include/FinnisSinclairPotential.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iostream>

FinnisSinclairPotential::FinnisSinclairPotential() {
    elementA = "A";
    elementB = "B";
    // Ĭ�ϳ�ʼ����ֹδ���ز���ʱ����
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

    // ���Զ�ȡ����Ԫ����������³����
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

        // ���ݴ���ǩ���У��� "Pt Cu"
        std::string content = line;
        size_t bracketPos = line.find(']');
        if (bracketPos != std::string::npos) {
            content = line.substr(bracketPos + 1);
        }

        std::istringstream iss(content);
        std::vector<double> row;
        std::string value;

        while (std::getline(iss, value, ',')) {
            // �����հ��ַ�
            value.erase(0, value.find_first_not_of(" \t\r\n"));
            value.erase(value.find_last_not_of(" \t\r\n") + 1);
            if (value.empty()) continue;

            try {
                row.push_back(std::stod(value));
            }
            catch (...) {
                continue;
            }
        }

        // FS����ͨ����7��: cutoff, c0, c1, c2, d, c, beta
        if (row.size() >= 7) {
            params.push_back(row);
        }
    }

    if (params.size() >= 3) {
        // 0: AA, 1: BB, 2: AB
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

void FinnisSinclairPotential::setElements(const std::string& elemA,
    const std::string& elemB) {
    elementA = elemA;
    elementB = elemB;
}

void FinnisSinclairPotential::computeDistanceMatrix(const BinaryAlloyCluster& cluster, std::vector<double>& dist) const {
    int n = cluster.getNumAtoms();
    dist.resize(n * n);

    const double* x = cluster.data();
    const double* y = x + n;
    const double* z = y + n;

    for (int i = 0; i < n - 1; ++i) {
        dist[i * n + i] = 0.0;
        for (int j = i + 1; j < n; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            dist[i * n + j] = r;
            dist[j * n + i] = r;
        }
    }
    dist[n * n - 1] = 0.0;
}

// ��ȡ�����ĸ�������
// ���ڶ���(Pair Potential)������� AB �ԣ���Ҫ���ػ�ϲ��� paramsAB
FinnisSinclairParameters FinnisSinclairPotential::getParameters(
    const BinaryAlloyCluster& cluster, int i, int j) const {
    int typeI = cluster.getAtomType(i);
    int typeJ = cluster.getAtomType(j);

    if (typeI == 0 && typeJ == 0) return paramsAA;
    else if (typeI == 1 && typeJ == 1) return paramsBB;
    else return paramsAB;
}

// Pair potential: V(r)
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

    // d/dr [(r-c)^2 * P(r)] = 2(r-c)P(r) + (r-c)^2 P'(r)
    return 2.0 * dr * poly + dr * dr * dPoly;
}

// Density function: phi(r)
double FinnisSinclairPotential::densityFunction(double r,
    const FinnisSinclairParameters& params) const {
    if (r >= params.d) return 0.0;

    double dr = r - params.d;
    // phi(r) = (r-d)^2 + beta/d * (r-d)^3
    return dr * dr + params.beta * dr * dr * dr / params.d;
}

// Density function derivative: dphi/dr
double FinnisSinclairPotential::densityFunctionDerivative(double r,
    const FinnisSinclairParameters& params) const {
    if (r >= params.d) return 0.0;

    double dr = r - params.d;
    // dphi/dr = 2(r-d) + 3*beta/d * (r-d)^2
    return 2.0 * dr + 3.0 * params.beta * dr * dr / params.d;
}

double FinnisSinclairPotential::calcEnergyWithDist(const BinaryAlloyCluster& cluster, const std::vector<double>& dist) const {
    int n = cluster.getNumAtoms();
    double totalEnergy = 0.0;

    std::vector<double> rho(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double r = dist[i * n + j];
            int typeJ = cluster.getAtomType(j);
            const FinnisSinclairParameters& paramsJ = (typeJ == 0) ? paramsAA : paramsBB;
            rho[i] += densityFunction(r, paramsJ);
        }
    }

    for (int i = 0; i < n; ++i) {
        double pairEnergy = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double r = dist[i * n + j];
            const FinnisSinclairParameters& paramsPair = getParameters(cluster, i, j);
            pairEnergy += pairPotential(r, paramsPair);
        }

        const FinnisSinclairParameters& paramsI = getParameters(cluster, i, i);
        double embeddingEnergy = 0.0;
        if (rho[i] > Constants::EPSILON) {
            embeddingEnergy = -paramsI.c * std::sqrt(rho[i]);
        }

        totalEnergy += 0.5 * pairEnergy + embeddingEnergy;
    }

    return totalEnergy;
}

void FinnisSinclairPotential::calcForcesWithDist(const BinaryAlloyCluster& cluster,
    std::vector<double>& f, const std::vector<double>& dist) const {
    int n = cluster.getNumAtoms();
    f.assign(3 * n, 0.0);

    std::vector<double> rho(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double r = dist[i * n + j];
            int typeJ = cluster.getAtomType(j);
            const FinnisSinclairParameters& paramsJ = (typeJ == 0) ? paramsAA : paramsBB;
            rho[i] += densityFunction(r, paramsJ);
        }
    }

    std::vector<double> dEdRho(n);
    for (int i = 0; i < n; ++i) {
        const FinnisSinclairParameters& paramsI = getParameters(cluster, i, i);
        if (rho[i] > Constants::EPSILON) {
            dEdRho[i] = -paramsI.c / (2.0 * std::sqrt(rho[i]));
        } else {
            dEdRho[i] = 0.0;
        }
    }

    const double* x = cluster.data();
    const double* y = x + n;
    const double* z = y + n;

    double* fx = f.data();
    double* fy = fx + n;
    double* fz = fy + n;

    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double r = dist[i * n + j];
            if (r < Constants::EPSILON) continue;

            const FinnisSinclairParameters& paramsPair = getParameters(cluster, i, j);
            double dV = pairPotentialDerivative(r, paramsPair);

            int typeJ = cluster.getAtomType(j);
            const FinnisSinclairParameters& paramsJ = (typeJ == 0) ? paramsAA : paramsBB;
            double dPhi_j = densityFunctionDerivative(r, paramsJ);

            int typeI = cluster.getAtomType(i);
            const FinnisSinclairParameters& paramsI = (typeI == 0) ? paramsAA : paramsBB;
            double dPhi_i = densityFunctionDerivative(r, paramsI);

            double embeddingForceTerm = dEdRho[i] * dPhi_j + dEdRho[j] * dPhi_i;

            double forceRadial = -dV - embeddingForceTerm;

            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
            double invR = 1.0 / r;

            double forceScale = forceRadial * invR;

            fx[i] += forceScale * dx;
            fx[j] -= forceScale * dx;
            fy[i] += forceScale * dy;
            fy[j] -= forceScale * dy;
            fz[i] += forceScale * dz;
            fz[j] -= forceScale * dz;
        }
    }
}

double FinnisSinclairPotential::calculateEnergy(const BinaryAlloyCluster& cluster) {
    std::vector<double> dist;
    computeDistanceMatrix(cluster, dist);
    return calcEnergyWithDist(cluster, dist);
}

double FinnisSinclairPotential::calculateEnergyWithForces(const BinaryAlloyCluster& cluster,
    std::vector<double>& f) {
    std::vector<double> dist;
    computeDistanceMatrix(cluster, dist);
    calcForcesWithDist(cluster, f, dist);
    return calcEnergyWithDist(cluster, dist);
}