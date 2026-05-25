#include "../include/SuttonChenPotential.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iostream>

// 1. 构造函数
SuttonChenPotential::SuttonChenPotential() {
    elementA = "Pt";
    elementB = "Cu";
    // Pt-Pt: ε, C, a, n, m
    paramsAA = SuttonChenParameters(0.0097894, 71.336, 3.9163, 11.0, 7.0);
    // Cu-Cu
    paramsBB = SuttonChenParameters(0.0057921, 84.843, 3.6030, 10.0, 5.0);
    // Pt-Cu (几何/算术平均)
    paramsAB = SuttonChenParameters(0.0075300, 78.0895, 3.75965, 10.5, 6.0);
    cutoff = 10.0;
}

SuttonChenPotential::SuttonChenPotential(const std::string& elemA, const std::string& elemB) {
    elementA = elemA;
    elementB = elemB;
    // Pt-Pt: ε, C, a, n, m
    paramsAA = SuttonChenParameters(0.0097894, 71.336, 3.9163, 11.0, 7.0);
    // Cu-Cu
    paramsBB = SuttonChenParameters(0.0057921, 84.843, 3.6030, 10.0, 5.0);
    // Pt-Cu (几何/算术平均)
    paramsAB = SuttonChenParameters(0.0075300, 78.0895, 3.75965, 10.5, 6.0);
    cutoff = 10.0;
}

// 2. 参数加载函数
bool SuttonChenPotential::loadParameters(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line;
    std::vector<std::vector<double>> params_list;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        // 处理可能的标签
        size_t bracketPos = line.find_last_of(']');
        std::string dataPart;
        if (bracketPos != std::string::npos) {
            dataPart = line.substr(bracketPos + 1);
        }
        else {
            dataPart = line;
        }

        std::istringstream iss(dataPart);
        std::vector<double> row;
        std::string value;

        while (std::getline(iss, value, ',')) {
            size_t first = value.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) continue;
            size_t last = value.find_last_not_of(" \t\r\n");
            value = value.substr(first, (last - first + 1));

            try {
                row.push_back(std::stod(value));
            }
            catch (...) {
                continue;
            }
        }

        if (row.size() >= 5) {
            params_list.push_back(row);
        }
    }

    if (params_list.size() >= 3) {
        paramsAA = SuttonChenParameters(params_list[0][0], params_list[0][1],
            params_list[0][2], params_list[0][3], params_list[0][4]);
        paramsBB = SuttonChenParameters(params_list[1][0], params_list[1][1],
            params_list[1][2], params_list[1][3], params_list[1][4]);
        paramsAB = SuttonChenParameters(params_list[2][0], params_list[2][1],
            params_list[2][2], params_list[2][3], params_list[2][4]);
        return true;
    }

    return false;
}

// 3. 势函数
double SuttonChenPotential::pairPotential(double r, const SuttonChenParameters& params) const {
    return std::pow(params.a / r, params.n);
}

double SuttonChenPotential::pairPotentialDerivative(double r, const SuttonChenParameters& params) const {
    return -params.n * std::pow(params.a / r, params.n) / r;
}

double SuttonChenPotential::densityFunction(double r, const SuttonChenParameters& params) const {
    return std::pow(params.a / r, params.m);
}

double SuttonChenPotential::densityFunctionDerivative(double r, const SuttonChenParameters& params) const {
    return -params.m * std::pow(params.a / r, params.m) / r;
}

// 4. 距离矩阵计算 - 与Gupta完全一致
void SuttonChenPotential::computeDistanceMatrix(const BinaryAlloyCluster& cluster, std::vector<double>& dist) const {
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

// 5. 能量计算 - 完全模仿Gupta的双重遍历结构
double SuttonChenPotential::calcEnergyWithDist(const BinaryAlloyCluster& cluster, const std::vector<double>& dist) const {
    int n = cluster.getNumAtoms();
    double totalEnergy = 0.0;

    for (int i = 0; i < n; ++i) {
        double repulsive = 0.0;
        double densitySum = 0.0;

        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double r = dist[i * n + j];
            SuttonChenParameters params = getParameters(cluster, i, j);
            repulsive += params.epsilon * pairPotential(r, params);
            densitySum += densityFunction(r, params);
        }

        SuttonChenParameters paramsI = getParameters(cluster, i, i);
        totalEnergy += 0.5 * repulsive;
        if (densitySum > 1e-15) {
            totalEnergy -= paramsI.epsilon * paramsI.c * std::sqrt(densitySum);
        }
    }

    return totalEnergy;
}

double SuttonChenPotential::calculateEnergy(const BinaryAlloyCluster& cluster) {
    std::vector<double> dist;
    computeDistanceMatrix(cluster, dist);
    return calcEnergyWithDist(cluster, dist);
}

// 6. 力计算 - 完全模仿Gupta的结构
void SuttonChenPotential::calcForcesWithDist(const BinaryAlloyCluster& cluster, std::vector<double>& f, const std::vector<double>& dist) const {
    int n = cluster.getNumAtoms();
    f.assign(3 * n, 0.0);

    std::vector<double> rho(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double r = dist[i * n + j];
            SuttonChenParameters params = getParameters(cluster, i, j);
            rho[i] += densityFunction(r, params);
        }
    }

    std::vector<double> coeff(n);
    for (int i = 0; i < n; ++i) {
        if (rho[i] > 1e-15) {
            SuttonChenParameters paramsI = getParameters(cluster, i, i);
            coeff[i] = paramsI.epsilon * paramsI.c / (2.0 * std::sqrt(rho[i]));
        } else {
            coeff[i] = 0.0;
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
            if (r < 1e-15) continue;

            SuttonChenParameters params = getParameters(cluster, i, j);

            double dPair = params.epsilon * pairPotentialDerivative(r, params);
            double dDens = densityFunctionDerivative(r, params);
            double force = -dPair + (coeff[i] + coeff[j]) * dDens;

            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];

            fx[i] += force * dx / r;
            fx[j] -= force * dx / r;
            fy[i] += force * dy / r;
            fy[j] -= force * dy / r;
            fz[i] += force * dz / r;
            fz[j] -= force * dz / r;
        }
    }
}

void SuttonChenPotential::setElements(const std::string& elemA, const std::string& elemB) {
    elementA = elemA;
    elementB = elemB;
}

SuttonChenParameters SuttonChenPotential::getParameters(const BinaryAlloyCluster& cluster,
    int i, int j) const {
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

double SuttonChenPotential::calculateEnergyWithForces(const BinaryAlloyCluster& cluster,
    std::vector<double>& f) {
    std::vector<double> dist;
    computeDistanceMatrix(cluster, dist);
    calcForcesWithDist(cluster, f, dist);
    return calcEnergyWithDist(cluster, dist);
}
