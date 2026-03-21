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
    paramsAA = SuttonChenParameters(0.01, 34.408, 3.9163, 10.0, 8.0);
    paramsBB = SuttonChenParameters(0.0057, 39.432, 3.61, 9.0, 6.0);
    paramsAB = SuttonChenParameters(0.00755, 36.920, 3.76315, 9.5, 7.0);
    cutoff = 10.0;
}

SuttonChenPotential::SuttonChenPotential(const std::string& elemA, const std::string& elemB) {
    elementA = elemA;
    elementB = elemB;
    paramsAA = SuttonChenParameters(0.01, 34.408, 3.9163, 10.0, 8.0);
    paramsBB = SuttonChenParameters(0.0057, 39.432, 3.61, 9.0, 6.0);
    paramsAB = SuttonChenParameters(0.00755, 36.920, 3.76315, 9.5, 7.0);
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
void SuttonChenPotential::computeDistanceMatrix(const BinaryAlloyCluster& cluster) const {
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

// 5. 能量计算 - 完全模仿Gupta的双重遍历结构
double SuttonChenPotential::calculateEnergy(const BinaryAlloyCluster& cluster) {
    int n = cluster.getNumAtoms();
    computeDistanceMatrix(cluster);

    double totalEnergy = 0.0;

    // 【关键】：完全模仿Gupta的结构
    for (int i = 0; i < n; ++i) {
        double repulsive = 0.0;      // 对势项（对应Gupta的repulsive）
        double densitySum = 0.0;     // 密度项（对应Gupta的attractive）

        for (int j = 0; j < n; ++j) {
            if (i == j) continue;

            double r = distanceMatrix[i * n + j];
            SuttonChenParameters params = getParameters(cluster, i, j);

            // 累加对势和密度
            repulsive += params.epsilon * pairPotential(r, params);
            densitySum += densityFunction(r, params);
        }

        // 每个原子的能量贡献
        // E_i = (1/2) * Σⱼ V(rᵢⱼ) - c*ε*√ρᵢ
        // 注意：对势要除以2，因为双重遍历会重复计数
        SuttonChenParameters paramsI = getParameters(cluster, i, i);
        totalEnergy += 0.5 * repulsive;
        if (densitySum > 1e-15) {
            totalEnergy -= paramsI.epsilon * paramsI.c * std::sqrt(densitySum);
        }
    }

    return totalEnergy;
}

// 6. 力计算 - 完全模仿Gupta的结构
void SuttonChenPotential::calculateForces(const BinaryAlloyCluster& cluster, std::vector<double>& f) {
    int n = cluster.getNumAtoms();
    computeDistanceMatrix(cluster);

    if (f.size() != static_cast<size_t>(3 * n)) {
        f.resize(3 * n);
    }
    std::fill(f.begin(), f.end(), 0.0);

    // 【步骤1】：计算密度 - 与Gupta一致
    std::vector<double> rho(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;

            double r = distanceMatrix[i * n + j];
            SuttonChenParameters params = getParameters(cluster, i, j);
            rho[i] += densityFunction(r, params);
        }
    }

    // 【步骤2】：计算系数 = 嵌入能对密度的导数 - 与Gupta一致
    std::vector<double> coeff(n);
    for (int i = 0; i < n; ++i) {
        if (rho[i] > 1e-15) {
            SuttonChenParameters paramsI = getParameters(cluster, i, i);
            // coeff[i] = -∂E/∂ρᵢ = ε*c/(2√ρᵢ)
            coeff[i] = paramsI.epsilon * paramsI.c / (2.0 * std::sqrt(rho[i]));
        }
        else {
            coeff[i] = 0.0;
        }
    }

    // 【步骤3】：只遍历i<j计算力 - 与Gupta一致
    const double* x = cluster.data();
    const double* y = x + n;
    const double* z = y + n;

    double* fx = f.data();
    double* fy = fx + n;
    double* fz = fy + n;

    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double r = distanceMatrix[i * n + j];
            if (r < 1e-15) continue;

            SuttonChenParameters params = getParameters(cluster, i, j);

            // 对势项的导数
            double dPair = params.epsilon * pairPotentialDerivative(r, params);

            // 密度项的导数
            double dDens = densityFunctionDerivative(r, params);

            // 总力 = -对势导数 + 密度导数的贡献
            // 注意：与Gupta不同，SC的密度函数导数不需要除以2
            double force = -dPair + (coeff[i] + coeff[j]) * dDens;

            // 分解到xyz方向
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

// 辅助函数
void SuttonChenPotential::setElements(const std::string& elemA, const std::string& elemB) {
    elementA = elemA;
    elementB = elemB;
}

void SuttonChenPotential::setParameters(const SuttonChenParameters& aa,
    const SuttonChenParameters& bb,
    const SuttonChenParameters& ab) {
    paramsAA = aa;
    paramsBB = bb;
    paramsAB = ab;
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
    calculateForces(cluster, f);
    return calculateEnergy(cluster);
}
