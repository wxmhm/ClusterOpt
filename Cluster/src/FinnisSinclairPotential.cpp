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
    // 默认初始化防止未加载参数时崩溃
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

    // 尝试读取首行元素名，增加鲁棒性
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

        // 兼容带标签的行，如 "Pt Cu"
        std::string content = line;
        size_t bracketPos = line.find(']');
        if (bracketPos != std::string::npos) {
            content = line.substr(bracketPos + 1);
        }

        std::istringstream iss(content);
        std::vector<double> row;
        std::string value;

        while (std::getline(iss, value, ',')) {
            // 清理空白字符
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

        // FS参数通常有7个: cutoff, c0, c1, c2, d, c, beta
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

// 获取参数的辅助函数
// 对于对势(Pair Potential)，如果是 AB 对，需要返回混合参数 paramsAB
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

double FinnisSinclairPotential::calculateEnergy(const BinaryAlloyCluster& cluster) {
    int n = cluster.getNumAtoms();
    computeDistanceMatrix(cluster);

    double totalEnergy = 0.0;

    // 1. 计算每个原子处的电子密度 rho_i
    std::vector<double> rho(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double r = distanceMatrix[i * n + j];

            // 【关键修正】：密度贡献取决于邻居原子 j 的类型
            // 理论依据：rho_i = sum( phi_j(r_ij) )
            // 即使 i 和 j 不同，密度函数也只使用 j 的参数
            int typeJ = cluster.getAtomType(j);
            const FinnisSinclairParameters& paramsJ = (typeJ == 0) ? paramsAA : paramsBB;

            rho[i] += densityFunction(r, paramsJ);
        }
    }

    // 2. 计算总能量
    for (int i = 0; i < n; ++i) {
        double pairEnergy = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double r = distanceMatrix[i * n + j];

            // 对势 V_ij 依然使用 AB 混合参数
            FinnisSinclairParameters paramsPair = getParameters(cluster, i, j);
            pairEnergy += pairPotential(r, paramsPair);
        }

        // 嵌入能 E_emb = -A * sqrt(rho)
        // 使用原子 i 自身的参数 A (即代码中的 c)
        FinnisSinclairParameters paramsI = getParameters(cluster, i, i);
        double embeddingEnergy = 0.0;
        if (rho[i] > Constants::EPSILON) {
            embeddingEnergy = -paramsI.c * std::sqrt(rho[i]);
        }

        // E_total = 1/2 * Sum(V_ij) + Sum(E_emb_i)
        totalEnergy += 0.5 * pairEnergy + embeddingEnergy;
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

    // 1. 计算密度 (与 calculateEnergy 逻辑必须完全一致)
    std::vector<double> rho(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double r = distanceMatrix[i * n + j];
            int typeJ = cluster.getAtomType(j);
            const FinnisSinclairParameters& paramsJ = (typeJ == 0) ? paramsAA : paramsBB;
            rho[i] += densityFunction(r, paramsJ);
        }
    }

    // 2. 预计算嵌入能对密度的导数 dE/drho
    // E = -c * sqrt(rho) => dE/drho = -c / (2 * sqrt(rho))
    std::vector<double> dEdRho(n);
    for (int i = 0; i < n; ++i) {
        FinnisSinclairParameters paramsI = getParameters(cluster, i, i);
        if (rho[i] > Constants::EPSILON) {
            dEdRho[i] = -paramsI.c / (2.0 * std::sqrt(rho[i]));
        }
        else {
            dEdRho[i] = 0.0;
        }
    }

    const double* x = cluster.data();
    const double* y = x + n;
    const double* z = y + n;

    double* fx = f.data();
    double* fy = fx + n;
    double* fz = fy + n;

    // 3. 计算力 (遍历每一对 i, j)
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double r = distanceMatrix[i * n + j];
            if (r < Constants::EPSILON) continue;

            // --- A. 对势部分的力 ---
            // F_pair = -dV/dr
            FinnisSinclairParameters paramsPair = getParameters(cluster, i, j);
            double dV = pairPotentialDerivative(r, paramsPair);

            // --- B. 嵌入能部分的力 (多体效应) ---
            // 链式法则：F_emb = (dE/drho_i * drho_i/dr) + (dE/drho_j * drho_j/dr)

            // 1. 原子 j 对原子 i 的密度贡献 rho_i 的导数 -> phi'_j(r)
            // 需要使用 j 的参数
            int typeJ = cluster.getAtomType(j);
            const FinnisSinclairParameters& paramsJ = (typeJ == 0) ? paramsAA : paramsBB;
            double dPhi_j = densityFunctionDerivative(r, paramsJ);

            // 2. 原子 i 对原子 j 的密度贡献 rho_j 的导数 -> phi'_i(r)
            // 需要使用 i 的参数
            int typeI = cluster.getAtomType(i);
            const FinnisSinclairParameters& paramsI = (typeI == 0) ? paramsAA : paramsBB;
            double dPhi_i = densityFunctionDerivative(r, paramsI);

            // 组合嵌入力项
            // 注意：dEdRho 通常是负的，dPhi 也是负的(密度随距离衰减)
            // 所以 dEdRho * dPhi 是正值
            double embeddingForceTerm = dEdRho[i] * dPhi_j + dEdRho[j] * dPhi_i;

            // 总径向力导数 F_radial = -dE_total/dr
            // = -(dV + embeddingForceTerm)
            // = -dV - embeddingForceTerm
            double forceRadial = -dV - embeddingForceTerm;

            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
            double invR = 1.0 / r;

            // 将径向力分解到 XYZ 分量
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

double FinnisSinclairPotential::calculateEnergyWithForces(const BinaryAlloyCluster& cluster,
    std::vector<double>& f) {
    calculateForces(cluster, f);
    return calculateEnergy(cluster);
}