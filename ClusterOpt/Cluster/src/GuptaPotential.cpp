#include "../include/GuptaPotential.h"

GuptaPotential::GuptaPotential() : elementA("A"), elementB("B") {
    paramsAA = GuptaParameters(2.5, 0.1, 1.0, -10.0, -4.0);
    paramsBB = GuptaParameters(2.5, 0.1, 1.0, -10.0, -4.0);
    paramsAB = GuptaParameters(2.5, 0.1, 1.0, -10.0, -4.0);
}

GuptaPotential::GuptaPotential(const std::string& elemA, const std::string& elemB)
    : elementA(elemA), elementB(elemB) {
    paramsAA = GuptaParameters(2.5, 0.1, 1.0, -10.0, -4.0);
    paramsBB = GuptaParameters(2.5, 0.1, 1.0, -10.0, -4.0);
    paramsAB = GuptaParameters(2.5, 0.1, 1.0, -10.0, -4.0);
}

bool GuptaPotential::loadParameters(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    std::string line;
    std::vector<std::vector<double>> params;
    
    // Read first line - might be element symbols
    if (std::getline(file, line)) {
        if (line.find(',') == std::string::npos && line.find('#') != 0) {
            std::istringstream iss(line);
            iss >> elementA >> elementB;
        } else if (line.find('#') != 0) {
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
            } catch (...) {
                continue;
            }
        }
        
        if (row.size() == 5) {
            params.push_back(row);
        }
    }
    
    if (params.size() >= 3) {
        paramsAA = GuptaParameters(params[0][0], params[0][1], params[0][2],
                                   params[0][3], params[0][4]);
        paramsBB = GuptaParameters(params[1][0], params[1][1], params[1][2],
                                   params[1][3], params[1][4]);
        paramsAB = GuptaParameters(params[2][0], params[2][1], params[2][2],
                                   params[2][3], params[2][4]);
        return true;
    }
    
    return false;
}

void GuptaPotential::setParameters(const GuptaParameters& aa,
                                   const GuptaParameters& bb,
                                   const GuptaParameters& ab) {
    paramsAA = aa;
    paramsBB = bb;
    paramsAB = ab;
}

void GuptaPotential::setElements(const std::string& elemA, const std::string& elemB) {
    elementA = elemA;
    elementB = elemB;
}

void GuptaPotential::computeDistanceMatrix(const BinaryAlloyCluster& cluster) const {
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

GuptaParameters GuptaPotential::getParameters(const BinaryAlloyCluster& cluster, int i, int j) const {
    int typeI = cluster.getAtomType(i);
    int typeJ = cluster.getAtomType(j);
    
    if (typeI == 0 && typeJ == 0) {
        return paramsAA;
    } else if (typeI == 1 && typeJ == 1) {
        return paramsBB;
    } else {
        return paramsAB;
    }
}

double GuptaPotential::calculateEnergy(const BinaryAlloyCluster& cluster) {
    int n = cluster.getNumAtoms();
    computeDistanceMatrix(cluster);
    
    double totalEnergy = 0.0;
    
    for (int i = 0; i < n; ++i) {
        double repulsive = 0.0;
        double attractive = 0.0;
        
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            
            double r = distanceMatrix[i * n + j];
            GuptaParameters params = getParameters(cluster, i, j);
            
            repulsive += params.A * std::exp(params.p * (r / params.r0 - 1.0));
            attractive += params.xi * params.xi * std::exp(2.0 * params.q * (r / params.r0 - 1.0));
        }
        
        totalEnergy += repulsive - std::sqrt(attractive);
    }
    
    return totalEnergy;
}

void GuptaPotential::calculateForces(const BinaryAlloyCluster& cluster, std::vector<double>& f) {
    int n = cluster.getNumAtoms();
    computeDistanceMatrix(cluster);
    
    if (f.size() != static_cast<size_t>(3 * n)) {
        f.resize(3 * n);
    }
    std::fill(f.begin(), f.end(), 0.0);
    
    std::vector<double> rho(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            
            double r = distanceMatrix[i * n + j];
            GuptaParameters params = getParameters(cluster, i, j);
            rho[i] += params.xi * params.xi * std::exp(2.0 * params.q * (r / params.r0 - 1.0));
        }
    }
    
    std::vector<double> coeff(n);
    for (int i = 0; i < n; ++i) {
        coeff[i] = (rho[i] > Constants::EPSILON) ? 1.0 / std::sqrt(rho[i]) : 0.0;
    }
    
    const double* x = cluster.data();
    const double* y = x + n;
    const double* z = y + n;
    
    double* fx = f.data();
    double* fy = fx + n;
    double* fz = fy + n;
    
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double r = distanceMatrix[i * n + j];
            if (r < Constants::EPSILON) continue;
            
            GuptaParameters params = getParameters(cluster, i, j);
            
            double expRep = std::exp(params.p * (r / params.r0 - 1.0));
            double expAtt = std::exp(2.0 * params.q * (r / params.r0 - 1.0));
            
            double dRep = params.A * params.p * expRep / params.r0;
            double dAtt = params.xi * params.xi * params.q * expAtt / params.r0;
            
            double force = -dRep + (coeff[i] + coeff[j]) * dAtt / 2.0;
            
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

double GuptaPotential::calculateEnergyWithForces(const BinaryAlloyCluster& cluster, std::vector<double>& f) {
    calculateForces(cluster, f);
    return calculateEnergy(cluster);
}
