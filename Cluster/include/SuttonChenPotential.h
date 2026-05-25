#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "PotentialBase.h"

struct SuttonChenParameters {
    double epsilon;  // Energy unit (eV)
    double c;        // Embedding energy coefficient
    double a;        // Lattice constant
    double n;        // Pair potential exponent
    double m;        // Density exponent

    SuttonChenParameters()
        : epsilon(1.0), c(39.432), a(3.61), n(10.0), m(8.0) {
    }

    SuttonChenParameters(double epsilon_, double c_, double a_, double n_, double m_)
        : epsilon(epsilon_), c(c_), a(a_), n(n_), m(m_) {
    }
};

class SuttonChenPotential : public PotentialBase {
private:
    SuttonChenParameters paramsAA;
    SuttonChenParameters paramsBB;
    SuttonChenParameters paramsAB;

    double cutoff;

    void computeDistanceMatrix(const BinaryAlloyCluster& cluster, std::vector<double>& dist) const;
    SuttonChenParameters getParameters(const BinaryAlloyCluster& cluster,
        int i, int j) const;
    double calcEnergyWithDist(const BinaryAlloyCluster& cluster, const std::vector<double>& dist) const;
    void calcForcesWithDist(const BinaryAlloyCluster& cluster, std::vector<double>& f, const std::vector<double>& dist) const;

    // Pair potential V(r) = (a/r)^n
    double pairPotential(double r, const SuttonChenParameters& params) const;

    // Pair potential derivative dV/dr
    double pairPotentialDerivative(double r, const SuttonChenParameters& params) const;

    // Density function phi(r) = (a/r)^m
    double densityFunction(double r, const SuttonChenParameters& params) const;

    // Density function derivative dphi/dr
    double densityFunctionDerivative(double r, const SuttonChenParameters& params) const;

public:
    SuttonChenPotential();
    SuttonChenPotential(const std::string& elemA, const std::string& elemB);

    // Implement virtual functions from PotentialBase
    bool loadParameters(const std::string& filename) override;
    void setElements(const std::string& elemA, const std::string& elemB) override;
    double calculateEnergy(const BinaryAlloyCluster& cluster) override;
    double calculateEnergyWithForces(const BinaryAlloyCluster& cluster, std::vector<double>& f) override;
    std::string getPotentialType() const override { return "Sutton-Chen"; }
};
