#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "PotentialBase.h"

struct GuptaParameters {
    double r0;
    double A;
    double xi;
    double p;
    double q;
    
    GuptaParameters() : r0(2.5), A(0.1), xi(1.0), p(-10.0), q(-4.0) {}
    GuptaParameters(double r0_, double A_, double xi_, double p_, double q_)
        : r0(r0_), A(A_), xi(xi_), p(p_), q(q_) {}
};

class GuptaPotential : public PotentialBase {
private:
    GuptaParameters paramsAA;
    GuptaParameters paramsBB;
    GuptaParameters paramsAB;
    
    mutable std::vector<double> distanceMatrix;
    mutable std::vector<double> forces;
    
    void computeDistanceMatrix(const BinaryAlloyCluster& cluster) const;
    GuptaParameters getParameters(const BinaryAlloyCluster& cluster, int i, int j) const;
    
public:
    GuptaPotential();
    GuptaPotential(const std::string& elemA, const std::string& elemB);
    
    // Implement virtual functions from PotentialBase
    bool loadParameters(const std::string& filename) override;
    void setElements(const std::string& elemA, const std::string& elemB) override;
    double calculateEnergy(const BinaryAlloyCluster& cluster) override;
    void calculateForces(const BinaryAlloyCluster& cluster, std::vector<double>& f) override;
    double calculateEnergyWithForces(const BinaryAlloyCluster& cluster, std::vector<double>& f) override;
    std::string getPotentialType() const override { return "Gupta"; }
    
    // Gupta-specific methods
    void setParameters(const GuptaParameters& aa, const GuptaParameters& bb, const GuptaParameters& ab);
    const GuptaParameters& getParamsAA() const { return paramsAA; }
    const GuptaParameters& getParamsBB() const { return paramsBB; }
    const GuptaParameters& getParamsAB() const { return paramsAB; }
};
