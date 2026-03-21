#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "PotentialBase.h"

struct FinnisSinclairParameters {
    double cutoff;
    double c0, c1, c2;
    double d;
    double c;
    double beta;  // Cubic term coefficient for density function
    
    FinnisSinclairParameters() 
        : cutoff(5.0), c0(1.0), c1(-0.5), c2(0.1), d(4.5), c(1.0), beta(0.0) {}
    
    FinnisSinclairParameters(double cutoff_, double c0_, double c1_, 
                            double c2_, double d_, double c_, double beta_)
        : cutoff(cutoff_), c0(c0_), c1(c1_), c2(c2_), d(d_), c(c_), beta(beta_) {}
};

class FinnisSinclairPotential : public PotentialBase {
private:
    FinnisSinclairParameters paramsAA;
    FinnisSinclairParameters paramsBB;
    FinnisSinclairParameters paramsAB;
    
    mutable std::vector<double> distanceMatrix;
    
    void computeDistanceMatrix(const BinaryAlloyCluster& cluster) const;
    FinnisSinclairParameters getParameters(const BinaryAlloyCluster& cluster, 
                                           int i, int j) const;
    

    double pairPotential(double r, const FinnisSinclairParameters& params) const;
    

    double pairPotentialDerivative(double r, const FinnisSinclairParameters& params) const;
    

    double densityFunction(double r, const FinnisSinclairParameters& params) const;
    

    double densityFunctionDerivative(double r, const FinnisSinclairParameters& params) const;
    
public:
    FinnisSinclairPotential();
    FinnisSinclairPotential(const std::string& elemA, const std::string& elemB);
    
    // Implement virtual functions from PotentialBase
    bool loadParameters(const std::string& filename) override;
    void setElements(const std::string& elemA, const std::string& elemB) override;
    double calculateEnergy(const BinaryAlloyCluster& cluster) override;
    void calculateForces(const BinaryAlloyCluster& cluster, std::vector<double>& f) override;
    double calculateEnergyWithForces(const BinaryAlloyCluster& cluster, std::vector<double>& f) override;
    std::string getPotentialType() const override { return "Finnis-Sinclair"; }
    
    // FinnisSinclair-specific methods
    void setParameters(const FinnisSinclairParameters& aa, 
                      const FinnisSinclairParameters& bb, 
                      const FinnisSinclairParameters& ab);
    const FinnisSinclairParameters& getParamsAA() const { return paramsAA; }
    const FinnisSinclairParameters& getParamsBB() const { return paramsBB; }
    const FinnisSinclairParameters& getParamsAB() const { return paramsAB; }
};
