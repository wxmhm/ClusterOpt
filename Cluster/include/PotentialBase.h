#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"

/**
 * @brief Abstract base class for all potential energy functions
 * 
 * This provides a common interface for different potential types
 * (Gupta, Finnis-Sinclair, Sutton-Chen, etc.)
 */
class PotentialBase {
protected:
    std::string elementA;
    std::string elementB;

public:
    // Default constructor
    PotentialBase() : elementA("A"), elementB("B") {}
    
    // Constructor with elements
    PotentialBase(const std::string& elemA, const std::string& elemB) 
        : elementA(elemA), elementB(elemB) {}
    
    virtual ~PotentialBase() = default;

    // Pure virtual functions that must be implemented by derived classes
    virtual bool loadParameters(const std::string& filename) = 0;
    virtual void setElements(const std::string& elemA, const std::string& elemB) = 0;
    virtual double calculateEnergy(const BinaryAlloyCluster& cluster) = 0;
    virtual void calculateForces(const BinaryAlloyCluster& cluster, std::vector<double>& f) = 0;
    virtual double calculateEnergyWithForces(const BinaryAlloyCluster& cluster, std::vector<double>& f) = 0;

    // Getter methods
    virtual std::string getElementA() const { return elementA; }
    virtual std::string getElementB() const { return elementB; }
    
    // Virtual method to get potential type name
    virtual std::string getPotentialType() const = 0;
};
