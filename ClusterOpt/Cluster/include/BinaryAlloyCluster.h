#pragma once
#include "Common.h"

class BinaryAlloyCluster {
private:
    std::vector<double> coordinates;
    std::vector<int> atomTypes;
    std::string elementA;
    std::string elementB;
    double energy;
    int numAtoms;
    
public:
    BinaryAlloyCluster(int n, const std::string& elemA = "A", const std::string& elemB = "B");
    BinaryAlloyCluster(int nA, int nB, const std::string& elemA, const std::string& elemB);
    
    // Element access
    void setElements(const std::string& elemA, const std::string& elemB);
    std::string getElementA() const { return elementA; }
    std::string getElementB() const { return elementB; }
    std::string getElementSymbol(int atomIndex) const;
    std::string getCompositionString() const;
    
    // Coordinate access
    double* data() { return coordinates.data(); }
    const double* data() const { return coordinates.data(); }
    std::vector<double>& getCoordinates() { return coordinates; }
    const std::vector<double>& getCoordinates() const { return coordinates; }
    
    // Atom access
    std::array<double, 3> getAtomPosition(int i) const;
    void setAtomPosition(int i, double x, double y, double z);
    
    // Properties
    int getNumAtoms() const { return numAtoms; }
    int getAtomType(int i) const { return atomTypes[i]; }
    double getEnergy() const { return energy; }
    void setEnergy(double e) { energy = e; }
    int getNumElementA() const;
    int getNumElementB() const;
    
    // Operations
    void randomInitialize(double boxSize);
    int countStructureFiles();
    bool loadStructureInitialize(int x, int expectedAtoms);
    void centerAtOrigin();
    double getRadius() const;
    BinaryAlloyCluster copy() const { return *this; }
};
