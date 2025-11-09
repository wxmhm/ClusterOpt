#pragma once
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "GuptaPotential.h"
#include "NELbfgs.h"

class LocalOptimizer {
private:
    GuptaPotential* potential;
    NELbfgs* lbfgs;
    int iterationCount;

public:
    LocalOptimizer(GuptaPotential* pot)
        : potential(pot), iterationCount(0) {
        lbfgs = new NELbfgs(pot);
    }

    ~LocalOptimizer() {
        delete lbfgs;
    }

    bool optimize(BinaryAlloyCluster& cluster) {
        double energy = lbfgs->optimize(cluster);
        iterationCount++;
        return true;
    }

    int getIterationCount() const { return iterationCount; }
};
