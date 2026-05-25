
#pragma once
#include "lbfgs.h"
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "PotentialBase.h"

class NELbfgs {
private:
    PotentialBase* potential;
    
    // Context for callback functions
    struct OptimizationContext {
        PotentialBase* potential;
        BinaryAlloyCluster* cluster;
        int numAtoms;
    };
    
    // Static callback functions
    static lbfgsfloatval_t evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    );
    
    static int progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
    );

public:
    NELbfgs(PotentialBase* pot);
    ~NELbfgs();
    
    // Local optimization
    double optimize(BinaryAlloyCluster& cluster);
};
