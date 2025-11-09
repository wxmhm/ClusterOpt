
#pragma once
#include "lbfgs.h"
#include "Common.h"
#include "BinaryAlloyCluster.h"
#include "GuptaPotential.h"

class NELbfgs {
private:
    GuptaPotential* potential;
    
    // 用于回调函数的上下文
    struct OptimizationContext {
        GuptaPotential* potential;
        BinaryAlloyCluster* cluster;
        int numAtoms;
        mutable std::vector<double> distanceMatrix;
    };
    
    // 静态回调函数
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
    NELbfgs(GuptaPotential* pot);
    ~NELbfgs();
    
    // 局部优化
    double optimize(BinaryAlloyCluster& cluster);
    
    // 兼容接口
    double local(double *coords, double *force, int N);
};
