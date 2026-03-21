// ==================== src/NELbfgs.cpp ====================
#include "../include/NELbfgs.h"

NELbfgs::NELbfgs(PotentialBase* pot) : potential(pot) {
}

NELbfgs::~NELbfgs() {
}

lbfgsfloatval_t NELbfgs::evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
) {
    OptimizationContext* ctx = static_cast<OptimizationContext*>(instance);
    int N = n / 3;
    

    std::vector<double> coords(n);
    for (int i = 0; i < n; i++) {
        coords[i] = x[i];
    }
    ctx->cluster->getCoordinates() = coords;
    

    std::vector<double> forces(n);
    lbfgsfloatval_t energy = ctx->potential->calculateEnergyWithForces(*ctx->cluster, forces);
    

    for (int i = 0; i < n; i++) {
        g[i] = -forces[i];
    }
    
    return energy;
}

int NELbfgs::progress(
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
) {

    if (k % 100 == 0) {
        //printf("Iteration %d: fx = %f, gnorm = %f\n", k, fx, gnorm);
    }
    return 0;
}

double NELbfgs::optimize(BinaryAlloyCluster& cluster) {
    int N = cluster.getNumAtoms();
    int n = 3 * N;
    

    const auto& coords = cluster.getCoordinates();
    lbfgsfloatval_t *x = lbfgs_malloc(n);
    if (x == nullptr) {
        std::cerr << "ERROR: Failed to allocate memory for L-BFGS." << std::endl;
        return 0.0;
    }
    
    for (int i = 0; i < n; i++) {
        x[i] = coords[i];
    }
    

    OptimizationContext ctx;
    ctx.potential = potential;
    ctx.cluster = &cluster;
    ctx.numAtoms = N;
    

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    

    lbfgsfloatval_t fx;
    int ret = lbfgs(n, x, &fx, evaluate, progress, &ctx, &param);
    

    if (ret == LBFGS_SUCCESS) {
        //std::cout << "L-BFGS: Optimization converged successfully." << std::endl;
    } else if (ret == LBFGS_ALREADY_MINIMIZED) {
        //std::cout << "L-BFGS: Already at minimum." << std::endl;
    } else if (ret == LBFGS_STOP) {
        //std::cout << "L-BFGS: Stopped by user." << std::endl;
    } else {

        //std::cout << "L-BFGS: Terminated with code " << ret << std::endl;
    }
    

    std::vector<double> optimizedCoords(n);
    for (int i = 0; i < n; i++) {
        optimizedCoords[i] = x[i];
    }
    cluster.getCoordinates() = optimizedCoords;
    cluster.setEnergy(fx);
    

    lbfgs_free(x);
    
    return fx;
}

double NELbfgs::local(double *coords, double *force, int N) {

    BinaryAlloyCluster tempCluster(N, "A", "B");
    std::vector<double> coordVec(coords, coords + 3*N);
    tempCluster.getCoordinates() = coordVec;
    
    double energy = optimize(tempCluster);
    

    const auto& newCoords = tempCluster.getCoordinates();
    for (int i = 0; i < 3*N; i++) {
        coords[i] = newCoords[i];
    }
    
    return energy;
}
