// ==================== src/NELbfgs.cpp ====================
#include "../include/NELbfgs.h"

NELbfgs::NELbfgs(GuptaPotential* pot) : potential(pot) {
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
    
    // 更新cluster坐标
    std::vector<double> coords(n);
    for (int i = 0; i < n; i++) {
        coords[i] = x[i];
    }
    ctx->cluster->getCoordinates() = coords;
    
    // 计算能量和力
    std::vector<double> forces(n);
    lbfgsfloatval_t energy = ctx->potential->calculateEnergyWithForces(*ctx->cluster, forces);
    
    // L-BFGS需要负梯度（力已经是负梯度，所以需要取负）
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
    // 可选：每100次迭代打印一次进度
    if (k % 100 == 0) {
        //printf("Iteration %d: fx = %f, gnorm = %f\n", k, fx, gnorm);
    }
    return 0;
}

double NELbfgs::optimize(BinaryAlloyCluster& cluster) {
    int N = cluster.getNumAtoms();
    int n = 3 * N;
    
    // 准备初始坐标
    const auto& coords = cluster.getCoordinates();
    lbfgsfloatval_t *x = lbfgs_malloc(n);
    if (x == nullptr) {
        std::cerr << "ERROR: Failed to allocate memory for L-BFGS." << std::endl;
        return 0.0;
    }
    
    for (int i = 0; i < n; i++) {
        x[i] = coords[i];
    }
    
    // 准备优化上下文
    OptimizationContext ctx;
    ctx.potential = potential;
    ctx.cluster = &cluster;
    ctx.numAtoms = N;
    
    // 设置L-BFGS参数
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    //param.m = 6;
    //param.epsilon = 1e-3;
    //param.max_iterations = 1000;
    //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
    /*param.max_iterations = 3000;
    param.epsilon = 1e-6;         // 梯度收敛阈值
    param.delta = 1e-6;           // 函数值收敛阈值
    param.past = 1;               // 用于判断收敛的历史窗口
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
    param.max_linesearch = 40;
    param.ftol = 1e-4;            // Armijo条件参数
    param.gtol = 0.9;             // Wolfe条件参数
    param.xtol = 1e-16;           // 机器精度
    param.orthantwise_c = 0;      // 不使用L1正则化
    param.orthantwise_start = 0;
    param.orthantwise_end = -1;
    */

    // 运行优化
    lbfgsfloatval_t fx;
    int ret = lbfgs(n, x, &fx, evaluate, progress, &ctx, &param);
    
    // 处理返回状态
    if (ret == LBFGS_SUCCESS) {
        //std::cout << "L-BFGS: Optimization converged successfully." << std::endl;
    } else if (ret == LBFGS_ALREADY_MINIMIZED) {
        //std::cout << "L-BFGS: Already at minimum." << std::endl;
    } else if (ret == LBFGS_STOP) {
        //std::cout << "L-BFGS: Stopped by user." << std::endl;
    } else {
        // 其他情况（如达到最大迭代次数）也是可接受的
        //std::cout << "L-BFGS: Terminated with code " << ret << std::endl;
    }
    
    // 更新cluster坐标
    std::vector<double> optimizedCoords(n);
    for (int i = 0; i < n; i++) {
        optimizedCoords[i] = x[i];
    }
    cluster.getCoordinates() = optimizedCoords;
    cluster.setEnergy(fx);
    
    // 释放内存
    lbfgs_free(x);
    
    return fx;
}

double NELbfgs::local(double *coords, double *force, int N) {
    // 兼容旧接口
    BinaryAlloyCluster tempCluster(N, "A", "B");
    std::vector<double> coordVec(coords, coords + 3*N);
    tempCluster.getCoordinates() = coordVec;
    
    double energy = optimize(tempCluster);
    
    // 复制回坐标
    const auto& newCoords = tempCluster.getCoordinates();
    for (int i = 0; i < 3*N; i++) {
        coords[i] = newCoords[i];
    }
    
    return energy;
}
