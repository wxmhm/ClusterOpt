#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <deque>
#include <filesystem>
#include <thread>

namespace fs = std::filesystem;

// Constants
namespace Constants {
    constexpr double EPSILON = 1e-10;
}

// Random number utilities
class RandomGenerator {
private:
    static std::mt19937& getEngine() {
        static thread_local std::mt19937 engine(
            static_cast<unsigned int>(
                std::chrono::steady_clock::now().time_since_epoch().count() ^
                std::hash<std::thread::id>{}(std::this_thread::get_id())
            )
        );
        return engine;
    }

public:
    static double uniform(double min = 0.0, double max = 1.0) {
        std::uniform_real_distribution<double> dist(min, max);
        return dist(getEngine());
    }
    
    static int uniformInt(int min, int max) {
        std::uniform_int_distribution<int> dist(min, max);
        return dist(getEngine());
    }
    
    static double normal(double mean, double stddev) {
        std::normal_distribution<double> dist(mean, stddev);
        return dist(getEngine());
    }
    
    static double cauchy(double location, double scale) {
        std::cauchy_distribution<double> dist(location, scale);
        return dist(getEngine());
    }
    
    static std::vector<int> permutation(int n, int k) {
        std::vector<int> indices(n);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), getEngine());
        if (k < n) indices.resize(k);
        return indices;
    }
};
