#pragma once

// Windows-specific includes
#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <windows.h>
#include <direct.h>  // for _mkdir
#define mkdir(path) _mkdir(path)
#else
#include <sys/stat.h>
#endif

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
#include <map>
#include <functional>
#include <filesystem>

namespace fs = std::filesystem;

// Constants
namespace Constants {
    constexpr double PI = 3.14159265358979323846;
    constexpr double EPSILON = 1e-10;
    constexpr double BOLTZMANN = 8.617333262145e-5; // eV/K
}

// Element information system
class ElementInfo {
private:
    static std::map<std::string, double> atomicMasses;
    static std::map<std::string, int> atomicNumbers;
    
public:
    static void initialize() {
        if (!atomicMasses.empty()) return;
        
        atomicMasses = {
            {"H", 1.008}, {"Li", 6.941}, {"C", 12.011}, {"N", 14.007}, {"O", 15.999},
            {"Al", 26.982}, {"Si", 28.086}, {"Ti", 47.867}, {"V", 50.942}, {"Cr", 51.996},
            {"Mn", 54.938}, {"Fe", 55.845}, {"Co", 58.933}, {"Ni", 58.693}, {"Cu", 63.546},
            {"Zn", 65.409}, {"Ga", 69.723}, {"Ge", 72.64}, {"Mo", 95.94}, {"Ru", 101.07},
            {"Rh", 102.906}, {"Pd", 106.42}, {"Ag", 107.868}, {"Cd", 112.411}, {"In", 114.818},
            {"Sn", 118.71}, {"W", 183.84}, {"Re", 186.207}, {"Os", 190.23}, {"Ir", 192.217},
            {"Pt", 195.084}, {"Au", 196.967}, {"Hg", 200.59}, {"Pb", 207.2}, {"Bi", 208.98}
        };
        
        atomicNumbers = {
            {"H", 1}, {"Li", 3}, {"C", 6}, {"N", 7}, {"O", 8},
            {"Al", 13}, {"Si", 14}, {"Ti", 22}, {"V", 23}, {"Cr", 24},
            {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29},
            {"Zn", 30}, {"Ga", 31}, {"Ge", 32}, {"Mo", 42}, {"Ru", 44},
            {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49},
            {"Sn", 50}, {"W", 74}, {"Re", 75}, {"Os", 76}, {"Ir", 77},
            {"Pt", 78}, {"Au", 79}, {"Hg", 80}, {"Pb", 82}, {"Bi", 83}
        };
    }
    
    static double getMass(const std::string& element) {
        auto it = atomicMasses.find(element);
        return (it != atomicMasses.end()) ? it->second : 1.0;
    }
    
    static int getAtomicNumber(const std::string& element) {
        auto it = atomicNumbers.find(element);
        return (it != atomicNumbers.end()) ? it->second : 0;
    }
};

// Static member definitions
//std::map<std::string, double> ElementInfo::atomicMasses;
//std::map<std::string, int> ElementInfo::atomicNumbers;

// Random number utilities
class RandomGenerator {
private:
    static std::mt19937& getEngine() {
        static thread_local std::mt19937 engine(
            static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count())
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
