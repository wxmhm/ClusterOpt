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
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <atomic>

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

// Simple thread pool with submit-and-wait pattern
class ThreadPool {
private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queueMutex;
    std::condition_variable cv;
    std::condition_variable doneCv;
    std::atomic<int> tasksRemaining{0};
    bool stop = false;

public:
    explicit ThreadPool(int numThreads) {
        for (int i = 0; i < numThreads; ++i) {
            workers.emplace_back([this]() {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(queueMutex);
                        cv.wait(lock, [this]() { return stop || !tasks.empty(); });
                        if (stop && tasks.empty()) return;
                        task = std::move(tasks.front());
                        tasks.pop();
                    }
                    task();
                    tasksRemaining--;
                    doneCv.notify_one();
                }
            });
        }
    }

    void submit(std::function<void()> task) {
        tasksRemaining++;
        {
            std::lock_guard<std::mutex> lock(queueMutex);
            tasks.push(std::move(task));
        }
        cv.notify_one();
    }

    void waitAll() {
        std::unique_lock<std::mutex> lock(queueMutex);
        doneCv.wait(lock, [this]() { return tasksRemaining == 0; });
    }

    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lock(queueMutex);
            stop = true;
        }
        cv.notify_all();
        for (auto& w : workers) {
            if (w.joinable()) w.join();
        }
    }
};
