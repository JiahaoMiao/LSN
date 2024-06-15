#include <iostream>
#include <chrono>
#include "random.h"          // Include the new random number generator
#include "random_old.h"      // Include the old random number generator

using namespace std::chrono;

// Function to compare Rannyu() outputs of two RNGs
bool compareRannyu(Random& rng_new, Random_old& rng_old, size_t iterations)
{
    for (size_t i = 0; i < iterations; ++i)
    {
        double new_val = rng_new.Rannyu();
        double old_val = rng_old.Rannyu();
        if (new_val != old_val)
        {
            std::cout << "Difference found at iteration " << i << ": new_val = " 
                      << new_val << ", old_val = " << old_val << std::endl;
            return false;
        }
    }
    return true;
}

int main()
{
    size_t iterations = 1000000000; // Number of iterations for benchmarking

    // Initialize new RNG
    Random rng_new;
    rng_new.initialize();
    
    // Initialize old RNG
    Random_old rng_old;
    rng_old.initialize();
    
    // Benchmark new RNG
    auto start = high_resolution_clock::now();
    
    for (size_t i = 0; i < iterations; ++i)
    {
        rng_new.Rannyu();
    }
    
    auto end = high_resolution_clock::now();
    duration<double> duration = end - start;
    std::cout << "Time taken by new Rannyu() for " << iterations << " iterations: " 
              << duration.count() << " seconds" << std::endl;
    
    // Benchmark old RNG
    start = high_resolution_clock::now();
    
    for (size_t i = 0; i < iterations; ++i)
    {
        rng_old.Rannyu();
    }
    
    end = high_resolution_clock::now();
    duration = end - start;
    std::cout << "Time taken by old Rannyu() for " << iterations << " iterations: " 
              << duration.count() << " seconds" << std::endl;
    

    // Compare RNGs
    size_t compare_iterations = 1000000; // Number of iterations for comparison
    bool same_sequence = compareRannyu(rng_new, rng_old, compare_iterations);
    std::cout << "The RNGs produce the same sequence for the first " << compare_iterations 
              << " iterations: " << (same_sequence ? "Yes" : "No") << std::endl;

    return 0;
}
