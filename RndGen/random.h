/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random_Sim__
#define __Random_Sim__

#include <cstdint>
#include <functional>
#include <string>

// The Random class provides various methods for generating random numbers
class Random
{
private:
    const uint64_t m_tot{34522712143931ull}; // Constant value for the RNG
    uint64_t l_tot, n_tot; // Internal state variables for the RNG

public:
    // Default constructor
    Random() = default;
    // Default destructor
    ~Random() = default;

    // Sets the seed for the RNG
    void SetRandom(int* seeds, int p1, int p2);

    // Saves the current seed to a file
    void SaveSeed(const std::string& filename = "seed.out") const;

    // Generates a random number in the range [0, 1)
    double Rannyu(void);

    // Generates a random number in the range [min, max)
    double Rannyu(const double min, const double max);

    // Generates a random integer in the range [0, 2**48)
    uint64_t Ranint();

    // Generates a random integer in the range [min, max)
    uint64_t Ranint(const uint64_t min, const uint64_t max);

    // Generates a random number with a Gaussian distribution
    double Gauss(const double mean, const double sigma);

    // Generates a random number with an Exponential distribution
    double Exponential(const double lambda);

    // Generates a random number with a Lorentzian distribution
    double Lorentzian(const double x_0, const double gamma);

    // Generates a random number using the Accept-Reject method
    double AcceptReject(const double a, const double b, const double max, std::function<double(double)>& PDF);

    // Generates a random number using an inverse cumulative distribution function
    double ExternalInvCum(std::function<double(double)>& ICDF);

    // Initializes the RNG
    void initialize();

    // Calculates the standard deviation of the mean
    double error(double av, double av2, int n);
};

#endif // __Random_Sim__


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
