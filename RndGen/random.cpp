/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "random.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

// Saves the current seed to a file
void Random::SaveSeed(const std::string& filename) const
{
    std::ofstream writeSeed(filename);
    if (writeSeed.is_open())
    {
        writeSeed << ((l_tot >> 36) & 4095) << " " << ((l_tot >> 24) & 4095) << " "
                  << ((l_tot >> 12) & 4095) << " " << (l_tot & 4095) << std::endl;
    }
    else
    {
        std::cerr << "PROBLEM: Unable to open " << filename << std::endl;
    }
}

// Generates a random number with a Gaussian distribution
double Random::Gauss(const double mean, const double sigma)
{
    double s = Rannyu();
    double t = Rannyu();
    double x = sqrt(-2. * std::log(1. - s)) * std::cos(2. * M_PI * t);
    return mean + x * sigma;
}

// Generates a random number in the range [min, max)
double Random::Rannyu(const double min, const double max)
{
    assert(max > min && "Max should be greater than min");
    return min + (max - min) * Rannyu();
}

// Generates a random number in the range [0, 1)
double Random::Rannyu(void)
{
    constexpr double twom48 = 1. / (1ull << 48);
    l_tot = l_tot * m_tot + n_tot;
    l_tot &= ((1ull << 48) - 1);
    double r = twom48 * l_tot;
    return r;
}

// Sets the seed for the RNG
void Random::SetRandom(int* s, int p1, int p2)
{
    l_tot = (static_cast<uint64_t>(s[0]) << (12 * 3)) + (static_cast<uint64_t>(s[1]) << (12 * 2)) + 
            (static_cast<uint64_t>(s[2]) << 12) + static_cast<uint64_t>(s[3]);

    uint16_t n1 = 0;
    uint16_t n2 = 0;
    uint16_t n3 = p1;
    uint16_t n4 = p2;

    n_tot = (static_cast<uint64_t>(n1) << (12 * 3)) + (static_cast<uint64_t>(n2) << (12 * 2)) + 
            (static_cast<uint64_t>(n3) << 12) + static_cast<uint64_t>(n4);
}

// Generates a random number with an Exponential distribution
double Random::Exponential(const double lambda)
{
    assert(lambda > 0 && "Lambda should be greater than zero");
    double y = Rannyu();
    double r = -std::log(1 - y) / lambda;
    return r;
}

// Generates a random number with a Lorentzian distribution
double Random::Lorentzian(const double x_0, const double gamma)
{
    assert(gamma > 0 && "Gamma should be greater than zero");
    double y = Rannyu();
    double r = gamma * std::tan(M_PI * (y - 0.5)) + x_0;
    return r;
}

// Generates a random number using the Accept-Reject method
double Random::AcceptReject(const double a, const double b, const double max, std::function<double(double)>& PDF)
{
    double x = 0, y = 0;
    do
    {
        x = Rannyu(a, b);
        y = Rannyu(0, max);
    } while (PDF(x) < y);
    return x;
}

// Generates a random number using an inverse cumulative distribution function
double Random::ExternalInvCum(std::function<double(double)>& ICDF)
{
    return ICDF(Rannyu());
}

// Generates a random integer in the range [0, 2**48)
uint64_t Random::Ranint()
{
    l_tot = l_tot * m_tot + n_tot;
    l_tot &= ((1ull << 48) - 1);
    uint64_t r = l_tot;
    return r;
}

// Generates a random integer in the range [min, max)
uint64_t Random::Ranint(const uint64_t min, const uint64_t max)
{
    assert((max > min) && "Supplied wrong value range");
    auto range = max - min;
    auto res = Ranint() % range + min;
    return res;
}

// Initializes the RNG
void Random::initialize()
{
    int seed[4]; // Array to store the seed
    int p1, p2; // Prime numbers for RNG initialization

    // Read prime numbers from file
    std::ifstream primes("RndGen/Primes");
    if (primes.is_open())
    {
        primes >> p1 >> p2;
    }
    else
    {
        std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
        exit(1); // Exit with error
    }
    primes.close();

    // Read the seed from file
    std::ifstream input("RndGen/seed.in");
    std::string property;
    if (input.is_open())
    {
        while (!input.eof())
        {
            input >> property;
            if (property == "RANDOMSEED")
            {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                SetRandom(seed, p1, p2);
            }
        }
        input.close();
    }
    else
    {
        std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;
        exit(1); // Exit with error
    }
}

// Calculates the standard deviation of the mean
double Random::error(double av, double av2, int n)
{
    return n == 0 ? 0 : sqrt((av2 - av * av) / n);
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
