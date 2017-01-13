#ifndef MONTECARLO_H
#define MONTECARLO_H
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <cmath>
#include "lattice.h"
#include "site.h"    // Are these really neccessary... ?
#include "bond.h"
#include "printing.h"

using namespace std;
using std::ofstream; using std::string;

class MonteCarlo
{
public:
    MonteCarlo(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, bool isotropic, bool sianisotropy, bool magfield, bool dm, char type_lattice, string filenamePrefix);

    int N, eqsteps, mcsteps_inbin, no_of_bins, no_of_neighbours;
    double energy_old, acceptancerate;
    bool isotropic, sianisotropy, magfield, dm;
    bool DEBUG;

    Lattice mylattice;
    Printing print;

    void debugmode(bool on);

    void latticetype(int L, char type_lattice);
    void initialize_energy();
    void reset_energy();

    // Standard Metropolis functions
    void runmetropolis(double beta); // Or should beta be a class variable?
    void mcstepf_metropolis(double beta, std::default_random_engine generator_u, std::default_random_engine generator_v, std::default_random_engine generator_n, std::default_random_engine generator_prob,  std::uniform_real_distribution<double> distribution_prob, std::uniform_real_distribution<double> distribution_u, std::uniform_real_distribution<double> distribution_v, std::uniform_int_distribution<int> distribution_n); //

};

#endif // MONTECARLO_H
