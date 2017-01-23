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
#include "gaussiandeviate.h"  // Just to test against another random number generator.
//#include "printing.h"

using namespace std;
using std::ofstream; using std::string;

class MonteCarlo
{
public:
    // Printing part
    ofstream    allFile;
    ofstream    bigFile;
    ofstream    compareFile;
    string      filenamePrefix;

    int N, eqsteps, mcsteps_inbin, no_of_bins, no_of_neighbours;
    long int seed;
    double energy_old, acceptancerate;
    bool isotropic, sianisotropy, magfield, dm;
    bool notperiodic;
    bool printeveryMCstep;
    bool DEBUG, MAJORDEBUG;

    // Random generators
    std::default_random_engine generator_u;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_u; //(-1,1)

    std::default_random_engine generator_v;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_v; //(-1,1)

    std::default_random_engine generator_prob;                    // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_prob; // (0,1)

    // For index. This is given helical boundary conditions, then I only need one index
    std::default_random_engine generator_n;
    std::uniform_int_distribution<int> distribution_n; // (0, N-1)

    Lattice mylattice;
    //Printing print;

    MonteCarlo();
    MonteCarlo(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, bool isotropic, bool sianisotropy, bool magfield, bool periodic, bool printeveryMCstep, bool dm, char type_lattice, string filenamePrefix);

    void chooseprintfile(string filenamePrefix);
    // debugging
    void debugmode(bool on);
    void majordebugtrue();
    void debug1d2p();

    void latticetype(int L, char type_lattice);
    void initialize_energy();
    void reset_energy();
    void setrandomgenerators();

    // Standard Metropolis functions
    void runmetropolis(double beta); // Or should beta be a class variable?
    void mcstepf_metropolis(double beta); //, std::default_random_engine generator_u, std::default_random_engine generator_v, std::default_random_engine generator_n, std::default_random_engine generator_prob,  std::uniform_real_distribution<double> distribution_prob, std::uniform_real_distribution<double> distribution_u, std::uniform_real_distribution<double> distribution_v, std::uniform_int_distribution<int> distribution_n); //

    //
    void endsims();
};

#endif // MONTECARLO_H
