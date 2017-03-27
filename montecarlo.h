#ifndef MONTECARLO_H
#define MONTECARLO_H
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <cmath>
#include <complex>
#include "lattice.h"
#include "site.h"    // Are these really neccessary... ?
#include "bond.h"
#include "gaussiandeviate.h"  // Just to test against another random number generator.


// For correlation function, need to access FFTW. Is this notation really neccessary? Could I just include?
extern "C"
{
#include "fftw3.h"
}


using namespace std;
using std::ofstream; using std::string;

class MonteCarlo
{
public:
    // Printing part
    ofstream    allFile;
    ofstream    bigFile;
    ofstream    compareFile;
    ofstream    spcorFile;
    ofstream    ftspcorFile;
    ofstream    randomtestFile;
    ofstream    qFile;
    string      filenamePrefix;

    int N, eqsteps, mcsteps_inbin, no_of_bins, no_of_neighbours;
    long int seed1, seed2, testseed;
    double energy_old, acceptancerate;
    bool isotropic, sianisotropy, magfield, dm;
    bool notperiodic;
    bool printeveryMCstep, calculatespincorrelationfunction, randomtest;
    bool DEBUG, MAJORDEBUG;

    fftw_plan p;
    fftw_plan pinv;

    // Random generators
    // A lot of these I don't use anymore
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
    MonteCarlo(int L1, int eqsteps, int mcsteps_inbin, int no_of_bins, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction,  char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
    MonteCarlo(int L1, int L2, int L3, int eqsteps, int mcsteps_inbin, int no_of_bins, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction,  char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);


    // Debugging/testing functions
    void debugmode(bool on);
    void majordebugtrue();
    void debug1d2p();
    void testFFTW();
    void compareFFTW_withmanual(double beta);
    double check_the_energy(); // Function that finds the energy by going through every spin
    void test_couplings_strengths();

    // Other initialization procedures
    void initialize_energy();
    void reset_energy();

    // Function making a plan for the FFT
    void giveplanforFFT(vector<double> &r, vector<complex<double> > &q);
    void giveplanforFFT_inverse(vector<double>& rout, vector<complex<double> >& q);


    // Standard Metropolis functions
    void runmetropolis(double beta);
    void mcstepf_metropolis(double beta);

    // Printing functions
    void writeallqstofile();
    // Closing the output files
    void endsims();
};

#endif // MONTECARLO_H
