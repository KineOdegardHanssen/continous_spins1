#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <string>
#include "bond.h"
#include "site.h"
#include "lattice.h"
#include "printing.h"
#include "montecarlo.h"

//#include "fftw3.h"  // Only use in MonteCarlo?
// Do not think I need these paths:
//INCLUDEPATH += "/usr/share/doc/libfftw3-3"
//INCLUDEPATH += "/home/ubu/Downloads/fftw-3.3.6-pl1"
// Have    extern "C"   somewhere

using namespace std;
using std::ofstream; using std::string;

void one_run(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction, char type_lattice, string filenamePrefix);
void run_for_several_betas(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, double betamin, double betamax, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, char type_lattice, string filenamePrefix);
void test_betagenerator(int beta_n, int betamin, int betamax);
void test_fftw(int L);

int main()   // main. Monte Carlo steps here?
{
    bool DEBUG = true;
    if(DEBUG)    cout << "In main" << endl;

    // Input parameters
    int L = 4; // The program is going to be slower than before as we have a 3D lattice

    // bools to determine system type
    bool isotropic    = true;
    bool sianisotropy = false;  // This one does not change its energy unless Dix, Diy and Diz are not all equal.
    bool magfield     = false;
    bool dm           = false;

    // Bool to determine periodicity
    bool periodic     = true; // To determine whether we have periodic boundary conditions or not

    // Selecting the lattice type
    char type_lattice = 'O';   // F: face-centered cubic; C: cubic; Q:quadratic; O: chain;
    // If periodic is false, that means we get a grid with open boundary conditions. Currently,
    // that is only implemented for the chain.


    // Bools to determine printing
    bool printeveryMCstep = false;
    bool calculatespincorrelationfunction = true; // test this

    // A beta value for one run
    double beta = 10000.0;

    // Run parameters
    int eqsteps = 1000; // Number of steps in the equilibration procedure
    int mcsteps_inbin = 1000; // MCsteps per bin. Do I need bins?
    int no_of_bins = 100;     // The number of bins.

    // Filenames (choose one to use or change slightly)
    string filenamePrefix = "printtestyo";
    //string filenamePrefix = "test_periodicchain_2p";
    //string filenamePrefix = "bigtest_periodicchain_2p_beta0p00001and4000_10000eqsteps_10000mcsteps_1000bins";
    //string filenamePrefix = "bigtest_openchain_5p_beta1em5and4000_10000eqsteps_10000mcsteps_100bins";
    //string filenamePrefix = "bigtest2";
    //string filenamePrefix = "chain2_periodic_iso1_beta0to4";

    //test_betagenerator(10, 0, 4);
    // Input parameters specifically for run_for_several_betas
    int beta_n = 2;
    double betamin = 1e-5;
    double betamax = 4000;

    // By default, run_for_several_betas do not calculate the correlation function
    //run_for_several_betas(L, eqsteps, mcsteps_inbin, no_of_bins, beta_n, betamin, betamax, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, type_lattice, filenamePrefix);

    test_fftw(L);

    //one_run(L, eqsteps, mcsteps_inbin, no_of_bins, beta, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix);
}

void one_run(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction, char type_lattice, string filenamePrefix)
{
    // Initializing Monte Carlo
    MonteCarlo mymc(L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix);
    mymc.debugmode(true);
    // Run Metropolis algorithm
    mymc.runmetropolis(beta);
    mymc.endsims();
}

void run_for_several_betas(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, double betamin, double betamax, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, char type_lattice, string filenamePrefix)
{
    bool calculatespincorrelationfunction = false;
    // Initializing Monte Carlo
    MonteCarlo mymc(L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix);
    //mymc.debugmode(true);
    //mymc.majordebugtrue();

    vector<double> betas = vector<double>(beta_n);
    double deltabeta = (betamax-betamin)/(beta_n-1.0);

    for(int i=0; i<beta_n; i++)
    {
        betas[i] = betamin + deltabeta*i;
        mymc.runmetropolis(betas[i]);
        mymc.reset_energy();
    }
    mymc.endsims();
}

void test_betagenerator(int beta_n,  int betamin, int betamax)
{
    vector<double> betas = vector<double>(beta_n);
    double deltabeta = (betamax-betamin)/(beta_n-1.0);

    cout << "betamin = " << betamin << "; betamax = " << betamax << "; number of values: " << beta_n << endl;
    cout << "deltabeta = " << deltabeta << endl;

    for(int i=0; i<beta_n; i++)
    {
        betas[i] = betamin + deltabeta*i;
        cout << "i = " << i << "; betas[i] = " << betas[i] << endl;
    }
}


void test_fftw(int L)
{
    // All these parameters are irrelevant for this simple test
    int eqsteps = 1000; int mcsteps_inbin = 1000; int no_of_bins= 100;
    bool isotropic = true; bool sianisotropy = false; bool magfield = false; bool dm = false;
    bool periodic = true;

    // These we want to turn off
    bool printeveryMCstep = false;
    bool calculatespincorrelationfunction = false;

    // Set these
    char type_lattice = 'O';  // Looking at a chain (at least for now)
    string filenamePrefix = "discard"; // Want to know that this file is unimportant

    // Initializing Monte Carlo
    MonteCarlo mymc(L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix);
    mymc.testFFTW(); // Test FFTW
    mymc.endsims();  // Close the file
}
