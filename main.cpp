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

using namespace std;
using std::ofstream; using std::string;

void one_run(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, char type_lattice, string filenamePrefix);
void run_for_several_betas(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, double betamin, double betamax, bool isotropic, bool sianisotropy, bool magfield, bool dm, char type_lattice, string filenamePrefix);
void test_betagenerator(int beta_n, int betamin, int betamax);

int main()   // main. Monte Carlo steps here?
{
    bool DEBUG = true;
    if(DEBUG)    cout << "In main" << endl;

    // Input parameters
    int L = 2; // The program is going to be slower than before as we have a 3D lattice
    // bools to determine system
    bool isotropic    = true;
    bool sianisotropy = false;  // This one does not change its energy unless Dix, Diy and Diz are not all equal.
    bool magfield     = false;
    bool dm           = false;
    char type_lattice = 'P';
    //char latticetype = 'FH'; // FH: face-centered cubic,helical; C: cubic, helical; Q:quadratic, helical
    double beta = 2.5; // Just setting a beta.

    // Run parameters
    int eqsteps = 1000; // Number of steps in the equilibration procedure
    int mcsteps_inbin = 1000; // MCsteps per bin. Do I need bins?
    int no_of_bins = 100;     // The number of bins.

    string filenamePrefix = "test";
    //string filenamePrefix = "chain2_periodic_iso1_beta0to4";

    //test_betagenerator(10, 0, 4);
    int beta_n = 100;
    double betamin = 0.01;
    double betamax = 4;

    run_for_several_betas(L, eqsteps, mcsteps_inbin, no_of_bins, beta_n, betamin, betamax, isotropic, sianisotropy, magfield, dm, type_lattice, filenamePrefix);

    //one_run(L, eqsteps, mcsteps_inbin, no_of_bins, beta, isotropic, sianisotropy, magfield, dm, type_lattice, filenamePrefix);
}

void one_run(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, char type_lattice, string filenamePrefix)
{
    // Initializing Monte Carlo
    MonteCarlo mymc(L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, type_lattice, filenamePrefix);
    mymc.debugmode(true);
    // Run Metropolis algorithm
    mymc.runmetropolis(beta);
    mymc.endsims();
}

void run_for_several_betas(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, double betamin, double betamax, bool isotropic, bool sianisotropy, bool magfield, bool dm, char type_lattice, string filenamePrefix)
{
    // Initializing Monte Carlo
    MonteCarlo mymc(L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, type_lattice, filenamePrefix);
    mymc.debugmode(true);
    //mymc.majordebugtrue();
    //mymc.latticetype(L, type_lattice);  // Runs for chain, quadratic and cubic lattice, not for fcc.

    //int beta_n = 100; // Or something. Could have this as input
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

