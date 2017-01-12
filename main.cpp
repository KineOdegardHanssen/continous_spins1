#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <string>
#include "bond.h"
#include "site.h"
#include "lattice.h"
#include "montecarlo.h"

using namespace std;
using std::ofstream; using std::string;

int main()   // main. Monte Carlo steps here?
{
    bool DEBUG = true;
    if(DEBUG)    cout << "In main" << endl;

    // Input parameters
    int L = 10; // The program is going to be slower than before as we have a 3D lattice
    // bools to determine system
    bool isotropic    = true;
    bool sianisotropy = false;  // This one does not change its energy unless Dix, Diy and Diz are not all equal.
    bool magfield     = false;
    bool dm           = false;
    char type_lattice = 'F';
    //char latticetype = 'FH'; // FH: face-centered cubic,helical; C: cubic, helical; Q:quadratic, helical
    double beta = 2.5; // Just setting a beta.

    // Run parameters
    int eqsteps = 1000; // Number of steps in the equilibration procedure
    int mcsteps_inbin = 1000; // MCsteps per bin. Do I need bins?
    int no_of_bins = 100;     // The number of bins.

    if(DEBUG)    cout << "Parameters set" << endl;
    if(DEBUG)    cout << "beta = " << beta << endl;

    // Setting up the lattice with site parameters and interactions

    /*
    start_clock = clock();
    // Initializing instance of class Lattice
    Lattice mylattice = Lattice(L, isotropic, sianisotropy, magfield, dm);
    if(DEBUG)    cout << "Instance of class Lattice initialized" << endl;

    // Choosing type of lattice
    mylattice.fcc_helical_initialize();
    //mylattice.cubic_helical_initialize();
    //mylattice.quadratic_helical_initialize();
    cout << "in main again" << endl;

    end_clock = clock();
    double total_time_initialize_lattice = (end_clock - start_clock)/(double) CLOCKS_PER_SEC;
    cout << "Time to initialize Lattice: " << total_time_initialize_lattice  << endl;
    int no_of_neighbours = mylattice.no_of_neighbours;

    if(DEBUG)     cout << "Lattice set up" << endl;
    if(DEBUG)     cout << "Number of neighbours: " << no_of_neighbours << endl;
    */

    string filenamePrefix = "fcc10t10t10_iso1_beta2p5_compareMCclass";

    MonteCarlo mymc = MonteCarlo(L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, type_lattice);
    mymc.debugmode(true);
    mymc.latticetype(L, type_lattice);
    mymc.runmetropolis(beta, filenamePrefix);

}

