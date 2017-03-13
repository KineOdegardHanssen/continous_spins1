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
//#include "gaussiandeviate.h"

using namespace std;
using std::ofstream; using std::string;

void one_run(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void run_for_several_betas(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, double betamin, double betamax, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void run_for_betasgiven(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, char type_lattice, string filenamePrefix, vector<double> betas, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void run_for_betasgiven2(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, char type_lattice, string filenamePrefix, vector<double> betas, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void run_for_betasgiven_diffdirs(int L1, int L2, int L3, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, char type_lattice, string filenamePrefix, vector<double> betas, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);

// Special functions
void extract_yline(int L, string ylinefilenamePrefix); // Extracting the lattice sites along (0,y,0) in the fcc


// Test functions
void test_betagenerator(int beta_n, int betamin, int betamax);
void test_fftw(int L, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void test_fcc_extended(int L, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void test_fcc_extended_diffdims(int L1, int L2, int L3, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void test_dm();
double dm_oneneighbour_fortest(double x1, double y1, double z1, double x2, double y2, double z2, double Dx, double Dy, double Dz);
double J_oneneighbour_fortest(double x1, double y1, double z1, double x2, double y2, double z2, double J);

int main()
{
    bool DEBUG = true;
    if(DEBUG)    cout << "In main" << endl;

    // Input parameters
    int L = 2; // The program is going to be slow if we run for many particles on a 3D lattice

    int L1 = 2;
    int L2 = 2;
    int L3 = 2;

    // bools to determine system type
    bool isotropic    = true;
    bool sianisotropy = false;  // This one does not change its energy unless Dix, Diy and Diz are not all equal.
    bool magfield     = false;
    bool dm           = false;

    // Bool to determine periodicity
    bool periodic     = false; // To determine whether we have periodic boundary conditions or not

    // Selecting the lattice type
    // F: face-centered cubic (fcc); E: fcc with different directions C: cubic; Q:quadratic; O: chain;
    char type_lattice = 'O';
    // If periodic is false, that means we get a grid with open boundary conditions. Currently,
    // that is only implemented for the chain.

    // Setting the strengths of the Hamiltonian terms
    // Magnetic field terms
    double hx = 1;    double hy = 1;    double hz = 1;
    // Single-ion anisotropy terms
    double Dix = 2;    double Diy = 2;    double Diz = 0;
    // Heisenberg term
    double J = -1;
    // Heisenberg terms with varying strengths (for fcc_initialize_extended E)
    double Jy  = 0;    double Jz  = 0;
    double Jxy = -1;    double Jxz = -1;    double Jyz = -1;
    // DM terms
    double Dx = 0;     double Dy = 0;    double Dz = 1;

    vector<double> sitestrengthsin = vector<double>(6);
    sitestrengthsin[0] = hx;    sitestrengthsin[1] = hy;    sitestrengthsin[2] = hz;
    sitestrengthsin[3] = Dix;   sitestrengthsin[4] = Diy;   sitestrengthsin[5] = Diz;

    vector<double> heisenbergin = vector<double>(6);
    heisenbergin[0] = J;
    heisenbergin[1] = Jy;       heisenbergin[2] = Jz;
    heisenbergin[3] = Jxy;      heisenbergin[4] = Jxz;      heisenbergin[5] = Jyz;

    vector<double> dm_in = vector<double>(3);
    dm_in[0] = Dx;              dm_in[1] = Dy;              dm_in[2] = Dz;

    // Bools to determine printing
    bool printeveryMCstep = false;
    bool calculatespincorrelationfunction = true;

    // A beta value for one run
    double beta = 1.0;

    // Run parameters
    int eqsteps = 10000; // Number of steps in the equilibration procedure
    int mcsteps_inbin = 100000; //100000; // MCsteps per bin.
    int no_of_bins = 100;     // The number of bins.

    // Filenames (choose one to use or change slightly)
    //string filenamePrefix = "0spincorrtest_fcc10x10x10_Dix2Diy2_DxDyDz1_Jyz1_beta0p5_eqsteps50000_mcsteps_inbin_1000_no_of_bins1000";
    //string filenamePrefix = "0fcc_2x2x2_thecomparebetas_eqsteps1000_mcstepsinbin1000_bins100_Jyz1_DixDiy2Diz0";
    //string filenamePrefix = "2pchainopen_J1Dz100_mcsteps_inbin100000_no_of_bins_100";
    //string filenamePrefix = "bigtest_periodicchain_2p_beta0p00001and4000_10000eqsteps_10000mcsteps_1000bins";
    //string filenamePrefix = "bigtest_openchain_5p_beta1em5and4000_10000eqsteps_10000mcsteps_100bins";
    //string filenamePrefix = "00a_fcc20x20x20_Jyz1_Jxy1_DMDxyz1_beta1_eqsteps10000_mcsteps_inbin_1000_no_of_bins100";
    //string filenamePrefix = "0aa0_openchain2_beta1_Jsm1_eqsteps10000_mcsteps_inbin_100000_no_of_bins100";
    string filenamePrefix = "0aa0_chain2x2x2_beta1_Jsm1_eqsteps10000_mcsteps_inbin_10000000_no_of_bins100";
    //string filenamePrefix = "bigtest_periodicchain_2p_beta1em6and7000_10000eqsteps_1e6mcsteps_1000bins";
    //string filenamePrefix = "000aa_test_somethingwrongwithcorrfunc_mcstepsinbin100000_bins1000";
    //test_betagenerator(10, 0, 4);
    // Input parameters specifically for run_for_several_betas
    int beta_n = 2;
    double betamin = 1e-6;
    double betamax = 100;
    int betanset = 5;
    vector<double> betas = vector<double>(betanset);
    betas[0] = 0.5; betas[1] = 1.0; betas[2] = 2.0; betas[3] = 10.0; betas[4] = 50.0;

    // By default, the run_for_several_betas-functions do not calculate the correlation function
    //run_for_several_betas(L, eqsteps, mcsteps_inbin, no_of_bins, beta_n, betamin, betamax, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    //run_for_betasgiven(L, eqsteps, mcsteps_inbin, no_of_bins, betanset, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, type_lattice, filenamePrefix, betas, sitestrengthsin, heisenbergin, dm_in);
    //run_for_betasgiven2(L, eqsteps, mcsteps_inbin, no_of_bins, betanset, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, type_lattice, filenamePrefix, betas, sitestrengthsin, heisenbergin, dm_in);

    one_run(L, eqsteps, mcsteps_inbin, no_of_bins, beta, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);

    //run_for_betasgiven_diffdirs(L1, L2, L3, eqsteps, mcsteps_inbin, no_of_bins, betanset, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, type_lattice, filenamePrefix1, betas, sitestrengthsin, heisenbergin, dm_in);

    // Test functions
    //test_fftw(L, sitestrengthsin, heisenbergin, dm_in);
    //test_fcc_extended(L, isotropic, sianisotropy, magfield, dm, periodic, sitestrengthsin, heisenbergin, dm_in);
    //test_fcc_extended_diffdims(L1, L2, L3, isotropic, sianisotropy, magfield, dm, periodic, sitestrengthsin, heisenbergin, dm_in);
    //test_dm();

    // We only need to find the line (0,y,0) once for each fcc L1xL2xL3 Lattice
    //string ylinefilenamePrefix = "fcc2x2x2";
    // Special functions
    //extract_yline(L, ylinefilenamePrefix);

}

void one_run(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    mymc.debugmode(true);
    // Run Metropolis algorithm
    mymc.runmetropolis(beta);
    mymc.endsims();
    //cout << "L = " << L << endl;
}

void run_for_several_betas(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, double betamin, double betamax, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    bool calculatespincorrelationfunction = false;
    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
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

void run_for_betasgiven(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, char type_lattice, string filenamePrefix, vector<double> betas, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    bool calculatespincorrelationfunction = false;
    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    //mymc.debugmode(true);
    //mymc.majordebugtrue();

    for(int i=0; i<beta_n; i++)
    {
        mymc.runmetropolis(betas[i]);
        mymc.reset_energy();
        cout << betas[i] << " done" << endl;
    }
    mymc.endsims();
}

void run_for_betasgiven2(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, char type_lattice, string filenamePrefix, vector<double> betas, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    bool calculatespincorrelationfunction = false;
    // Initializing Monte Carlo
    MonteCarlo mymc(L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    //mymc.debugmode(true);
    //mymc.majordebugtrue();

    for(int i=0; i<beta_n; i++)
    {
        mymc.runmetropolis(betas[i]);
        mymc.reset_energy();
    }
    mymc.endsims();
}

void run_for_betasgiven_diffdirs(int L1, int L2, int L3, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, char type_lattice, string filenamePrefix, vector<double> betas, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    bool calculatespincorrelationfunction = false;
    // Initializing Monte Carlo
    MonteCarlo mymc(L1, L2, L3, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    //mymc.debugmode(true);
    //mymc.majordebugtrue();

    for(int i=0; i<beta_n; i++)
    {
        mymc.runmetropolis(betas[i]);
        mymc.reset_energy();
    }
    mymc.endsims();
}

void extract_yline(int L, string ylinefilenamePrefix)
{
    // Getting the Lattice sites along (0,L,0)
    Lattice mylattice = Lattice(L, false, false, false, false); // We are only interested in the sites
    vector<int> ysites = mylattice.fccyline();

    ofstream ylineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_yline.txt", ylinefilenamePrefix.c_str() );   // Create filename with prefix and ending
    ylineFile.open(filename);
    delete filename;

    int N = ysites.size();

    for(int i=0; i<N; i++)
    {
        ylineFile << ysites[i]  << endl;
    }
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


void test_fftw(int L, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    // All these parameters are irrelevant for this simple test
    int eqsteps = 1000; int mcsteps_inbin = 1000; int no_of_bins= 100;
    bool isotropic = true; bool sianisotropy = false; bool magfield = false; bool dm = false;
    bool periodic = true;

    // These we want to turn off
    bool printeveryMCstep = false;
    bool calculatespincorrelationfunction = true; // To get the q-file

    // Set these
    char type_lattice = 'E';
    string filenamePrefix = "discard"; // Want to know that this file is unimportant

    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    mymc.testFFTW(); // Test FFTW
    mymc.endsims();  // Close the file
}

void test_fcc_extended(int L, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    Lattice mylattice(L, isotropic, sianisotropy, magfield, dm);
    mylattice.setstrengths(sitestrengthsin, heisenbergin, dm_in);
    mylattice.fcc_helical_initialize_extended();

    // Test different sites and their neighbours
    // Neighbours OK
    /*
    int n = 3;
    int no_of_neighbours = mylattice.no_of_neighbours;
    vector<Bond> neighbours =  mylattice.sites[n].bonds;
    vector<double> siteposn = mylattice.sitepositions[n];
    vector<int> sitecoordn = mylattice.sitecoordinates[n];
    cout << "Position , site " << n << " : [" <<  siteposn[0] << "," << siteposn[1] << "," << siteposn[2] << "]" << endl;
    cout << "Coordinates , site " << n << " : [" <<  sitecoordn[0] << "," << sitecoordn[1] << "," << sitecoordn[2] << "]" << endl;

    for(int i=0;i<no_of_neighbours;i++)
    {
        int l = neighbours[i].siteindex2;
        vector<double> siteposl = mylattice.sitepositions[l];
        vector<int> sitecoordl = mylattice.sitecoordinates[l];
        cout << "Position , neighbour " << i << " (site " << l << ") : [" <<  siteposl[0] << "," << siteposl[1] << "," << siteposl[2] << "]" << endl;
        cout << "Coordinates , neighbour " << i << "(site " << l << ") [" <<  sitecoordl[0] << "," << sitecoordl[1] << "," << sitecoordl[2] << "]" << endl;
    }
    */

    // Testing the implementation of J and the sian terms.
    /*
    int no_of_neighbours = mylattice.no_of_neighbours;
    // n = 5
    cout << "For n = 5" << endl;
    vector<Bond> site5bonds = mylattice.sites[5].bonds;
    cout << "Single-ion anisotropy terms: Dix = " << mylattice.sites[5].Dix << "; Diy = " << mylattice.sites[5].Diy << "; Diz = " << mylattice.sites[5].Diz << endl;
    for(int i=0; i<no_of_neighbours; i++)        cout << "i = " << i << "; J = " << site5bonds[i].J  << endl;

    int n;
    if(L<3)    n = 7;
    else       n = 21;
    cout << "Now for " << n << endl;
    vector<Bond> sitebonds = mylattice.sites[n].bonds;
    cout << "Single-ion anisotropy terms: Dix = " << mylattice.sites[n].Dix << "; Diy = " << mylattice.sites[n].Diy << "; Diz = " << mylattice.sites[n].Diz << endl;
    for(int i=0; i<no_of_neighbours; i++)        cout << "i = " << i << "; J = " << sitebonds[i].J  << endl;
    */

    // Testing the position of each point
    int N = mylattice.N;
    for(int i=0; i<N;i++)
    {
        vector<double> pos = mylattice.sitepositions[i];
        vector<int>    crd = mylattice.sitecoordinates[i];
        cout << "Site " << i << endl;
        cout << "Coordinates: [" << crd[0] << " , " << crd[1] << " , " << crd[2] << "]" << endl;
        cout << "Position:    [" << pos[0] << " , " << pos[1] << " , " << pos[2] << "]" << endl << endl;
    }

    // Finding the neighbours of each point:
    int no_of_neighbours = mylattice.no_of_neighbours;
    for(int i=0; i<N;i++)
    {
        cout << "Site " << i << endl;
        for(int j=0; j<no_of_neighbours; j++)
        {
            cout << "Neighbour " << j << ": " << mylattice.sites[i].bonds[j].siteindex2 << endl;
        }
    }

    for(int i=0; i<N; i++)
    {
        for(int j=0; j<no_of_neighbours; j++)
        {
            int n1 = mylattice.sites[i].bonds[j].siteindex1;
            int n2 = mylattice.sites[i].bonds[j].siteindex2;
            cout << "Spin " << i <<  ", bond no. " << j << endl;
            cout << "From site " << n1 << " to site " << n2 << endl;
            cout << "The Heisenberg interaction is: " << mylattice.sites[i].bonds[j].J << endl << endl;

        }
    }
}

void test_fcc_extended_diffdims(int L1, int L2, int L3, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    Lattice mylattice(L1, L2, L3, isotropic, sianisotropy, magfield, dm);
    mylattice.setstrengths(sitestrengthsin, heisenbergin, dm_in);
    mylattice.fcc_helical_initialize_extended();

    // Test different sites and their neighbours
    // Neighbours OK
    /*
    int n = 3;
    int no_of_neighbours = mylattice.no_of_neighbours;
    vector<Bond> neighbours =  mylattice.sites[n].bonds;
    vector<double> siteposn = mylattice.sitepositions[n];
    vector<int> sitecoordn = mylattice.sitecoordinates[n];
    cout << "Position , site " << n << " : [" <<  siteposn[0] << "," << siteposn[1] << "," << siteposn[2] << "]" << endl;
    cout << "Coordinates , site " << n << " : [" <<  sitecoordn[0] << "," << sitecoordn[1] << "," << sitecoordn[2] << "]" << endl;

    for(int i=0;i<no_of_neighbours;i++)
    {
        int l = neighbours[i].siteindex2;
        vector<double> siteposl = mylattice.sitepositions[l];
        vector<int> sitecoordl = mylattice.sitecoordinates[l];
        cout << "Position , neighbour " << i << " (site " << l << ") : [" <<  siteposl[0] << "," << siteposl[1] << "," << siteposl[2] << "]" << endl;
        cout << "Coordinates , neighbour " << i << "(site " << l << ") [" <<  sitecoordl[0] << "," << sitecoordl[1] << "," << sitecoordl[2] << "]" << endl;
    }
    */

    // Testing the implementation of J and the sian terms.
    /*
    int no_of_neighbours = mylattice.no_of_neighbours;
    // n = 5
    cout << "For n = 5" << endl;
    vector<Bond> site5bonds = mylattice.sites[5].bonds;
    cout << "Single-ion anisotropy terms: Dix = " << mylattice.sites[5].Dix << "; Diy = " << mylattice.sites[5].Diy << "; Diz = " << mylattice.sites[5].Diz << endl;
    for(int i=0; i<no_of_neighbours; i++)        cout << "i = " << i << "; J = " << site5bonds[i].J  << endl;

    int n;
    if(L<3)    n = 7;
    else       n = 21;
    cout << "Now for " << n << endl;
    vector<Bond> sitebonds = mylattice.sites[n].bonds;
    cout << "Single-ion anisotropy terms: Dix = " << mylattice.sites[n].Dix << "; Diy = " << mylattice.sites[n].Diy << "; Diz = " << mylattice.sites[n].Diz << endl;
    for(int i=0; i<no_of_neighbours; i++)        cout << "i = " << i << "; J = " << sitebonds[i].J  << endl;
    */

    // Testing the position of each point
    int N = mylattice.N;
    for(int i=0; i<N;i++)
    {
        vector<double> pos = mylattice.sitepositions[i];
        vector<int>    crd = mylattice.sitecoordinates[i];
        cout << "Site " << i << endl;
        cout << "Coordinates: [" << crd[0] << " , " << crd[1] << " , " << crd[2] << "]" << endl;
        cout << "Position:    [" << pos[0] << " , " << pos[1] << " , " << pos[2] << "]" << endl << endl;
    }

    // Finding the neighbours of each point:
    int no_of_neighbours = mylattice.no_of_neighbours;
    for(int i=0; i<N;i++)
    {
        cout << "Site " << i << endl;
        for(int j=0; j<no_of_neighbours; j++)
        {
            cout << "Neighbour " << j << ": " << mylattice.sites[i].bonds[j].siteindex2 << endl;
        }
    }
}

void test_dm()
{

    double x1, y1, z1, x2, y2, z2;
    double Dx, Dy, Dz, J;

    // Input variables:

    // Spin (unnormalized)
    x1 = 1;
    y1 = 2;
    z1 = 3;

    x2 = 3;
    y2 = 5;
    z2 = 1;

    // Direction of DM
    Dx = 1.0;
    Dy = 1.0;
    Dz = 1.0;

    // Heisenberg coupling
    J = 1.0;

    double den1 = sqrt(x1*x1+y1*y1+z1*z1);
    double den2 = sqrt(x2*x2+y2*y2+z2*z2);

    double spin1x = x1/den1;
    double spin1y = y1/den1;
    double spin1z = z1/den1;

    double spin2x = x2/den2;
    double spin2y = y2/den2;
    double spin2z = z2/den2;

    double onetotwo = dm_oneneighbour_fortest(spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, Dx, Dy, Dz);
    double twotoone = dm_oneneighbour_fortest(spin2x, spin2y, spin2z, spin1x, spin1y, spin1z, Dx, Dy, Dz);
    double orthogonalitytest = dm_oneneighbour_fortest(spin1x, spin1y, spin1z, spin1x, spin1y, spin1z, Dx, Dy, Dz);

    cout << "D*(S1 x S2): " << std::setprecision(std::numeric_limits<double>::digits10 + 1) <<  onetotwo << endl;
    cout << "D*(S2 x S1): " << twotoone << endl;
    cout << "Orthogonality test: D*(S1 x S1): " << orthogonalitytest << endl;
    cout << "D*(S1 x S2)+D*(S2 x S1): " << onetotwo+twotoone << endl;

    double Jterm = J_oneneighbour_fortest(spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, J);

    cout << "J*S1*S2 + D*(S1 x S2): " << Jterm+onetotwo << endl;

    // Energy difference, spin flipped
    //Renaming spins to fit the formula
    double sx = spin1x;
    double sy = spin1y;
    double sz = spin1z;

    double sxk = spin2x;
    double syk = spin2y;
    double szk = spin2z;

    // Changing a spin arbitrarily
    double xt, yt, zt, dent;
    xt = 1;    yt = 3;    zt = 1;
    // Normalizing
    dent = sqrt(xt*xt+yt*yt+zt*zt);
    xt = xt/dent; yt = yt/dent; zt = zt/dent;
    cout << "xt: " << xt << "; yt = " << yt << "; zt = " << zt << endl;
    double sx_t = xt; double sy_t= yt; double sz_t = zt;

    double energy_new = dm_oneneighbour_fortest(xt, yt, zt, spin2x, spin2y, spin2z, Dx, Dy, Dz);

    cout << "Energy_new = " << energy_new << endl;
    double energy_diff = -(Dx*((sy-sy_t)*szk-syk*(sz-sz_t))+Dy*((sz-sz_t)*sxk-szk*(sx-sx_t))+Dz*((sx-sx_t)*syk-(sy-sy_t)*sxk));
    cout << "energy difference, dm, according to class MonteCarlo: " << energy_diff << endl;
    cout << "energy difference dm, according to our energies before and after: " << energy_new-onetotwo << endl;
    cout << "Difference between the two approaches: " << energy_diff-(energy_new-onetotwo) << endl;

    cout << "energy_new, extrapolated from MonteCarlo: " << energy_diff+onetotwo << endl;

}

double dm_oneneighbour_fortest(double x1, double y1, double z1, double x2, double y2, double z2, double Dx, double Dy, double Dz)
{
    return Dx*(y1*z2-y2*z1) + Dy*(z1*x2-z2*x1) +Dz*(x1*y2-y1*x2);
}

double J_oneneighbour_fortest(double x1, double y1, double z1, double x2, double y2, double z2, double J)
{
    return J*(x1*x2+y1*y2+z1*z2);
}
