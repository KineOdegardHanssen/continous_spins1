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

// Lattice functions
void extract_yline(int L, string ylinefilenamePrefix); // Extracting the lattice sites along (0,y,0) in the fcc
void extract_xline(int L, string ylinefilenamePrefix);
void extract_zline(int L, string ylinefilenamePrefix);
void lattice_coordinates_straightforward(int L, char type_lattice, string latticefilenamePrefix);
void lattice_coordinates_xyz_lines(int L, string latticefilenamePrefix);
void reciprocallattice_coordinates(int L, char type_lattice, string latticefilenamePrefix);
void reciprocallattice_coordinates_xyzline(int L, char type_lattice, string latticefilenamePrefix);
//void lattice_coordinates_generatingys(int L, string latticefilenamePrefix);



// Test functions
void test_betagenerator(int beta_n, int betamin, int betamax);
void test_fftw(int L, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void test_fftw_againstsims(int L, int eqsteps, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, char type_lattice, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void test_fcc_extended(int L, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void test_fcc_extended_diffdims(int L1, int L2, int L3, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void test_dm();
double dm_oneneighbour_fortest(double x1, double y1, double z1, double x2, double y2, double z2, double Dx, double Dy, double Dz);
double J_oneneighbour_fortest(double x1, double y1, double z1, double x2, double y2, double z2, double J);

// Find out if points are on a plane
vector<double> cross(vector<double> vec1, vector<double> vec2);
vector<double> vecdiffforcross(vector<double> vec1, vector<double> vec2);
void isitonplane(string latticefilenamePrefix);
void isitonline(string latticefilenamePrefix);
void findlinethroughmax(int L, int maxindex, string latticefilenamePrefix);


int main()
{
    bool DEBUG = true;
    if(DEBUG)    cout << "In main" << endl;

    // Input parameters
    int L = 6; // The program is going to be slow if we run for many particles on a 3D lattice

    int L1 = 2;
    int L2 = 2;
    int L3 = 2;

    // bools to determine system type
    bool isotropic    = true;
    bool sianisotropy = false;  // This one does not change its energy unless Dix, Diy and Diz are not all equal.
    bool magfield     = false;
    bool dm           = false;

    // Bool to determine periodicity
    bool periodic     = true; // To determine whether we have periodic boundary conditions or not

    // Selecting the lattice type
    // F: face-centered cubic (fcc); E: fcc with different directions C: cubic; Q:quadratic; O: chain;
    char type_lattice = 'E';
    // If periodic is false, that means we get a grid with open boundary conditions. Currently,
    // that is only implemented for the chain.

    // Setting the strengths of the Hamiltonian terms
    // Magnetic field terms
    double hx = 1;    double hy = 1;    double hz = 1;
    // Single-ion anisotropy terms
    double Dix = 1;    double Diy = 1;    double Diz = 0;
    // Heisenberg term
    double J = 1;
    // Heisenberg terms with varying strengths (for fcc_initialize_extended E)
    double Jy  = 0;    double Jz  = 0;
    double Jxy = -1;    double Jxz = 0;    double Jyz = 1;
    // DM terms
    double Dx = 0;     double Dy = 0;    double Dz = 1;

    vector<double> sitestrengthsin = vector<double>(6);
    sitestrengthsin[0] = hx;    sitestrengthsin[1] = hy;    sitestrengthsin[2] = hz;
    sitestrengthsin[3] = Dix;
    sitestrengthsin[4] = Diy;   sitestrengthsin[5] = Diz;

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
    double beta = 5.0;

    // Run parameters
    int eqsteps = 10000; // Number of steps in the equilibration procedure
    int mcsteps_inbin = 10000; //100000; // MCsteps per bin.
    int no_of_bins = 100;     // The number of bins.

    // Filenames (choose one to use or change slightly)
    //string filenamePrefix = "000aa_test";
    //string filenamePrefix = "0spincorrtest_fcc10x10x10_Dix2Diy2_DxDyDz1_Jyz1_beta0p5_eqsteps50000_mcsteps_inbin_1000_no_of_bins1000";
    //string filenamePrefix = "0fcc_2x2x2_thecomparebetas_eqsteps1000_mcstepsinbin1000_bins100_Jyz1_DixDiy2Diz0";
    //string filenamePrefix = "2pchainopen_J1Dz100_mcsteps_inbin100000_no_of_bins_100";
    //string filenamePrefix = "bigtest_periodicchain_2p_beta0p00001and4000_10000eqsteps_10000mcsteps_1000bins";
    //string filenamePrefix = "bigtest_openchain_5p_beta1em5and4000_10000eqsteps_10000mcsteps_100bins";
    //string filenamePrefix = "00a_fcc20x20x20_Jyz1_Jxy1_DMDxyz1_beta1_eqsteps10000_mcsteps_inbin_1000_no_of_bins100";
    //string filenamePrefix = "0aa0_fcc6x6x6_Jxy1_Jyz1_sianDx1_sianDy1_beta1p025_eqsteps10000_mcsteps_inbin_1000_no_of_bins100";
    string filenamePrefix = "0aa0_fcc6x6x6_Jxym1_Jxz0_Jyz1_beta5_eqsteps10000_mcsteps_inbin_10000_no_of_bins100";
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
    //--------------------------------Running for several betas--------------------------------------//
    //run_for_several_betas(L, eqsteps, mcsteps_inbin, no_of_bins, beta_n, betamin, betamax, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    //run_for_betasgiven(L, eqsteps, mcsteps_inbin, no_of_bins, betanset, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, type_lattice, filenamePrefix, betas, sitestrengthsin, heisenbergin, dm_in);
    //run_for_betasgiven2(L, eqsteps, mcsteps_inbin, no_of_bins, betanset, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, type_lattice, filenamePrefix, betas, sitestrengthsin, heisenbergin, dm_in);
    //run_for_betasgiven_diffdirs(L1, L2, L3, eqsteps, mcsteps_inbin, no_of_bins, betanset, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, type_lattice, filenamePrefix1, betas, sitestrengthsin, heisenbergin, dm_in);

    //----------------------------------Running for one beta-----------------------------------------//
    one_run(L, eqsteps, mcsteps_inbin, no_of_bins, beta, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);


    //-------------------------------------Test functions--------------------------------------------//
    //L = 2;
    //test_fftw(L, sitestrengthsin, heisenbergin, dm_in);
    //test_fftw_againstsims(L, eqsteps, beta, isotropic, sianisotropy, magfield, dm, periodic, type_lattice, sitestrengthsin, heisenbergin, dm_in);
    //test_fcc_extended(L, isotropic, sianisotropy, magfield, dm, periodic, sitestrengthsin, heisenbergin, dm_in);
    //test_fcc_extended_diffdims(L1, L2, L3, isotropic, sianisotropy, magfield, dm, periodic, sitestrengthsin, heisenbergin, dm_in);
    //test_dm();

    // We only need to find the line (0,y,0) once for each fcc L1xL2xL3 Lattice
    //L = 6;
    string latticefilenamePrefix = "fcc6x6x6";
    //string latticefilenamePrefix = "test";
    // Special functions
    //extract_yline(L, latticefilenamePrefix);
    //extract_xline(L, latticefilenamePrefix);
    //extract_zline(L, latticefilenamePrefix);
    //lattice_coordinates_straightforward(L, type_lattice, latticefilenamePrefix);
    //reciprocallattice_coordinates(L, type_lattice, latticefilenamePrefix);
    //lattice_coordinates_xyz_lines(L, latticefilenamePrefix);
    //reciprocallattice_coordinates_xyzline(L, type_lattice, latticefilenamePrefix);
    //isitonplane(filenamePrefix);
    //isitonline(filenamePrefix);
    string indexfilenamePrefix = "fcc6x6x6_index27";
    int maxindex = 27;
    //findlinethroughmax(L, maxindex, indexfilenamePrefix);

}

void one_run(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    mymc.debugmode(true);
    // Run Metropolis algorithm
    mymc.runmetropolis(beta);
    mymc.endsims();
    //mymc.test_couplings_strengths(); // Just to test it. Can remove later
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

void extract_xline(int L, string ylinefilenamePrefix)
{
    // Getting the Lattice sites along (0,L,0)
    Lattice mylattice = Lattice(L, false, false, false, false); // We are only interested in the sites
    vector<int> xsites = mylattice.fccxline();

    ofstream xlineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_xline.txt", ylinefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xlineFile.open(filename);
    delete filename;

    int N = xsites.size();

    for(int i=0; i<N; i++)
    {
        xlineFile << xsites[i]  << endl;
    }
}

void extract_zline(int L, string ylinefilenamePrefix)
{
    // Getting the Lattice sites along (0,L,0)
    Lattice mylattice = Lattice(L, false, false, false, false); // We are only interested in the sites
    vector<int> zsites = mylattice.fcczline();

    ofstream zlineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_zline.txt", ylinefilenamePrefix.c_str() );   // Create filename with prefix and ending
    zlineFile.open(filename);
    delete filename;

    int N = zsites.size();

    for(int i=0; i<N; i++)
    {
        zlineFile << zsites[i]  << endl;
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

void lattice_coordinates_straightforward(int L, char type_lattice, string latticefilenamePrefix)
{
    bool isotropic = false; bool sianisotropy = false; bool magfield = false; bool dm = false;

    Lattice mylattice = Lattice(L,isotropic, sianisotropy, magfield, dm);

    // Type of Lattice
    if(type_lattice=='F')           mylattice.fcc_helical_initialize();          // F for fcc
    else if(type_lattice=='E')      mylattice.fcc_helical_initialize_extended(); // E for extended
    else if(type_lattice=='C')      mylattice.cubic_helical_initialize();        // C for cubic
    else if(type_lattice=='Q')      mylattice.quadratic_helical_initialize();    // Q for quadratic
    else if(type_lattice=='O')      mylattice.chain_periodic_initialize();       // O for one-dimensional

    // Printing information about the q-vectors straight away
    ofstream xyzFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_xyz.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyzFile.open(filename);
    delete filename;

    xyzFile << "xyz-positions are listed units of grid lengths" << endl;

    cout << "File for printing xyz-coord to file is initiated" << endl;

    double x;
    double y;
    double z;

    int N  = mylattice.N;

    if(mylattice.dim==1) // Not really neccessary. Not double tested, either
    {
        for(int i=0; i<N; i++)    xyzFile << "Running index: " << i << "; Index : i: " << i << "; Position : x = " << i << endl;
    }
    if(mylattice.dim==2)
    {

        cout << "For quadratic, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns = mylattice.sitecoordinates[i];
            vector<double> pos = mylattice.sitepositions[i];
            x = ns[0];
            y = ns[1];
            xyzFile << "Running index: " << i << "; Index : i = " << x << ";   j = " << y << ";   Position : x = ";
            x = pos[0];
            y = pos[1];
            xyzFile << x << ";   y = " << y << endl;
        }
        cout << "Done printing to xyzFile" << endl;

    }

    else if(mylattice.dim==3)
    {
        cout << "For cubic or fcc, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns     = mylattice.sitecoordinates[i];
            vector<double> pos = mylattice.sitepositions[i];
            //cout << "ns retrieved" << endl;
            // These must be changed if we change into possibly setting L1, L2, L3 different
            // Don't really need b1, b2, b3, could just use a1, a2, a3 multiplied by 2*M_PI...
            // Could be more general, but we don't need to play around with our lattices that much...
            x = ns[0];
            y = ns[1];
            z = ns[2];
            // Print to file. Site number, qx, qy, qz.
            xyzFile << "Running index: " << i << ";  Index :   i = " << x << ";   j = " << y <<  ";   k = " << z <<"; Position :   x = ";
            x = pos[0];
            y = pos[1];
            z = pos[2];
            xyzFile << x << ";   y = " << y << ";   z = " << z << endl;
        }
        //cout << "Done printing to qFile" << endl;
    }
    xyzFile.close();
}

void lattice_coordinates_xyz_lines(int L, string latticefilenamePrefix)
{
    bool isotropic = false; bool sianisotropy = false; bool magfield = false; bool dm = false;
    Lattice mylattice = Lattice(L,isotropic, sianisotropy, magfield, dm);
    mylattice.fcc_helical_initialize_extended();

    ofstream xyzFilex;
    char *filenamex = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamex, "%s_xyzxline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyzFilex.open(filenamex);
    delete filenamex;

    ofstream xyzFiley;
    char *filenamey = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamey, "%s_xyzyline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyzFiley.open(filenamey);
    delete filenamey;

    ofstream xyzFilez;
    char *filenamez = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamez, "%s_xyzzline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyzFilez.open(filenamez);
    delete filenamez;

    xyzFilex << "xyz-positions are listed units of grid lengths" << endl;
    xyzFiley << "xyz-positions are listed units of grid lengths" << endl;
    xyzFilez << "xyz-positions are listed units of grid lengths" << endl;

    cout << "File for printing xyz-coord a la x-, y- and z-line to file is initiated" << endl;

    int N  = mylattice.N;

    int i,j,k;
    for(int n=0; n<N; n++)
    {
        vector<int> ns = mylattice.sitecoordinates[n];
        i = ns[0];
        j = ns[1];
        k = ns[2];

        vector<double> posx = mylattice.giveposition_fcc_lines(i,j,k,'x');
        vector<double> posy = mylattice.giveposition_fcc_lines(i,j,k,'y');
        vector<double> posz = mylattice.giveposition_fcc_lines(i,j,k,'z');

        // Print to file. Site number, qx, qy, qz.
        double xlx = posx[0];    double ylx = posy[0];    double zlx = posz[0];
        double xly = posx[1];    double yly = posy[1];    double zly = posz[1];
        double xlz = posx[2];    double ylz = posz[2];    double zlz = posz[2];

        xyzFilex << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        xyzFilex << "; Position: x = " << xlx << "; y = " << xly << "; z = " << xlz << endl;
        xyzFiley << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        xyzFiley << "; Position: x = " << ylx << "; y = " << yly << "; z = " << ylz << endl;
        xyzFilez << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        xyzFilez << "; Position: x = " << zlx << "; y = " << zly << "; z = " << zlz << endl;
    }
    xyzFilex.close();
    xyzFiley.close();
    xyzFilez.close();
}

void reciprocallattice_coordinates(int L, char type_lattice, string latticefilenamePrefix)
{
    bool isotropic = false; bool sianisotropy = false; bool magfield = false; bool dm = false;

    Lattice mylattice = Lattice(L,isotropic, sianisotropy, magfield, dm);

    if(type_lattice=='F')           mylattice.fcc_helical_initialize();          // F for fcc
    else if(type_lattice=='E')      mylattice.fcc_helical_initialize_extended(); // E for extended
    else if(type_lattice=='C')      mylattice.cubic_helical_initialize();        // C for cubic
    else if(type_lattice=='Q')      mylattice.quadratic_helical_initialize();    // Q for quadratic
    else if(type_lattice=='O')      mylattice.chain_periodic_initialize();       // O for one-dimensional

    // Printing information about the q-vectors straight away
    ofstream qFile;
    char *filenameq = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenameq, "%s_verbose_qs.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    qFile.open(filenameq);
    delete filenameq;

    cout << "File for printing q-vector to file is initiated" << endl;

    double qx;
    double qy;
    double qz;

    int L1 = mylattice.L1;
    int L2 = mylattice.L2;
    int L3 = mylattice.L3;
    int N  = mylattice.N;

    if(mylattice.dim==1) // Not really neccessary. Not double tested, either
    {
        for(int i=0; i<N; i++)    qFile << "Running index = " << i << ";  qx = " << 2*M_PI/L1*i << endl;
    }
    if(mylattice.dim==2)
    {

        cout << "For quadratic, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns = mylattice.sitecoordinates[i];
            cout << "ns retrieved" << endl;
            qx = ns[0]*mylattice.b1[0]/L1; // Not so sure about this, just need it to compile
            qy = ns[1]*mylattice.b2[1]/L2; // Double check
            cout << "qvec set" << endl;
            qFile << "Running index = " << i << "; Indices: i = " << ns[0] << "; j =  " << ns[1] << "; Positions in q-space:  qx = " << qx << ";  qy =  " << qy << " " << endl;

        }
        cout << "Done printing to qFile" << endl;

    }

    else if(mylattice.dim==3)
    {
        cout << "For cubic or fcc, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns = mylattice.sitecoordinates[i];
            //cout << "ns retrieved" << endl;
            // These must be changed if we change into possibly setting L1, L2, L3 different
            // Don't really need b1, b2, b3, could just use a1, a2, a3 multiplied by 2*M_PI...
            // Could be more general, but we don't need to play around with our lattices that much...
            qx = ns[0]*mylattice.b1[0]/L1 + ns[2]*mylattice.b3[0]/L3;
            qy = ns[0]*mylattice.b1[1]/L1 + ns[1]*mylattice.b2[1]/L2;
            qz = ns[1]*mylattice.b2[2]/L2 + ns[2]*mylattice.b3[2]/L3;
            // Print to file. Site number, qx, qy, qz.
            qFile << "Running index = " << i << "; Indices: i = " << ns[0] << ";  j = " << ns[1] << "; k = " << ns[2] << "; Positions in q-space:   qx = " << qx << ";   qy  = " << qy << ";    qz = " << qz << endl;
        }
        //cout << "Done printing to qFile" << endl;
    }
    qFile.close();
}

void reciprocallattice_coordinates_xyzline(int L, char type_lattice, string latticefilenamePrefix)
{
    bool isotropic = false; bool sianisotropy = false; bool magfield = false; bool dm = false;
    Lattice mylattice = Lattice(L,isotropic, sianisotropy, magfield, dm);
    mylattice.fcc_helical_initialize_extended(); // We default this as well.

    ofstream qFilex;
    char *filenamex = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamex, "%s_xlineq_verbose.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    qFilex.open(filenamex);
    delete filenamex;

    ofstream qFiley;
    char *filenamey = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamey, "%s_ylineq_verbose.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    qFiley.open(filenamey);
    delete filenamey;

    ofstream qFilez;
    char *filenamez = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamez, "%s_zlineq_verbose.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    qFilez.open(filenamez);
    delete filenamez;

    qFilex << "q-vectors are listed units of grid lengths" << endl;
    qFiley << "q-vectors are listed units of grid lengths" << endl;
    qFilez << "q-vectors are listed units of grid lengths" << endl;

    cout << "Files for printing q-vector a la x-, y- and z-line to file are initiated" << endl;

    int N  = mylattice.N;

    int i, j, k;
    for(int n=0; n<N; n++)
    {
        vector<int> ns = mylattice.sitecoordinates[n];
        i = ns[0];
        j = ns[1];
        k = ns[2];

        vector<double> qvecx = mylattice.giveqvector_fcc_lines(i,j,k,'x');
        vector<double> qvecy = mylattice.giveqvector_fcc_lines(i,j,k,'y');
        vector<double> qvecz = mylattice.giveqvector_fcc_lines(i,j,k,'z');

        // Print to file. Site number, qx, qy, qz.
        double qxlx = qvecx[0];    double qylx = qvecy[0];    double qzlx = qvecz[0];
        double qxly = qvecx[1];    double qyly = qvecy[1];    double qzly = qvecz[1];
        double qxlz = qvecx[2];    double qylz = qvecy[2];    double qzlz = qvecz[2];

        qFilex << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        qFilex << "; Vector: qx = " << qxlx << "; qy = " << qxly << "; qz = " << qxlz << endl;
        qFiley << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        qFiley << "; Vector: qx = " << qylx << "; qy = " << qyly << "; qz = " << qylz << endl;
        qFilez << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        qFilez << "; Vector: qx = " << qzlx << "; qy = " << qzly << "; qz = " << qzlz << endl;
    }
    qFilex.close();
    qFiley.close();
    qFilez.close();

}


void test_fftw(int L, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    // All these parameters are irrelevant for this simple test
    int eqsteps = 1000; int mcsteps_inbin = 1000; int no_of_bins= 100;
    bool isotropic = true; bool sianisotropy = false; bool magfield = false; bool dm = false;
    bool periodic = true;

    // These we want to turn off
    bool printeveryMCstep = false;
    bool calculatespincorrelationfunction = false;

    // Set these
    char type_lattice = 'C';
    string filenamePrefix = "discard"; // Want to know that this file is unimportant

    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    mymc.testFFTW(); // Test FFTW
    mymc.endsims();  // Close the file
}

void test_fftw_againstsims(int L, int eqsteps, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, char type_lattice, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    // Input irrelevant for this test
    int mcsteps_inbin = 20000;
    int no_of_bins    = 1;

    // These we want to turn off
    bool printeveryMCstep = false;
    bool calculatespincorrelationfunction = false; // We don't need this for the testing function

    // So it's easy to throw it away afterwards
    string filenamePrefix = "discard"; // Want to know that this file is unimportant

    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    mymc.compareFFTW_withmanual(beta); // Test FFTW
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

// Plane test function
vector<double> cross(vector<double> vec1, vector<double> vec2)
{
    vector<double> returnvec(4);
    returnvec[0] = 0; // To give it the same dimension as the others
    returnvec[1] = vec1[2]*vec2[3]-vec2[2]*vec1[3];
    returnvec[2] = vec1[3]*vec2[1]-vec2[3]*vec1[1];
    returnvec[3] = vec1[1]*vec2[2]-vec2[1]*vec1[2];
    return returnvec;
}

vector<double> vecdiffforcross(vector<double> vec1, vector<double> vec2)
{
    vector<double> returnvec(4);
    returnvec[0] = 0; // To give it the same dimension as the others
    for(int i=1; i<4; i++)
    {
        returnvec[1] = vec1[1]-vec2[1];
        returnvec[2] = vec1[2]-vec2[2];
        returnvec[3] = vec1[3]-vec2[3];
    }
    return returnvec;
}

void isitonplane(string latticefilenamePrefix)
{
    ofstream Fileforplanes;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_investigateplanes.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    Fileforplanes.open(filename);
    delete filename;

    vector<vector<double> > storage;
    // Specify q-vectors here.
    vector<double> v1(4), v2(4), v3(4), v4(4), v5(4), v6(4), v7(4), v8(4), v9(4), v10(4), v11(4), v12(4), v13(4), v14(4), v15(4);
    v1[0] = 21;   v1[1] = -1./6*M_PI;  v1[2] = 1./6*M_PI;  v1[3] =          0;
    v2[0] = 58;   v2[1] = -1./6*M_PI;  v2[2] = 2./3*M_PI;  v2[3] =  1./6*M_PI;
    v3[0] = 63;   v3[1] = -1./3*M_PI;  v3[2] = 5./6*M_PI;  v3[3] =  1./6*M_PI;
    v4[0] = 95;   v4[1] =  1./6*M_PI;  v4[2] = 5./6*M_PI;  v4[3] =  1./3*M_PI;
    v5[0] = 105;  v5[1] = -1./6*M_PI;  v5[2] = 7./6*M_PI;  v5[3] =  1./3*M_PI;
    v6[0] = 111;  v6[1] =          0;  v6[2] = 1./2*M_PI;  v6[3] = -1./2*M_PI;
    v7[0] = 118;  v7[1] =  1./6*M_PI;  v7[2] = 2./3*M_PI;  v7[3] = -1./6*M_PI;
    v8[0] = 125;  v8[1] =  1./3*M_PI;  v8[2] = 5./6*M_PI;  v8[3] =  1./6*M_PI;
    v9[0] = 126;  v9[1] =  1./2*M_PI;  v9[2] =      M_PI;  v9[3] =  1./2*M_PI;
    v10[0] = 133; v10[1] = -1./3*M_PI; v10[2] = 7./6*M_PI; v10[3] = -1./6*M_PI;
    v11[0] = 140; v11[1] = -1./6*M_PI; v11[2] = 4./3*M_PI; v11[3] =  1./6*M_PI;
    v12[0] = 153; v12[1] =  1./6*M_PI; v12[2] = 5./6*M_PI; v12[3] = -1./3*M_PI;
    v13[0] = 163; v13[1] = -1./2*M_PI; v13[2] = 7./6*M_PI; v13[3] = -1./3*M_PI;
    v14[0] = 195; v14[1] =  1./3*M_PI; v14[2] = 7./6*M_PI; v14[3] = -1./6*M_PI;
    v15[0] = 200; v15[1] =  1./6*M_PI; v15[2] = 4./3*M_PI; v15[3] = -1./6*M_PI;

    storage.push_back(v1);
    storage.push_back(v2);
    storage.push_back(v3);
    storage.push_back(v4);
    storage.push_back(v5);
    storage.push_back(v6);
    storage.push_back(v7);
    storage.push_back(v8);
    storage.push_back(v9);
    storage.push_back(v10);
    storage.push_back(v11);
    storage.push_back(v12);
    storage.push_back(v13);
    storage.push_back(v14);
    storage.push_back(v15);

    for(int i=1; i<15; i++)
    {
        for(int j=0; j<i; j++)
        {
            for(int k=1; k<15; k++)
            {
                for(int l=0; l<k; l++)
                {
                    if(k!=i && k!=j && l!=i && l!=j)
                    {
                        vector<double> vdfc1 = vecdiffforcross(storage[i], storage[j]);
                        vector<double> vdfc2 = vecdiffforcross(storage[k], storage[l]);

                        vector<double> thenormal = cross(vdfc1, vdfc2);
                        vector<double> pinp1 = storage[i];
                        vector<double> pinp2 = storage[j];
                        vector<double> pinp3 = storage[k];
                        vector<double> pinp4 = storage[l];
                        double a = thenormal[1]; // zero is for the index number
                        double b = thenormal[2];
                        double c = thenormal[3];
                        double d = a*pinp1[1]+b*pinp1[2]+c*pinp1[3];
                        // Shouls probably print to file.
                        Fileforplanes << "Plane spanned by points " << pinp1[0]<< " , " << pinp2[0] << " , " << pinp3[0]<< " , " << pinp4[0] << ":" << endl;
                        for(int m=0; m<15;m++)
                        {
                            if(m!=i && m!=j && m!=k && m!=l)
                            {
                                vector<double> snake = storage[m];
                                double phew = a*snake[1]+b*snake[2]+c*snake[3]-d;
                                if(abs(phew)<1e-14)    Fileforplanes << "Point " << snake[0] << " on plane. Diff from 0: " << phew << endl;
                            } // End if-test "must be a point we haven't tried yet"
                        } // End check all points
                    } // End "all points in the vectors must be different"
                    Fileforplanes << endl;
                } // End loop over l
            } // End loop over k
        } // End loop over j
    } // End loop over i
    Fileforplanes.close();
}

void isitonline(string latticefilenamePrefix)
{
    ofstream Fileforlines;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_investigatelines.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    Fileforlines.open(filename);
    delete filename;

    vector<vector<double> > storage;
    // Specify q-vectors here.
    // There is a factor of pi everywhere, but we drop it to ease computations
    vector<double> v1(4), v2(4), v3(4), v4(4), v5(4), v6(4), v7(4), v8(4), v9(4), v10(4), v11(4), v12(4), v13(4), v14(4), v15(4);
    v1[0] = 21;   v1[1] = -1./6;  v1[2] = 1./6;  v1[3] =     0;
    v2[0] = 58;   v2[1] = -1./6;  v2[2] = 2./3;  v2[3] =  1./6;
    v3[0] = 63;   v3[1] = -1./3;  v3[2] = 5./6;  v3[3] =  1./6;
    v4[0] = 95;   v4[1] =  1./6;  v4[2] = 5./6;  v4[3] =  1./3;
    v5[0] = 105;  v5[1] = -1./6;  v5[2] = 7./6;  v5[3] =  1./3;
    v6[0] = 111;  v6[1] =     0;  v6[2] = 1./2;  v6[3] = -1./2;
    v7[0] = 118;  v7[1] =  1./6;  v7[2] = 2./3;  v7[3] = -1./6;
    v8[0] = 125;  v8[1] =  1./3;  v8[2] = 5./6;  v8[3] =  1./6;
    v9[0] = 126;  v9[1] =  1./2;  v9[2] =    1;  v9[3] =  1./2;
    v10[0] = 133; v10[1] = -1./3; v10[2] = 7./6; v10[3] = -1./6;
    v11[0] = 140; v11[1] = -1./6; v11[2] = 4./3; v11[3] =  1./6;
    v12[0] = 153; v12[1] =  1./6; v12[2] = 5./6; v12[3] = -1./3;
    v13[0] = 163; v13[1] = -1./2; v13[2] = 7./6; v13[3] = -1./3;
    v14[0] = 195; v14[1] =  1./3; v14[2] = 7./6; v14[3] = -1./6;
    v15[0] = 200; v15[1] =  1./6; v15[2] = 4./3; v15[3] = -1./6;

    storage.push_back(v1);
    storage.push_back(v2);
    storage.push_back(v3);
    storage.push_back(v4);
    storage.push_back(v5);
    storage.push_back(v6);
    storage.push_back(v7);
    storage.push_back(v8);
    storage.push_back(v9);
    storage.push_back(v10);
    storage.push_back(v11);
    storage.push_back(v12);
    storage.push_back(v13);
    storage.push_back(v14);
    storage.push_back(v15);


    vector<double> vec1(4);
    vector<double> vec2(4);
    vector<double> linevec(3);
    vector<double> drawnvec(3);
    vector<double> comparevec(4);

    for(int i=1; i<15;i++)
    {
        vec1 = storage[i];
        for(int j=0; j<i; j++)
        {
            vec2 = storage[j];
            // Making linevec:
            for(int k=1; k<4;k++)    linevec[k-1] = vec1[k]-vec2[k];

            Fileforlines << "Points on the same line as point " << vec2[0] << " and " << vec1[0] << ":";
            for(int k=0; k<15; k++)
            {
                comparevec = storage[k];
                if(k!=i && k!=j)
                {
                    for(int k=1; k<4;k++)    drawnvec[k-1] = comparevec[k]-vec1[k];

                    double yo = drawnvec[0]/linevec[0];
                    double ya = drawnvec[1]/linevec[1];
                    double yi = drawnvec[2]/linevec[2];
                    if(abs(yo-ya)<1e-14)
                    {
                        if(abs(yo-yi)<1e-14)
                        {
                            if(abs(yi-ya)<1e-14)
                            {
                                Fileforlines << " " << comparevec[0] << " ,";
                                //cout << "abs(yo-ya): " << abs(yo-ya) << "; abs(yo-yi): " << abs(yo-yi) << "; abs(yi-ya): " << abs(yi-ya) << endl;
                            } // End iftest yiya
                        } // End if-test yoyi
                    } // End if-test yoya
                } // End if-test to check if k is already drawn
            } // End loop over points k that my be on the line
            Fileforlines << endl;
        } // End loop over j
    } // End loop over i

    Fileforlines.close();
}

void findlinethroughmax(int L, int maxindex, string latticefilenamePrefix)
{   // This function is only for 3D lattices (Made with fcc in mind)
    ofstream linethroughmaxFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_linesthroughmax.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    linethroughmaxFile.open(filename);
    delete filename;

    Lattice mylattice = Lattice(L, false, false, false, false);
    mylattice.fcc_helical_initialize_extended();

    int N = mylattice.N;
    /*
    int L1 = mylattice.L1;
    int L2 = mylattice.L2;
    int L3 = mylattice.L3;
    */

    vector<vector<double> > storage;
    vector<double> qs(3);
    for(int i=0; i<N; i++)
    {   // Finding all points
        // Unshifted lattice:

        vector<int> ns = mylattice.sitecoordinates[i];
        /*
        qs[0] = ns[0]*mylattice.b1[0]/L1 + ns[2]*mylattice.b3[0]/L3;
        qs[1] = ns[0]*mylattice.b1[1]/L1 + ns[1]*mylattice.b2[1]/L2;
        qs[2] = ns[1]*mylattice.b2[2]/L2 + ns[2]*mylattice.b3[2]/L3;
        */

        // Or we could use the shifted lattice:
        qs = mylattice.giveqvector_fcc_lines(ns[0], ns[1], ns[2], 'y');

        storage.push_back(qs);
        // Print to file. Site number, qx, qy, qz.
    }

    vector<double> vec1(3);
    vector<double> vec2(3);
    vector<double> linevec(3);
    vector<double> vec3(3);
    vector<double> comparevec(3);
    vec1 = storage[maxindex];
    for(int i=0; i<N; i++)
    {
        // Making a line of all the points with our max value point.
        if(i!=maxindex)
        {
            vec2 = storage[i];

            // Making linevec:
            for(int k=0; k<3;k++)    linevec[k] = vec1[k]-vec2[k];

            linethroughmaxFile << "Points on the same line as point " << maxindex << " and " << i << ":";
            for(int k=0; k<N; k++)
            {
                vec3 = storage[k];
                if(k!=i && k!=maxindex)
                {
                    for(int k=0; k<3;k++)    comparevec[k] = vec3[k]-vec1[k];

                    double yo = comparevec[0]/linevec[0];
                    double ya = comparevec[1]/linevec[1];
                    double yi = comparevec[2]/linevec[2];
                    if(abs(yo-ya)<1e-14)
                    {
                        if(abs(yo-yi)<1e-14)
                        {
                            if(abs(yi-ya)<1e-14)
                            {
                                linethroughmaxFile << " " << k << " ,";
                                cout << "abs(yo-ya): " << abs(yo-ya) << "; abs(yo-yi): " << abs(yo-yi) << "; abs(yi-ya): " << abs(yi-ya) << endl;
                            } // End iftest yiya
                        } // End if-test yoyi
                    } // End if-test yoya
                } // End if-test to check if k is already drawn
            } // End loop over points k that my be on the line
            linethroughmaxFile << endl;
        }
    }
    linethroughmaxFile.close();

}
