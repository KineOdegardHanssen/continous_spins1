#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <string>
#include <bond.h>
#include <site.h>
#include <lattice.h>


using namespace std;
using std::ofstream; using std::string;

// Energy functions
double sianisotropy_energy(int i, double sx, double sy, double sz, Lattice mylattice);   // Should probably send no of neighbours in. Or attach it to Site, in case I will look at open boundary conditions.
double magfield_energy(int i, double sx, double sy, double sz, Lattice mylattice);
double isotropic_energy(int i, double sx, double sy, double sz, Lattice mylattice);
double dm_energy(int i, double sxi, double syi, double szi, Lattice mylattice);


void MonteCarlo();
double mcstepf(int no_of_neighbours, double N, double beta, double energy_old, bool sianisotropy, bool magfield, bool isotropic, bool dm, Lattice &mylattice, std::default_random_engine generator_u, std::default_random_engine generator_v, std::default_random_engine generator_n, std::default_random_engine generator_prob,  std::uniform_real_distribution<double> distribution_prob, std::uniform_real_distribution<double> distribution_u, std::uniform_real_distribution<double> distribution_v, std::uniform_int_distribution<int> distribution_n);
// Is the next one neccessary?
double energy(bool isotropic, bool sianisotropy, bool magfield, bool dm, Lattice mylattice);

int main()   // main. Monte Carlo steps here?
{
    bool DEBUG = true;
    bool BEDBUG = false; // The bedbug is defeated
    bool HUMBUG = false; // The humbug is defeated
    bool LADYBUG = false; // The ladybug is always cute and dandy
    if(DEBUG)    cout << "In main" << endl;

    // Input parameters
    int L = 10; // The program is going to be slower than before as we have a 3D lattice
    bool isotropic = true;
    bool sianisotropy = true;
    bool magfield = false;
    bool dm = false;
    double beta = 0.1; // Just setting a beta.

    // Run parameters
    int eqsteps = 10000; // Number of steps in the equilibration procedure
    int mcsteps_inbin = 1000; // MCsteps per bin. Do I need bins?
    int no_of_bins = 100;     // The number of bins.

    // Other variables
    double start_clock, end_clock;

    // Opening file to print to
    ofstream printFile;
    //string filenamePrefix = "test10x10x10_cubic";
    string filenamePrefix = "alejandro";
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_cspinMC.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    printFile.open(filename);
    delete filename;

    ofstream bigFile;
    //string filenamePrefix = "test10x10x10_cubic";
    string filenamePrefixb = "alicia";
    char *filenameb = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenameb, "%s_dev_energyav.txt", filenamePrefixb.c_str() );   // Create filename with prefix and ending
    bigFile.open(filenameb);
    delete filenameb;

    if(DEBUG)    cout << "Parameters set" << endl;

    // Setting up the lattice with site parameters and interactions

    start_clock = clock();
    // Initializing instance of class Lattice
    Lattice mylattice = Lattice(L, isotropic, sianisotropy, magfield, dm);
    if(DEBUG)    cout << "Instance of class Lattice initialized" << endl;
    // Choosing type of lattice
    mylattice.cubic_helical_initialize();
    // Should I be storing the number of neighbours in Lattice or Site? Probably a good idea. But for know:
    end_clock = clock();
    double total_time_initialize_lattice = (end_clock - start_clock)/(double) CLOCKS_PER_SEC;
    cout << "Time to initialize Lattice: " << total_time_initialize_lattice  << endl;
    int no_of_neighbours = mylattice.no_of_neighbours;

   if(DEBUG)     cout << "Lattice set up" << endl;

    // Setting the initial energy    // Should we have an own function for calculating the energy?...Probably not
    int N = mylattice.N;
    double energy_old = 0;
    double energy_contribution_sites = 0;
    double energy_contribution_bonds = 0;

    // Put these in functions to be called in the Monte Carlo procedure as well
    // -- Seems to be a lot slower
    start_clock = clock();
    for(int i=0; i<N; i++)
    {
        if(BEDBUG) cout << "In loop to set initial energy" << endl;
        // Contribution from sites
        if(sianisotropy)
        {
            if(BEDBUG)    cout << "In sianisotropy" << endl;
            double Dix = mylattice.sites[i].Dix;
            double Diy = mylattice.sites[i].Diy;
            double Diz = mylattice.sites[i].Diz;
            double sx = mylattice.sites[i].spinx;
            double sy = mylattice.sites[i].spiny;
            double sz = mylattice.sites[i].spinz;
            energy_contribution_sites += Dix*sx*sx + Diy*sy*sy+ Diz*sz*sz;
            //energy_contribution_sites += sianisotropy_energy(i, mylattice);
        }
        if(magfield)
        {
            if(BEDBUG)    cout << "In magfield" << endl;
            double hx = mylattice.sites[i].hx;
            double hy = mylattice.sites[i].hy;
            double hz = mylattice.sites[i].hz;
            double sx = mylattice.sites[i].spinx;
            double sy = mylattice.sites[i].spiny;
            double sz = mylattice.sites[i].spinz;
            energy_contribution_sites -= hx*sx + hy*sy + hz*sz;
        }
        // Contribution from bonds
        if(isotropic)
        {
            if(BEDBUG)   cout << "In isotropic" << endl;
            // Declare no_of_neighbours here in case
            double partnerspinx = 0;
            double partnerspiny = 0;
            double partnerspinz = 0;
            for(int j=0; j<no_of_neighbours; j++)
            {
                if(BEDBUG)    cout << "in loop in isotropic, j = " << j << endl;
                int k = mylattice.sites[i].bonds[j].siteindex2; // Hope I can actually get to this value.
                //                                          // I could, alternatively, just store the index
                //                                          // of bonds. But then I need to organize the bonds
                //                                          // and that is such a hassle.
                if(BEDBUG)    cout << "Have set k in loop, k = " << k << "; N = " << N << endl;
                double J = mylattice.sites[i].bonds[j].J;
                if(BEDBUG)    cout << "Have accessed J in bond between spins" << endl;
                double sx = mylattice.sites[k].spinx;
                double sy = mylattice.sites[k].spiny;
                double sz = mylattice.sites[k].spinz;
                if(BEDBUG)    cout << "Have accessed the components of the spin on the other end" << endl;
                partnerspinx += J*sx;
                partnerspiny += J*sy;
                partnerspinz += J*sz;
                if(BEDBUG)    cout << "Have gathered this contribution into partnerspin" << endl;
            }
            if(BEDBUG)   cout << "Done with the loop in isotropic" << endl;
            double sx = mylattice.sites[i].spinx;
            double sy = mylattice.sites[i].spiny;
            double sz = mylattice.sites[i].spinz;
            energy_contribution_bonds += partnerspinx*sx + partnerspiny*sy + partnerspinz*sz;
            // half this thing. Or find a reasonable way to not double count.
            if(BEDBUG)     cout << "Done with isotropic" << endl;
        }
        if(dm)
        {
            if(BEDBUG)    cout << "In dm" << endl;
            // Double loops and stuff. Could maybe make this more efficient
            double sxi = mylattice.sites[i].spinx;
            double syi = mylattice.sites[i].spiny;
            double szi = mylattice.sites[i].spinz;
            for(int j=0; j<no_of_neighbours; j++)
            {

                int k = mylattice.sites[i].bonds[j].siteindex2; // Hope I can actually get to this value.
                double Dx = mylattice.sites[i].bonds[j].Dx;
                double Dy = mylattice.sites[i].bonds[j].Dy;
                double Dz = mylattice.sites[i].bonds[j].Dz;

                double sxk = mylattice.sites[k].spinx;
                double syk = mylattice.sites[k].spiny;
                double szk = mylattice.sites[k].spinz;

                energy_contribution_bonds -= Dx*(syi*szk-syk*szi)+Dy*(szi*sxk-szk*sxi)+Dz*(sxi*syk-syi*sxk);
            }
        }
        if(BEDBUG) cout << "Done with one, onto the others" << endl;
    }

    if(BEDBUG)   cout << "Done with all" << endl;

    end_clock = clock();
    double total_time_firstenergy = (end_clock - start_clock)/(double) CLOCKS_PER_SEC;
    cout << "Time to initialize energy: " << total_time_firstenergy  << endl;

    energy_old = energy_contribution_sites + energy_contribution_bonds/2.0;

    if(DEBUG)    cout << energy_old << endl;


    // Random number generators
    // For the spin
    // Should I pick random numbers for angles, or components?
    std::default_random_engine generator_u;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_u(0,1);

    std::default_random_engine generator_v;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_v(0,1);

    std::default_random_engine generator_prob;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_prob(0,1);

    // For index. This is given helical boundary conditions, then I only need one index
    std::default_random_engine generator_n;
    std::uniform_int_distribution<int> distribution_n(0,N-1);
    // This should work... Have defined N from the relevant function in instance of class Lattice

    if(DEBUG)    cout << "Done creating random generators" << endl;


    // Equilibration steps
    start_clock = clock();
    for(int i=0; i<eqsteps; i++)
    {
        if(HUMBUG)    cout << "In equilibration steps loop, i = " << i << endl;
        //cout << "In equilibration loop" << endl;
        //double start_clock_i  = clock();
        energy_old = mcstepf(no_of_neighbours, N, beta, energy_old, sianisotropy, magfield, isotropic, dm, mylattice, generator_u, generator_v, generator_n, generator_prob, distribution_prob, distribution_u, distribution_v, distribution_n);
        //double end_clock_i = clock();
        //double equilibration_comptime_i = (end_clock_i - start_clock_i)/(double) CLOCKS_PER_SEC;
        //cout << "i: " << i << "; time to compute mcstep i: " << equilibration_comptime_i << endl;
        if(HUMBUG)    cout << "Done with the equilibration step" << endl;
    }
    end_clock = clock();
    double equilibration_comptime = (end_clock - start_clock)/(double) CLOCKS_PER_SEC;
    cout << "Time to equilibrate: " << equilibration_comptime   << endl;

    start_clock = clock();
    // Monte Carlo steps and measurements
    for(int i=0; i<no_of_bins; i++)
    {   // For each bin
        double energy_av = 0;
        if(LADYBUG)    cout << "i = " << i << "; energy_av before loop: " << energy_av << endl;
        std::vector<double> energies = std::vector<double>(mcsteps_inbin);
        for(int j=0; j<mcsteps_inbin; j++)
        {    // For each mcstep
            energy_old = mcstepf(no_of_neighbours, N, beta, energy_old, sianisotropy, magfield, isotropic, dm, mylattice, generator_u, generator_v, generator_n, generator_prob, distribution_prob, distribution_u, distribution_v, distribution_n);

            // Measurements
            // energy
            energies[j] = energy_old;    // Storing to get the standard deviation
            energy_av +=energy_old;
            bigFile << i << " " << energy_av/(j+1) << endl;

            // Some sort of measurement of the magnetization... How to do this when we have a continuous spin?
        }
        // Energy
        energy_av = energy_av/mcsteps_inbin;

        // Error in the energy //
        double E_stdv = 0;
        for(int i=0; i<no_of_bins; i++)    E_stdv += (energies[i]-energy_av)*(energies[i]-energy_av);
        E_stdv = sqrt(E_stdv/(no_of_bins*(no_of_bins-1)));

        printFile << energy_av << " " << E_stdv << endl;
        // Print to file
    }
    end_clock = clock();
    double mc_comptime = (end_clock - start_clock)/(double) CLOCKS_PER_SEC;
    cout << "Time spent on MC-procedure: " << mc_comptime   << endl;

    if(DEBUG)
    {
        ofstream bondsatsiteFile;
        string filenamePrefix1 = "test10x10x10_cubic";
        char *filename1 = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename1, "%s_bondsatsite.txt", filenamePrefix1.c_str() );   // Create filename with prefix and ending
        bondsatsiteFile.open(filename1);
        delete filename1;

        // Opening file to print to
        ofstream sitesatbondFile;
        string filenamePrefix2 = "test10x10x10_cubic";
        char *filename2 = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename2, "%s_sitesofbond.txt", filenamePrefix2.c_str() );   // Create filename with prefix and ending
        sitesatbondFile.open(filename2);
        delete filename2;

        int N = mylattice.N;

        for(int i=0; i<N; i++)
        {
            bondsatsiteFile << i << " ";
            vector<int> linkstobond = mylattice.bondsofsites[i];
            for(int j=0; j<no_of_neighbours; j++) bondsatsiteFile << linkstobond[j] << " ";
            bondsatsiteFile << endl;
        }

        for(int i=0; i<0.5*no_of_neighbours*N; i++)
        {
            vector<int> linkstosites = mylattice.sitesofbonds[i];
            sitesatbondFile << i << " " << linkstosites[0] << " " << linkstosites[1] << endl;
        }
        bondsatsiteFile.close();
        sitesatbondFile.close();
        printFile.close();
        bigFile.close();
    }
    else
    {
        printFile.close();
        bigFile.close();
    }
}

// Functions which determine the Hamiltonian, feeding in the right interactions.

// Send in all neccessary quantities instead...
double sianisotropy_energy(int i, double sx, double sy, double sz, Lattice mylattice)   // This seemed to be quite slow, unfortunately
{
    double Dix = mylattice.sites[i].Dix;
    double Diy = mylattice.sites[i].Diy;
    double Diz = mylattice.sites[i].Diz;
    double sia_energy_contr = -(Dix*sx*sx + Diy*sy*sy+ Diz*sz*sz);
    return sia_energy_contr;

}

double magfield_energy(int i, double sx, double sy, double sz, Lattice mylattice)
{
    double mf_energy_contr = 0;
    double hx = mylattice.sites[i].hx;
    double hy = mylattice.sites[i].hy;
    double hz = mylattice.sites[i].hz;
    mf_energy_contr -= hx*sx + hy*sy + hz*sz;
    return mf_energy_contr;
}

double isotropic_energy(int i, double sx, double sy, double sz, Lattice mylattice)
{
    // Declare no_of_neighbours here in case
    int no_of_neighbours = 12;
    double iso_energy_contr;
    double partnerspinx = 0;
    double partnerspiny = 0;
    double partnerspinz = 0;
    for(int j=0; j<no_of_neighbours; j++)
    {
        int k = mylattice.sites[i].bonds[j].siteindex2;
        // I could, alternatively, just store the index
        // of bonds. But then I need to organize the bonds
        // and coordinate them with the sites. List of
        //sites that points to the bonds and vice versa.
        double J = mylattice.sites[i].bonds[j].J;
        double sx = mylattice.sites[k].spinx;
        double sy = mylattice.sites[k].spiny;
        double sz = mylattice.sites[k].spinz;
        partnerspinx += J*sx;
        partnerspiny += J*sy;
        partnerspinz += J*sz;
    }
    iso_energy_contr= partnerspinx*sx + partnerspiny*sy + partnerspinz*sz;
    return iso_energy_contr;
}

double dm_energy(int i, double sxi, double syi, double szi, Lattice mylattice)
{
    // Double loops and stuff. Could maybe make this more efficient
    double dm_energy_contr = 0;
    int no_of_neighbours = 12;
    for(int j=0; j<no_of_neighbours; j++)
    {
        int k = mylattice.sites[i].bonds[j].siteindex2; // Hope I can actually get to this value.
        double Dx = mylattice.sites[i].bonds[j].Dx;
        double Dy = mylattice.sites[i].bonds[j].Dy;
        double Dz = mylattice.sites[i].bonds[j].Dz;

        double sxk = mylattice.sites[k].spinx;
        double syk = mylattice.sites[k].spiny;
        double szk = mylattice.sites[k].spinz;

        dm_energy_contr -= Dx*(syi*szk-syk*szi)+Dy*(szi*sxk-szk*sxi)+Dz*(sxi*syk-syi*sxk);

    }
    return dm_energy_contr;
}


// Monte Carlo function
// Or, rather a Simulation class
/*
void MonteCarlo()
{

}
*/

// Should rather call Metropolis
// But have to make sure that mylattice is changed.
// make it double to return the energy?
double mcstepf(int no_of_neighbours, double N, double beta, double energy_old, bool sianisotropy, bool magfield, bool isotropic, bool dm, Lattice &mylattice, std::default_random_engine generator_u, std::default_random_engine generator_v, std::default_random_engine generator_n, std::default_random_engine generator_prob, std::uniform_real_distribution<double> distribution_prob, std::uniform_real_distribution<double> distribution_u, std::uniform_real_distribution<double> distribution_v, std::uniform_int_distribution<int> distribution_n)
{   // Include a counter that measures how many 'flips' are accepted. But what to do with it? Write to file?
    bool HUMBUG = false;  // The humbug is defeated. I think...
    if(HUMBUG)   cout << "In mcstepf. Looping over spins now" << endl;
    for(int n=0; n<N; n++)
    {
        if(HUMBUG)    cout << "Inside loop in mcstepf. n = " << n << endl;
        double energy_new = energy_old;
        double energy_diff = 0;

        int k = distribution_n(generator_n);
        if(HUMBUG)    cout << "Random spin k drawn. k = " << k << endl;

        double sx = mylattice.sites[k].spinx;
        double sy = mylattice.sites[k].spiny;
        double sz = mylattice.sites[k].spinz;

        if(HUMBUG)    cout << "Components of spin " << k << " accessed" << endl;



        //Energy of relevant spin before flip
        // Should I have the random generator here? Or send it in?
        if(sianisotropy)
        {
            if(HUMBUG)    cout << "In sianisotropy in mcstepf" << endl;
            double Dix = mylattice.sites[k].Dix;
            double Diy = mylattice.sites[k].Diy;
            double Diz = mylattice.sites[k].Diz;
            energy_diff += Dix*sx*sx + Diy*sy*sy+ Diz*sz*sz; // This is - originally
            //energy_diff -= sianisotropy_energy(k, sx, sy, sz, mylattice);
        }
        if(magfield)
        {
            if(HUMBUG)    cout << "In magfield in mcstepf" << endl;
            double hx = mylattice.sites[k].hx;
            double hy = mylattice.sites[k].hy;
            double hz = mylattice.sites[k].hz;
            energy_diff += hx*sx + hy*sy + hz*sz;
            //energy_diff -= magfield_energy(k, sx, sy, sz, mylattice);
        }
        if(isotropic)
        {
            if(HUMBUG)    cout << "In isotropic in mcstepf" << endl;
            if(HUMBUG)    cout << "no_of_neighbours = " << no_of_neighbours << endl;
            double partnerspinx = 0;
            double partnerspiny = 0;
            double partnerspinz = 0;
            for(int j=0; j<no_of_neighbours; j++)
            {
                int l = mylattice.sites[k].bonds[j].siteindex2;
                if(HUMBUG)    cout << "Spin no. " << l << " chosen." << endl;
                // I could, alternatively, just store the index
                // of bonds. But then I need to organize the bonds
                // and coordinate them with the sites. List of
                //sites that points to the bonds and vice versa.
                double J = mylattice.sites[k].bonds[j].J;
                double sx = mylattice.sites[l].spinx;
                double sy = mylattice.sites[l].spiny;
                double sz = mylattice.sites[l].spinz;
                partnerspinx += J*sx;
                partnerspiny += J*sy;
                partnerspinz += J*sz;
                if(HUMBUG)    cout << "Why it no work?!" << endl;
            }
            if(HUMBUG)    cout << "Out of that blasted loop!" << endl;
            energy_diff -= partnerspinx*sx + partnerspiny*sy + partnerspinz*sz;
            //energy_diff -= isotropic_energy(k, sx, sy, sz, mylattice);
        }
        if(dm)
        {
            if(HUMBUG)    cout << "In dm in mcstepf" << endl;
            for(int j=0; j<no_of_neighbours; j++)
            {
                int l = mylattice.sites[k].bonds[j].siteindex2; // Hope I can actually get to this value.
                if(HUMBUG)    cout << "Spin no. " << l << " chosen." << endl;


                double Dx = mylattice.sites[k].bonds[j].Dx;
                double Dy = mylattice.sites[k].bonds[j].Dy;
                double Dz = mylattice.sites[k].bonds[j].Dz;
                if(HUMBUG)    cout << "Bonds accessed" << endl;

                double sxk = mylattice.sites[l].spinx;
                double syk = mylattice.sites[l].spiny;
                double szk = mylattice.sites[l].spinz;
                if(HUMBUG)    cout << "Components of spin no. " << l << " accessed." << endl;

                energy_diff += Dx*(sy*szk-syk*sz)+Dy*(sz*sxk-szk*sx)+Dz*(sx*syk-sy*sxk);
            }
            if(HUMBUG)    cout << "Done with the loop in dm in mcstepf" << endl;
            //energy_diff -= dm_energy(k, sx, sy, sz, mylattice);
        }
        if(HUMBUG)    cout << "Done with dm in mcstepf" << endl;

        // Changing the spin (tentatively):
        double u = distribution_u(generator_u);
        double v = distribution_v(generator_v);
        if(HUMBUG)    cout << "Have drawn random numbers in mcstepf" << endl;

        double theta = acos(1.0-2.0*u);
        double phi = 2.0*M_PI*v;

        double sx_t = sin(theta)*cos(phi);
        double sy_t = sin(theta)*sin(phi);
        double sz_t = cos(theta);

        if(HUMBUG)    cout << "Have made a uniform spherical distribution using them" << endl;

        /*
        if(sianisotropy)         energy_diff += sianisotropy_energy(k, sx_t, sy_t, sz_t, mylattice);
        if(magfield)             energy_diff += magfield_energy(k, sx_t, sy_t, sz_t, mylattice);
        if(isotropic)            energy_diff += isotropic_energy(k, sx_t, sy_t, sz_t, mylattice);
        if(dm)                   energy_diff += dm_energy(k, sx_t, sy_t, sz_t, mylattice);
        */

        // Energy contribution after spin change
        if(sianisotropy)
        {
            if(HUMBUG)    cout << "Finding the energy difference from sianisotropy" << endl;
            double Dix = mylattice.sites[k].Dix;
            double Diy = mylattice.sites[k].Diy;
            double Diz = mylattice.sites[k].Diz;
            energy_diff -= Dix*sx_t*sx_t + Diy*sy_t*sy_t+ Diz*sz_t*sz_t; // This is - originally
            //energy_diff -= sianisotropy_energy(k, sx, sy, sz, mylattice);
        }
        if(magfield)
        {
            if(HUMBUG)    cout << "Finging the energy difference from magfield" << endl;
            double hx = mylattice.sites[k].hx;
            double hy = mylattice.sites[k].hy;
            double hz = mylattice.sites[k].hz;
            energy_diff -= hx*sx_t + hy*sy_t + hz*sz_t;
            //energy_diff -= magfield_energy(k, sx, sy, sz, mylattice);
        }
        if(isotropic)
        {
            if(HUMBUG)    cout << "Finging the energy difference from isotropic" << endl;
            double partnerspinx = 0;
            double partnerspiny = 0;
            double partnerspinz = 0;
            for(int j=0; j<no_of_neighbours; j++)
            {
                int l = mylattice.sites[k].bonds[j].siteindex2;
                // I could, alternatively, just store the index
                // of bonds. But then I need to organize the bonds
                // and coordinate them with the sites. List of
                //sites that points to the bonds and vice versa.
                double J = mylattice.sites[k].bonds[j].J;
                double sx = mylattice.sites[l].spinx;
                double sy = mylattice.sites[l].spiny;
                double sz = mylattice.sites[l].spinz;
                partnerspinx += J*sx;
                partnerspiny += J*sy;
                partnerspinz += J*sz;
            }
            energy_diff += partnerspinx*sx_t + partnerspiny*sy_t + partnerspinz*sz_t;
            //energy_diff -= isotropic_energy(k, sx, sy, sz, mylattice);
        }
        if(dm)
        {
            if(HUMBUG)    cout << "Finging the energy difference from dm" << endl;
            for(int j=0; j<no_of_neighbours; j++)
            {
                int l = mylattice.sites[k].bonds[j].siteindex2; // Hope I can actually get to this value.
                double Dx = mylattice.sites[k].bonds[j].Dx;
                double Dy = mylattice.sites[k].bonds[j].Dy;
                double Dz = mylattice.sites[k].bonds[j].Dz;

                double sxk = mylattice.sites[l].spinx;
                double syk = mylattice.sites[l].spiny;
                double szk = mylattice.sites[l].spinz;

                energy_diff -= Dx*(sy_t*szk-syk*sz_t)+Dy*(sz_t*sxk-szk*sx_t)+Dz*(sx_t*syk-sy_t*sxk);
            }
            //energy_diff -= dm_energy(k, sx, sy, sz, mylattice);
        }

        // This should work, but there is probably some error here...
        energy_new = energy_old + energy_diff;
        // Or should I just test if energy_new < 0? ... May have to find energy_new anyways...

        // Updating the energy and the state according to Metropolis
        if(energy_new <= energy_old)
        {
            mylattice.sites[k].spinx = sx_t;
            mylattice.sites[k].spiny = sy_t;
            mylattice.sites[k].spinz = sz_t;
            energy_old = energy_new; // Update energy
        }
        else
        {
            double prob = exp(-beta*(energy_new-energy_old));
            double drawn = distribution_prob(generator_prob);
            if(drawn<prob)
            {
                mylattice.sites[k].spinx = sx_t;
                mylattice.sites[k].spiny = sy_t;
                mylattice.sites[k].spinz = sz_t;
                energy_old = energy_new;  // Update energy
            }
        }
    }
    return energy_old;

}
