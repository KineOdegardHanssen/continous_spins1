#include <iostream>
#include <random>
#include <bond.h>
#include <site.h>
#include <lattice.h>

using namespace std;

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
    if(DEBUG)    cout << "In main" << endl;

    // Input parameters
    int L = 10; // The program is going to be slower than before as we have a 3D lattice
    bool isotropic = true;
    bool sianisotropy = true;
    bool magfield = false;
    bool dm = false;
    double beta = 0.1; // Just setting a beta.

    // Run parameters
    int eqsteps = 1000; // Number of steps in the equilibration procedure
    int mcsteps_inbin = 1000; // MCsteps per bin. Do I need bins?
    int no_of_bins = 100;     // The number of bins.

    if(DEBUG)    cout << "Parameters set" << endl;

    // Setting up the lattice with site parameters and interactions

    // Initializing instance of class Lattice
    Lattice mylattice = Lattice(L, isotropic, sianisotropy, magfield, dm);
    if(DEBUG)    cout << "Instance of class Lattice initialized" << endl;
    // Choosing type of lattice
    mylattice.fcc_helical_initialize();
    // Should I be storing the number of neighbours in Lattice or Site? Probably a good idea. But for know:
    int no_of_neighbours = 12; // For fcc

   if(DEBUG)     cout << "Lattice set up" << endl;

    // Setting the initial energy    // Should we have an own function for calculating the energy?...Probably not
    int N = mylattice.N;
    double energy_old = 0;
    double energy_contribution_sites = 0;
    double energy_contribution_bonds = 0;

    // Put these in functions to be called in the Monte Carlo procedure as well
    // -- Seems to be a lot slower
    for(int i=0; i<N; i++)
    {
        // Contribution from sites
        if(sianisotropy)
        {
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
            // Declare no_of_neighbours here in case
            double partnerspinx = 0;
            double partnerspiny = 0;
            double partnerspinz = 0;
            for(int j=0; j<no_of_neighbours; j++)
            {
                int k = mylattice.sites[i].bonds[j].siteindex2; // Hope I can actually get to this value.
                //                                          // I could, alternatively, just store the index
                //                                          // of bonds. But then I need to organize the bonds
                //                                          // and that is such a hassle.
                double J = mylattice.sites[i].bonds[j].J;
                double sx = mylattice.sites[k].spinx;
                double sy = mylattice.sites[k].spiny;
                double sz = mylattice.sites[k].spinz;
                partnerspinx += J*sx;
                partnerspiny += J*sy;
                partnerspinz += J*sz;
            }
            double sx = mylattice.sites[i].spinx;
            double sy = mylattice.sites[i].spiny;
            double sz = mylattice.sites[i].spinz;
            energy_contribution_bonds += partnerspinx*sx + partnerspiny*sy + partnerspinz*sz;
            // half this thing. Or find a reasonable way to not double count.
        }
        if(dm)
        {
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
    }

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


    // Equilibration steps
    for(int i=0; i<eqsteps; i++)
    {
        energy_old = mcstepf(no_of_neighbours, N, beta, energy_old, sianisotropy, magfield, isotropic, dm, mylattice, generator_u, generator_v, generator_n, generator_prob, distribution_prob, distribution_u, distribution_v, distribution_n);
    }

    // Monte Carlo steps and measurements
    for(int i=0; i<no_of_bins; i++)
    {   // For each bin
        double energy_av = 0;
        for(int j=0; j<mcsteps_inbin; j++)
        {    // For each mcstep
            energy_old = mcstepf(no_of_neighbours, N, beta, energy_old, sianisotropy, magfield, isotropic, dm, mylattice, generator_u, generator_v, generator_n, generator_prob, distribution_prob, distribution_u, distribution_v, distribution_n);

            // Measurements
            // energy
            energy_av +=energy_old;

            // Some sort of measurement of the magnetization... How to do this when we have a continuous spin?
        }
        energy_av = energy_av/mcsteps_inbin;
        // Print to file
    }
    // Printing to file in some way
}

// Functions which determine the Hamiltonian, feeding in the right interactions.

// Send in all neccessary quantities instead...
double sianisotropy_energy(int i, double sx, double sy, double sz, Lattice mylattice)   // This seemed to be quite slow, unfortunately
{
    double Dix = mylattice.sites[i].Dix;
    double Diy = mylattice.sites[i].Diy;
    double Diz = mylattice.sites[i].Diz;
    double sia_energy_contr = Dix*sx*sx + Diy*sy*sy+ Diz*sz*sz;
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
        return dm_energy_contr;
    }
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
    for(int n=0; n<N; n++)
    {
        double energy_new = energy_old;
        double energy_diff = 0;

        int k = distribution_n(generator_n);

        double sx = mylattice.sites[k].spinx;
        double sy = mylattice.sites[k].spiny;
        double sz = mylattice.sites[k].spinz;

        //Energy of relevant spin before flip
        // Should I have the random generator here? Or send it in?
        if(sianisotropy)         energy_diff -= sianisotropy_energy(k, sx, sy, sz, mylattice);
        if(magfield)             energy_diff -= magfield_energy(k, sx, sy, sz, mylattice);
        if(isotropic)            energy_diff -= isotropic_energy(k, sx, sy, sz, mylattice);
        if(dm)                   energy_diff -= dm_energy(k, sx, sy, sz, mylattice);

        // Changing the spin (tentatively):
        double u = distribution_u(generator_u);
        double v = distribution_v(generator_v);

        double theta = acos(1.0-2.0*u);
        double phi = 2.0*M_PI*v;

        double sx_t = sin(theta)*cos(phi);
        double sy_t = sin(theta)*sin(phi);
        double sz_t = cos(theta);

        if(sianisotropy)         energy_diff += sianisotropy_energy(k, sx_t, sy_t, sz_t, mylattice);
        if(magfield)             energy_diff += magfield_energy(k, sx_t, sy_t, sz_t, mylattice);
        if(isotropic)            energy_diff += isotropic_energy(k, sx_t, sy_t, sz_t, mylattice);
        if(dm)                   energy_diff += dm_energy(k, sx_t, sy_t, sz_t, mylattice);

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
