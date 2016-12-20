#include <iostream>
#include <random>
#include <bond.h>
#include <site.h>
#include <lattice.h>

using namespace std;

// Energy functions
double sianisotropy_energy(int i, Lattice mylattice);
double magfield_energy(int i, Lattice mylattice);
double isotropic_energy(int i, Lattice mylattice);
double dm_energy(int i, Lattice mylattice);


void MonteCarlo();
double mcstepf(int no_of_neighbours, double N, double beta, double energy_old, bool sianisotropy, bool magfield, bool isotropic, bool dm, Lattice mylattice, std::default_random_engine generator_x, std::default_random_engine generator_y, std::default_random_engine generator_z, std::default_random_engine generator_n,std::default_random_engine generator_prob,  std::uniform_real_distribution<double> distribution_prob, std::uniform_real_distribution<double> distribution_x, std::uniform_real_distribution<double> distribution_y, std::uniform_real_distribution<double> distribution_z, std::uniform_int_distribution<int> distribution_n);
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

    int eqsteps = 1000; // Number of steps in the equilibration procedure

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
    std::default_random_engine generator_x;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_x(0,1);

    std::default_random_engine generator_y;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_y(0,1);

    std::default_random_engine generator_z;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_z(0,1);

    std::default_random_engine generator_prob;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_prob(0,1);

    // For index. This is given helical boundary conditions, then I only need one index
    std::default_random_engine generator_n;
    std::uniform_int_distribution<int> distribution_n(0,N-1);


    // Equilibration steps
    for(int i=0; i<eqsteps; i++)
    {
        energy_old = mcstepf(no_of_neighbours, N, beta, energy_old, sianisotropy, magfield, isotropic, dm, mylattice, generator_x, generator_y, generator_z, generator_n, generator_prob, distribution_prob, distribution_x, distribution_y, distribution_z, distribution_n);
    }

    // Monte Carlo steps and measurements

    // Printing to file in some way

}

// Functions which determine the Hamiltonian, feeding in the right interactions.

double sianisotropy_energy(int i, Lattice mylattice)   // This seemed to be quite slow, unfortunately
{
    double sia_energy_contr;
    double Dix = mylattice.sites[i].Dix;
    double Diy = mylattice.sites[i].Diy;
    double Diz = mylattice.sites[i].Diz;
    double sx = mylattice.sites[i].spinx;
    double sy = mylattice.sites[i].spiny;
    double sz = mylattice.sites[i].spinz;
    sia_energy_contr = Dix*sx*sx + Diy*sy*sy+ Diz*sz*sz;
    return sia_energy_contr;

}

double magfield_energy(int i, Lattice mylattice)
{
    double mf_energy_contr = 0;
    double hx = mylattice.sites[i].hx;
    double hy = mylattice.sites[i].hy;
    double hz = mylattice.sites[i].hz;
    double sx = mylattice.sites[i].spinx;
    double sy = mylattice.sites[i].spiny;
    double sz = mylattice.sites[i].spinz;
    mf_energy_contr -= hx*sx + hy*sy + hz*sz;
    return mf_energy_contr;
}

double isotropic_energy(int i, Lattice mylattice)
{
    // Declare no_of_neighbours here in case
    int no_of_neighbours = 12;
    double iso_energy_contr;
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
    iso_energy_contr= partnerspinx*sx + partnerspiny*sy + partnerspinz*sz;
    return iso_energy_contr;
}

double dm_energy(int i, Lattice mylattice)
{
    // Double loops and stuff. Could maybe make this more efficient
    double dm_energy_contr = 0;
    int no_of_neighbours = 12;
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

        dm_energy_contr -= Dx*(syi*szk-syk*szi)+Dy*(szi*sxk-szk*sxi)+Dz*(sxi*syk-syi*sxk);
        return dm_energy_contr;
    }
}


// Monte Carlo function
/*
void MonteCarlo()
{

}
*/

// Should rather call Metropolis
// But have to make sure that mylattice is changed.
// make it double to return the energy?
double mcstepf(int no_of_neighbours, double N, double beta, double energy_old, bool sianisotropy, bool magfield, bool isotropic, bool dm, Lattice mylattice, std::default_random_engine generator_x, std::default_random_engine generator_y, std::default_random_engine generator_z, std::default_random_engine generator_n, std::default_random_engine generator_prob, std::uniform_real_distribution<double> distribution_prob, std::uniform_real_distribution<double> distribution_x, std::uniform_real_distribution<double> distribution_y, std::uniform_real_distribution<double> distribution_z, std::uniform_int_distribution<int> distribution_n)
{
    for(int n=0; n<N; n++)
    {
        double energy_new = energy_old;
        double energy_diff = 0;
        // Are these really neccessary?
        //double sia_diff = 0;
        //double mf_diff  = 0;
        //double iso_diff = 0;
        //double dm_diff  = 0;

        int k = distribution_n(generator_n);
        //Energy of relevant spin before flip
        // Should I have the random generator here? Or send it in?
        if(sianisotropy)
        {   // I should probably just send these into functions
            double Dix = mylattice.sites[k].Dix;
            double Diy = mylattice.sites[k].Diy;
            double Diz = mylattice.sites[k].Diz;
            double sx = mylattice.sites[k].spinx;
            double sy = mylattice.sites[k].spiny;
            double sz = mylattice.sites[k].spinz;
            energy_diff -= Dix*sx*sx + Diy*sy*sy+ Diz*sz*sz; // This will be changed
            // Or just subtract from energy_old? No, have to test against that...
        }
        if(magfield)
        {
            double hx = mylattice.sites[k].hx;
            double hy = mylattice.sites[k].hy;
            double hz = mylattice.sites[k].hz;
            double sx = mylattice.sites[k].spinx;
            double sy = mylattice.sites[k].spiny;
            double sz = mylattice.sites[k].spinz;
            energy_diff -= hx*sx + hy*sy + hz*sz;
        }
        if(isotropic)
        {   // Should really make a function for these, I guess.
            // Declare no_of_neighbours here in case
            double partnerspinx = 0;
            double partnerspiny = 0;
            double partnerspinz = 0;
            for(int j=0; j<no_of_neighbours; j++)
            {
                int l = mylattice.sites[k].bonds[j].siteindex2; // Hope I can actually get to this value.
                //                                          // I could, alternatively, just store the index
                //                                          // of bonds. But then I need to organize the bonds
                //                                          // and that is such a hassle.
                double J = mylattice.sites[k].bonds[j].J;
                double sx = mylattice.sites[l].spinx;
                double sy = mylattice.sites[l].spiny;
                double sz = mylattice.sites[l].spinz;
                partnerspinx += J*sx;
                partnerspiny += J*sy;
                partnerspinz += J*sz;
            }
            double sx = mylattice.sites[k].spinx;
            double sy = mylattice.sites[k].spiny;
            double sz = mylattice.sites[k].spinz;
            energy_diff -= partnerspinx*sx + partnerspiny*sy + partnerspinz*sz;
        }
        if(dm)
        {   // Or send in to functions
            double sxi = mylattice.sites[k].spinx;
            double syi = mylattice.sites[k].spiny;
            double szi = mylattice.sites[k].spinz;
            for(int j=0; j<no_of_neighbours; j++)
            {
                int i = mylattice.sites[k].bonds[j].siteindex2; // Hope I can actually get to this value.
                double Dx = mylattice.sites[k].bonds[j].Dx;
                double Dy = mylattice.sites[k].bonds[j].Dy;
                double Dz = mylattice.sites[k].bonds[j].Dz;

                double sxk = mylattice.sites[i].spinx;
                double syk = mylattice.sites[i].spiny;
                double szk = mylattice.sites[i].spinz;

                energy_diff -= Dx*(syi*szk-syk*szi)+Dy*(szi*sxk-szk*sxi)+Dz*(sxi*syk-syi*sxk);
            }
        }


        // Updating the energy and the state according to Metropolis
        if(energy_new <= energy_old)
        {
            //state(x,y) = -state(x,y); // This line has to be changed
            energy_old = energy_new; // Update energy
        }
        else
        {
            double prob = exp(-beta*(energy_new-energy_old));
            double drawn = distribution_prob(generator_prob);
            if(drawn<prob)
            {

                //state(x,y) = -state(x,y);  // This line has to be changed
                energy_old = energy_new;  // Update energy
            }
        }
    }

}
