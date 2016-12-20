#include <iostream>
#include <bond.h>
#include <site.h>
#include <lattice.h>

using namespace std;

void MonteCarlo();
// Is the next one neccessary?
double energy(bool isotropic, bool sianisotropy, bool magfield, bool dm, Lattice mylattice);

int main()   // main. Monte Carlo steps here?
{
    bool DEBUG = true;
    if(DEBUG)    cout << "In main" << endl;
    // Input parameters
    int L = 10; // The program is going to be slower than before as we have a 3D lattice
    bool isotropic = false;
    bool sianisotropy = false;
    bool magfield = false;
    bool dm = false;

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

                energy_contribution_bonds += Dx*(syi*szk-syk*szi)+Dy*(szi*sxk-szk*sxi)+Dz*(sxi*syk-syi*sxk);
            }
        }
    }

    energy_old = energy_contribution_sites + energy_contribution_bonds;

    if(DEBUG)    cout << energy_old << endl;

    /*
    // These may be more computationally efficient, but uglier. Efficiency is perhaps not as important, as this
    // is just a start-up procedure. This is not the bottleneck
    // Add a lot of stuff.
    if(sianisotropy && magfield)
    {
        for(int i=0; i<N; i++)
        {
            double Dix = mylattice.sites[i].Dix;
            double Diy = mylattice.sites[i].Diy;
            double Diz = mylattice.sites[i].Diz;
            double sx = mylattice.sites[i].spinx;
            double sy = mylattice.sites[i].spiny;
            double sz = mylattice.sites[i].spinz;
            double hx = mylattice.sites[i].hx;
            double hy = mylattice.sites[i].hy;
            double hz = mylattice.sites[i].hz;
            energy_contribution_sites += Dix*sx*sx + Diy*sy*sy+ Diz*sz*sz;
            energy_contribution_sites -= hx*sx + hy*sy + hz*sz;
        }
    }
    else if(sianisotropy && (!magfield))
    {
        for(int i=0; i<N; i++)
        {
            double Dix = mylattice.sites[i].Dix;
            double Diy = mylattice.sites[i].Diy;
            double Diz = mylattice.sites[i].Diz;
            double sx = mylattice.sites[i].spinx;
            double sy = mylattice.sites[i].spiny;
            double sz = mylattice.sites[i].spinz;
            double hx = mylattice.sites[i].hx;
            double hy = mylattice.sites[i].hy;
            double hz = mylattice.sites[i].hz;
            energy_contribution_sites += Dix*sx*sx + Diy*sy*sy+ Diz*sz*sz;
            energy_contribution_sites -= hx*sx + hy*sy + hz*sz;
        }
    }
    */



    // Equilibration steps

    // Monte Carlo steps and measurements

    // Printing to file in some way

}

// Functions which determine the Hamiltonian, feeding in the right interactions.


// Monte Carlo function
