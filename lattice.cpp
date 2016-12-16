#include "lattice.h"

Lattice::Lattice(int L)
{
    this->L = L;
}

Lattice::fcc_helical_initialize()
{
    N = L*(L+1)*(L+1); // Look this up!
    // Setting up the sites
    // We set up the matrix by having all spins in the same direction. Or maybe draw at random?
    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    // Interactions. Should have a way of choosing which terms we look at. Maybe different initialization
    // functions in site? Send in a char for that instead of having all these ones. Quickly get a lot of
    // unneccessary calculations.
    double hx = 1;
    double hy = 1;
    double hz = 1;
    double Dix = 1;
    double Diy = 1;
    double Diz = 1;

    double J = 1;
    double Dx = 1;
    double Dy = 1;
    double Dz = 1;
    // Could have these inside the loop and add randomness.

    for(int n=0; n<N; n++)
    {
        // Finding the neighbours to n
        // NB!: So far, I have only added neighbours that are a distance 1 apart, 2 being the length of the cell.
        // So no linking two cell corner atoms together as of now.
        // This should only be done once. And that is exactly what we are doing.
        // Doing modulo operations, as suggested in Newman & Barkema
        // Should probably verify this bit
        np1 = (n+1)%N;
        npL = (n+L)%N;
        npL2 = (n+L*L)%N;
        npLm1 = (n+L-1)%N;
        npL2m1 = (n+L*L-1)%N;
        npL2mL = (n+L*L-L)%N;

        // And should DEFINITELY verify this:
        nm1 = (n+N-1)%N;
        nmL = (n+N-L)%N;
        nmL2 = (n+N-L*L)%N;
        nmLm1 = (n+N-L+1)%N;
        nmL2m1 = (n+N-L*L+1)%N;
        nmL2mL = (n+N-L*L+L)%N;

        std::vector<Bond> bonds;

        // Making a lot of bond classes to be added to bonds.
        bonds.push_back(Bond(J, Dx, Dy, Dz, n, np1));  // Do I really need to send in n?
        bonds.push_back(Bond(J, Dx, Dy, Dz, n, nm1));
        bonds.push_back(Bond(J, Dx, Dy, Dz, n, npL));
        bonds.push_back(Bond(J, Dx, Dy, Dz, n, nmL));
        bonds.push_back(Bond(J, Dx, Dy, Dz, n, npL2));
        bonds.push_back(Bond(J, Dx, Dy, Dz, n, nmL2));
        bonds.push_back(Bond(J, Dx, Dy, Dz, n, npLm1));
        bonds.push_back(Bond(J, Dx, Dy, Dz, n, nmLm1));
        bonds.push_back(Bond(J, Dx, Dy, Dz, n, npL2m1));
        bonds.push_back(Bond(J, Dx, Dy, Dz, n, nmL2m1));
        bonds.push_back(Bond(J, Dx, Dy, Dz, n, npL2mL));
        bonds.push_back(Bond(J, Dx, Dy, Dz, n, nmL2mL));

        // Is it too nested to make Site inherit Bond? ... Seems fair?
        sites.push_back(Site(n, spinx, spiny, spinz, hx, hy, hz, Dix, Diy, Diz, bonds));
    }


    // Setting up the bonds
    // Setting up the sites
    // Should I use helical boundary conditions or periodic?
}
