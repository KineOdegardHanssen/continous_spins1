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
    // Input parameters
    int L = 10; // The program is going to be slower than before as we have a 3D lattice
    bool isotropic = true;
    bool sianisotropy = false;
    bool magfield = false;
    bool dm = false;

    // Setting up the lattice with site parameters and interactions

    // Initializing instance of class Lattice
    Lattice mylattice = Lattice(L, isotropic, sianisotropy, magfield, dm);
    // Choosing type of lattice
    mylattice.fcc_helical_initialize();

    // Setting the initial energy    // Should we have an own function for calculating the energy?...Probably not

    // Equilibration steps

    // Monte Carlo steps and measurements

    // Printing to file in some way

}

// Functions which determine the Hamiltonian, feeding in the right interactions.


// Monte Carlo function
