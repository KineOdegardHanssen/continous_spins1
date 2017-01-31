#ifndef LATTICE_H
#define LATTICE_H
#include <iostream>
#include <math.h>
#include "bond.h"
#include "site.h"
#include "gaussiandeviate.h"
//#include <armadillo>  // in this case, enter LIBS += -larmadillo -llapack -lblas in .pro-file

using namespace std; // Now, I can remove all std's.

class Lattice
{
public:
    bool isotropic,  dm;         // Bools for n.n. terms
    bool sianisotropy, magfield; // Bools for site terms
    bool notperiodic;            // Bool to keep track of neighbours

    int dim, L, N, no_of_neighbours;
    long int seed;
    // Typedefs
    typedef vector<int> vecint;
    typedef vector<vecint> intmatrix;

    //std::vector<Bond> bonds;
    std::vector<Site> sites;
    //std::vector<std::vector<double> > sitepositions_chain; // In case we introduce a grid length
    std::vector<std::vector<double> > sitepositions;
    std::vector<std::vector<int> >    sitecoordinates;

    // Initialization
    Lattice();
    Lattice(int L, bool isotropic, bool sianisotropy, bool magfield, bool dm);


    //Lattice grid functions
    //void chain_2p_periodic_initialize();
    void chain_periodic_initialize();
    void chain_closed_initialize();
    void quadratic_helical_initialize();
    void cubic_helical_initialize();
    void fcc_helical_initialize();

    // Feed interaction functions
    std::vector<double> givethesiteints(double Dix, double Diy, double Diz, double hx, double hy, double hz, bool sianisotropy, bool magfield);
    std::vector<double> givethebondints(double J, double Dx, double Dy, double Dz, bool isotropic, bool dm);

};

#endif // LATTICE_H
