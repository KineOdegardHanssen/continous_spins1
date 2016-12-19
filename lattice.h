#ifndef LATTICE_H
#define LATTICE_H
#include <bond.h>
#include <site.h>

class Lattice
{
public:
    Lattice(int L);
    bool isotropic,  dm;         // Bools for n.n. terms
    bool sianisotropy, magfield; // Bools for site terms

    int L,N;

    //std::vector<Bond> bonds;
    std::vector<Site> sites;

    //Initialization functions
    fcc_helical_initialize();
};

#endif // LATTICE_H
