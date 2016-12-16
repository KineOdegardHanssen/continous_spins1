#ifndef LATTICE_H
#define LATTICE_H
#include <bond.h>
#include <site.h>

class Lattice
{
public:
    Lattice(int L);

    int L,N;

    //std::vector<Bond> bonds;
    std::vector<Site> sites;

    //Initialization functions
    fcc_helical_initialize();
};

#endif // LATTICE_H
