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

    //Lattice grid functions
    fcc_helical_initialize();

    // Feed interaction functions
    std::vector<double> givethesiteints(double Dix, double Diy, double Diz, double hx, double hy, double hz, bool sianisotropy, bool magfield);
    std::vector<double> givethebondints(double J, double Dx, double Dy, double Dz, bool isotropic, bool dm);

};

#endif // LATTICE_H
