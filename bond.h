#ifndef BOND_H
#define BOND_H
#include<iostream>
#include <vector>
//#include <site.h> // Or do I?

class Bond
{
public:

    Bond();
    Bond(double J, double Dx, double Dy, double Dz, std::vector site1index, std::vector site2index);

    // Bond quantities
    double J, Dx, Dy, Dz;

    // Bond type
    char bondtype; // Or something
    double bondstrength; // feed this in by Lattice.

    // The sites which the bond connects
    int siteindex1, siteindex2;
    std::vector<int> site1indexvec, site2indexvec; // or something along these lines
};

#endif // BOND_H
