#ifndef BOND_H
#define BOND_H
#include<iostream>
#include <vector>
//#include <site.h> // Or do I?

class Bond
{
public:
    // Bond quantities
    double J, Dx, Dy, Dz;

    // Bond type
    char bondtype; // Or something
    double bondstrength; // feed this in by Lattice.
    std::vector<double> bondints; // use this

    // The sites which the bond connects
    int siteindex1, siteindex2;
    //std::vector<int> site1indexvec, site2indexvec; // or something along these lines

    Bond();
    Bond(int siteindex1, int siteindex2, bool isotropic, bool dm, std::vector<double> bondints);
    //Bond(double J, double Dx, double Dy, double Dz, std::vector<int> site1indexvec, std::vector<int> site2indexvec);
};

#endif // BOND_H
