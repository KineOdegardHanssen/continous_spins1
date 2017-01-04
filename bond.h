#ifndef BOND_H
#define BOND_H
#include<iostream>
#include <vector>
//#include <site.h> // Or do I?

class Bond
{
public:

    Bond();
    Bond(bool isotropic, bool dm, std::vector<double> tJs, std::vector<double> tDxes, std::vector<double> tDys, std::vector<double> tDzs);
    //Bond(double J, double Dx, double Dy, double Dz, std::vector<int> site1indexvec, std::vector<int> site2indexvec);


    // Bond quantities
    double J, Dx, Dy, Dz;

    // Bond type. Could use this to simplify
    //char bondtype; // Or something
    //double bondstrength; // feed this in by Lattice.

    //std::vector<double> bondints; // for loading in
    std::vector<double> Js;       // All the Js, indexed like the bonds in Lattice
    std::vector<double> Dxes;
    std::vector<double> Dys;
    std::vector<double> Dzs;

    // The sites which the bond connects
    //This is now taken care of by
    //int siteindex1, siteindex2;
    //std::vector<int> site1indexvec, site2indexvec; // or something along these lines
};

#endif // BOND_H
