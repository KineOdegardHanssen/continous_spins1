#ifndef SITE_H
#define SITE_H
#include <iostream>
#include <vector>
#include <bond.h>


class Site
{
public:

    int index; // is this neccessary?
    int lenint; // Length of interaction array
    double hx, hy, hz, Dix, Diy, Diz; // This would lead to a lot of zero elements, possibly
    double spinx, spiny, spinz;

    // Should I have neighbours here?
    //std::vector<int> bondtype1; // Not here
    //std::vector<int> bondtype2;
    //std::vector<int> bondtype3;
    //std::vector<double> spin;   // maybe store in each direction instead
    //std::vector<int> bondindexes; // Or some sort of pointer...
    std::vector<double> siteint; // Array for site interaction
    std::vector<Bond> bonds; // I guess I do have to add a pointer to Bond. But then both classes inherit each other. Weird.

    // Initializers
    Site();
    //Site(int n,  bool sianisotropy, bool magfield, double spinx, double spiny, double spinz, double hx, double hy, double hz, double Dix, double Diy, double Diz, std::vector<Bond> bonds, std::vector<bool> boolvec);
    Site(int n,  bool sianisotropy, bool magfield, double spinx, double spiny, double spinz, std::vector<double> siteint, std::vector<Bond> bonds);
};

#endif // SITE_H
