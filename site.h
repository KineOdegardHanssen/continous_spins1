#ifndef SITE_H
#define SITE_H
#include <iostream>
#include <vector>
#include <bond.h>


class Site
{
public:

    int index; // is this neccessary?
    int lenint; // Length of interaction array.
    int n1, n2, n3;
    double hx, hy, hz, Dix, Diy, Diz; // This would lead to a lot of zero elements, possibly
    double spinx, spiny, spinz;

    // Should I have neighbours here?
    //std::vector<int> bondtype1; // Not here
    //std::vector<int> bondtype2;
    //std::vector<int> bondtype3;
    std::vector<double> siteint; // Array for site interaction
    std::vector<Bond> bonds; // I guess I do have to add a pointer to Bond. But then both classes inherit each other. Weird.
    // Should I have some bondlength parameter?

    // For closed boundary conditions
    int no_of_neighbours_site;

    // Initializer
    // For periodic boundary conditions
    Site(int n, bool sianisotropy, bool magfield, double spinx, double spiny, double spinz, std::vector<double> siteint, std::vector<Bond> bonds);
    // For closed boundary conditions
    Site(int n, int no_of_neighbours_site, bool sianisotropy, bool magfield, double spinx, double spiny, double spinz, std::vector<double> siteint, std::vector<Bond> bonds);
};

#endif // SITE_H
