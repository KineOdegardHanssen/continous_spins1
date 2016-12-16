#ifndef SITE_H
#define SITE_H
#include <iostream>
#include <vector>
#include <bond.h>


class Site
{
public:
    Site();
    int index; // is this neccessary?
    double hx, hy, hz, Dix, Diy, Diz;
    double spinx, spiny, spinz;

    // Should I have neighbours here?
    //std::vector<int> bondtype1; // Not here
    //std::vector<int> bondtype2;
    //std::vector<int> bondtype3;
    //std::vector<double> spin;   // maybe store in each direction instead
    //std::vector<int> bondindexes; // Or some sort of pointer...
    std::vector<Bond> bonds; // I guess I do have to add a pointer to Bond. But then both classes inherit each other. Weird.

};

#endif // SITE_H
