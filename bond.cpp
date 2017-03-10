#include "bond.h"

Bond::Bond()
{
}

// Standard initializer
Bond::Bond(int siteindex1, int siteindex2, bool isotropic, bool dm, bool increasing, std::vector<double> bondints)
{
    //std::cout << "Initializing instance of class Bond" << std::endl;
    this->bondints = bondints; // Do I really need this?
    // Or function feed bondints

    //std::cout << "Bondints set" << std::endl;

    this->siteindex1 = siteindex1;
    this->siteindex2 = siteindex2;

    this->increasing = increasing;

    //std::cout << "Siteindex set" << std::endl;

    if(isotropic==true)
    {
        J = bondints[0];
        if(dm==true)
        {
            Dx = bondints[1];
            Dy = bondints[2];
            Dz = bondints[3];
        }
    }
    else
    {
        if(dm==true)
        {
            Dx = bondints[0];
            Dy = bondints[1];
            Dz = bondints[2];
        }
    }
    //std::cout << "Interactions unpacked, exiting Bond" << std::endl;
}

// Bond for strong interactions. Should maybe change?
Bond::Bond(int siteindex1, int siteindex2, bool strong, bool isotropic, bool dm, bool increasing, std::vector<double> bondints)
{
    // bool strong should either cause J or the D's to be multiplied by something, or be a class variable to be
    // tested. Probably more efficient to implement the former.
    //std::cout << "Initializing instance of class Bond" << std::endl;
    this->bondints = bondints; // Do I really need this?
    // Or function feed bondints

    //std::cout << "Bondints set" << std::endl;

    this->siteindex1 = siteindex1;
    this->siteindex2 = siteindex2;

    this->increasing = increasing;

    //std::cout << "Siteindex set" << std::endl;

    if(isotropic==true)
    {
        //double strength = 1; // Default strength
        //if(strong)    strength = a; input variable?
        //J = a*bondints;   // This is general enough, isn't it? We only want to look at one lattice at a time
        J = bondints[0];
        if(dm==true)
        {
            // Dx = a*bondints[1];
            // Dy = a*bondints[1];
            // Dz = a*bondints[1];
            Dx = bondints[1];
            Dy = bondints[2];
            Dz = bondints[3];
        }
    }
    else
    {
        if(dm==true)
        {
            Dx = bondints[0];
            Dy = bondints[1];
            Dz = bondints[2];
        }
    }
    //std::cout << "Interactions unpacked, exiting Bond" << std::endl;
}

// Allowing for interactions of different strengths in different directions
Bond::Bond(int siteindex1, int siteindex2, double J, bool increasing, std::vector<double> bondints)
{
    this->J = J;     //
    this->bondints = bondints;
    this->siteindex1 = siteindex1;
    this->siteindex2 = siteindex2;
    this->increasing = increasing;

    Dx = bondints[0]; // Should have an option to not include all these. But I guess it is not neccessary...
    Dy = bondints[1]; //
    Dz = bondints[2]; //
}

// Initializing next nearest neighbour Bond
Bond::Bond(int siteindex1, int siteindex2, double J, bool increasing)
{
    this->siteindex1 = siteindex1;
    this->siteindex2 = siteindex2;
    this->J = J;
    this->increasing = increasing; // Is this neccessary here?
}
