#include "bond.h"

Bond::Bond()
{
}

Bond::Bond(int siteindex1, int siteindex2, bool isotropic, bool dm, std::vector<double> bondints)
{
    //std::cout << "Initializing instance of class Bond" << std::endl;
    this->bondints = bondints; // Do I really need this?
    // Or function feed bondints

    //std::cout << "Bondints set" << std::endl;

    this->siteindex1 = siteindex1;
    this->siteindex2 = siteindex2;

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


Bond::Bond(int siteindex1, int siteindex2, bool strong, bool isotropic, bool dm, std::vector<double> bondints)
{
    // bool strong should either cause J or the D's to be multiplied by something, or be a class variable to be
    // tested. Probably more efficient to implement the former.
    //std::cout << "Initializing instance of class Bond" << std::endl;
    this->bondints = bondints; // Do I really need this?
    // Or function feed bondints

    //std::cout << "Bondints set" << std::endl;

    this->siteindex1 = siteindex1;
    this->siteindex2 = siteindex2;

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

/*
Bond::Bond(int siteindex1, int siteindex2, double J, double Dx, double Dy, double Dz)
{ //or std::vector<int> site1indexvec, std::vector<int> site2indexvec)
    this->J = J;
    this->Dx = Dx;
    this->Dy = Dy;
    this->Dz = Dz;
    // Or function feed bondints

    this->siteindex1 = siteindex1;
    this->siteindex2 = siteindex2;

    //this->site1indexvec=site1indexvec; // Or do something else?
    //this->site2indexvec=site2indexvec; // Have as

}
*/
