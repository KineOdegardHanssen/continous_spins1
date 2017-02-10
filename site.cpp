#include "site.h"

Site::Site(int n, bool sianisotropy, bool magfield, double spinx, double spiny, double spinz, std::vector<double> siteint, std::vector<Bond> bonds)
{
    index = n;
    this->spinx = spinx;   // Have spinx in a new class, State?s
    this->spiny = spiny;
    this->spinz = spinz;
    this->siteint = siteint;
    this->bonds = bonds;

    // Only including the relevant terms.
    if(sianisotropy==true)
    {
        Dix = siteint[0];
        Diy = siteint[1];
        Diz = siteint[2];
        if(magfield==true)
        {
            hx = siteint[3];
            hy = siteint[4];
            hz = siteint[5];
        }
    }
    else
    {
        if(magfield==true)
        {
            hx = siteint[0];
            hy = siteint[1];
            hz = siteint[2];
        }
    }
}

Site::Site(int n, int no_of_neighbours_site, bool sianisotropy, bool magfield, double spinx, double spiny, double spinz, std::vector<double> siteint, std::vector<Bond> bonds)
{
    index = n;
    this->no_of_neighbours_site = no_of_neighbours_site;
    this->spinx = spinx;   // Have spinx in a new class, States?
    this->spiny = spiny;
    this->spinz = spinz;
    this->siteint = siteint;
    this->bonds = bonds;

    // Only including the relevant terms.
    if(sianisotropy==true)
    {
        Dix = siteint[0];
        Diy = siteint[1];
        Diz = siteint[2];
        if(magfield==true)
        {
            hx = siteint[3];
            hy = siteint[4];
            hz = siteint[5];
        }
    }
    else
    {
        if(magfield==true)
        {
            hx = siteint[0];
            hy = siteint[1];
            hz = siteint[2];
        }
    }
}




Site::Site(int n, bool sianisotropy, bool magfield, double spinx, double spiny, double spinz, std::vector<double> siteint, std::vector<Bond> bonds, std::vector<Bond> nextnearesty, std::vector<Bond> nextnearestz)
{
    index = n;
    this->spinx = spinx;   // Have spinx in a new class, States?
    this->spiny = spiny;
    this->spinz = spinz;
    this->siteint = siteint;
    this->bonds = bonds;
    this->nextnearesty = nextnearesty;
    this->nextnearestz = nextnearestz;

    // Only including the relevant terms.
    if(sianisotropy==true)
    {
        Dix = siteint[0];
        Diy = siteint[1];
        Diz = siteint[2];
        if(magfield==true)
        {
            hx = siteint[3];
            hy = siteint[4];
            hz = siteint[5];
        }
    }
    else
    {
        if(magfield==true)
        {
            hx = siteint[0];
            hy = siteint[1];
            hz = siteint[2];
        }
    }
}
