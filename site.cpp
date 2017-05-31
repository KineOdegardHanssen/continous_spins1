#include "site.h"

Site::Site(int n, double spinx, double spiny, double spinz, std::vector<Bond> bonds)
{
    index = n;
    this->spinx = spinx;   // Have spinx in a new class, State?s
    this->spiny = spiny;
    this->spinz = spinz;
    this->bonds = bonds;
}

Site::Site(int n, int no_of_neighbours_site, double spinx, double spiny, double spinz, std::vector<Bond> bonds)
{
    index = n;
    this->no_of_neighbours_site = no_of_neighbours_site;
    this->spinx = spinx;   // Have spinx in a new class, States?
    this->spiny = spiny;
    this->spinz = spinz;
    this->bonds = bonds;
}


Site::Site(int n, double spinx, double spiny, double spinz, std::vector<Bond> bonds, std::vector<Bond> nextnearesty, std::vector<Bond> nextnearestz)
{
    index = n;
    this->spinx = spinx;   // Have spinx in a new class, States?
    this->spiny = spiny;
    this->spinz = spinz;
    this->bonds = bonds;
    this->nextnearesty = nextnearesty;
    this->nextnearestz = nextnearestz;
}
