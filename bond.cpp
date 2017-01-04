#include "bond.h"

Bond::Bond()
{
}

Bond::Bond(bool isotropic, bool dm, std::vector<double> tJs, std::vector<double> tDxes, std::vector<double> tDys, std::vector<double> tDzs)
{
    if(isotropic)    Js = tJs;
    if(dm)
    {
        Dxes = tDxes;
        Dys = tDys;
        Dzs = tDzs;
    }
}


/*
Bond::Bond(bool isotropic, bool dm, std::vector bondints)
{
    this->bondints = bondints; // Do I really need this
    // Or function feed bondints

    this->siteindex1 = siteindex1;
    this->siteindex2 = siteindex2;

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
}
*/

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
