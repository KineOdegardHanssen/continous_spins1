#include "bond.h"

Bond::Bond()
{
}

Bond::Bond(double J, double Dx, double Dy, double Dz, std::vector site1indexvec, std::vector site2indexvec)
{
    this->J = J;
    this->Dx = Dx;
    this->Dy = Dy;
    this->Dz = Dz;

    this->site1indexvec=site1indexvec;
    this->site2indexvec=site2indexvec;
}

