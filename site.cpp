#include "site.h"

Site::Site(int n, int L1, int L2, bool sianisotropy, bool magfield, double spinx, double spiny, double spinz, std::vector<double> siteint)
{
    index = n;
    this->spinx = spinx;   // Have spinx in a new class, State?s
    this->spiny = spiny;
    this->spinz = spinz;
    this->siteint = siteint;
    //this->bonds = bonds;

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
    // Giving the position
    n1 = n%L1;
    n2 = n/L1 - n/(L1*L2)*L2;
    n3 = n/(L1*L2);

    xpos = 0.5*(n1+n3);
    ypos = 0.5*(n1+n2);
    zpos = 0.5*(n2+n3);

}




/*// This is probably too inefficient
Site::Site(int n, double hx, double hy, double hz, double Dix, double Diy, double Diz, std::vector<Bond> bonds)
{
    index = n;
    this->hx = hx;
    this->hy = hy;
    this->hz = hz;
    this->Dix= Dix;
    this->Diy= Diy;
    this->Diz= Diz;
    this->bonds = bonds;
}
*/
