#include "site.h"

Site::Site()
{
}


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
