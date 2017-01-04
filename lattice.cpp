#include "lattice.h"

Lattice::Lattice(int L, bool isotropic, bool sianisotropy, bool magfield, bool dm)
{   // Ta inn chars eller bools og bestemme hvilken funksjon som skal kalles?
    this->L = L;
    // Do these really need to be a part of Lattice?
    // Quite neat, though not neccessary. No, I guess these come in handy.
    this->isotropic = isotropic;
    this->sianisotropy = sianisotropy;
    this->magfield = magfield;
    this->dm = dm;
}

void Lattice::quadratic_helical_initialize()
{   // This one is primarily for testing.
    //N = L*(L+1); // Look this up!
    N = L*L;
    no_of_neighbours = 4;
    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    // Interactions. Should have a way of choosing which terms we look at. Maybe different initialization
    // functions in site? Send in a char for that instead of having all these ones. Quickly get a lot of
    // unneccessary calculations.
    // have some function for doing this:
    double hx = 1;
    double hy = 1;
    double hz = 1;
    double Dix = 1;
    double Diy = 1;
    double Diz = 1;

    // And, in the future, have it in the loop.

    //double J = 1;
    //double Dx = 1;
    //double Dy = 1;
    //double Dz = 1;

    // Move these when neccessary
    std::vector<double> siteint = givethesiteints(Dix, Diy, Diz, hx, hy, hz, sianisotropy, magfield);
    //std::vector<double> bondints = givethebondints(J, Dx, Dy, Dz, isotropic, dm);

    // Alternative implementation: Have Lattice tell which bond belongs to which sites. Initialization
    int bondcounter = 0;
    vecint bondsatsite = vecint(no_of_neighbours);
    vecint sitesatbond = vecint(2);


    // Could have these inside the loop and add randomness.

    for(int n=0; n<N; n++)
    {
        // Finding the neighbours to n
        // NB!: So far, I have only added neighbours that are a distance 1 apart, 2 being the length of the cell.
        // So no linking two cell corner atoms together as of now.
        // This should only be done once. And that is exactly what we are doing.
        // Doing modulo operations, as suggested in Newman & Barkema
        // These neighbours are consistent with the sketch in Newman & Barkema
        int np1 = (n+1)%N;
        int npL = (n+L)%N;
        int nm1 = (n+N-1)%N;
        int nmL = (n+N-L)%N;

        //std::vector<Bond> bonds;

        // Making a lot of bond classes to be added to bonds.
        //bonds.push_back(Bond(n, np1, isotropic, dm, bondints));  // Do I really need to send in n?
        //bonds.push_back(Bond(n, nm1, isotropic, dm, bondints));
        //bonds.push_back(Bond(n, npL, isotropic, dm, bondints));
        //bonds.push_back(Bond(n, nmL, isotropic, dm, bondints));
        // or
        //bonds.push_back(Bond(n, np1, bondints));

        // Is it too nested to make Site inherit Bond? ... Seems fair?
        // Send in bools
        sites.push_back(Site(n, sianisotropy, magfield, spinx, spiny, spinz, siteint));
        // or
        //sites.push_back(Site(n, spinx, spiny, spinz, hx, hy, hz, Dix, Diy, Diz, bonds));
        vecint neighbours;
        neighbours.push_back(np1);
        neighbours.push_back(npL);
        neighbours.push_back(nm1);
        neighbours.push_back(nmL);
        for(int i = 0; i<no_of_neighbours; i++)                // Traversing over all the neighbours
        {
            //cout << "in loop over neighbours" << endl;
            //cout << "n = " << n << endl;
            //cout << "neighbours[i]" << neighbours[i] << endl;

            if(n<neighbours[i])                  // If the bond is not set before
            {
                sitesatbond[0] = n;              // Store the lowest number first
                sitesatbond[1] = neighbours[i];
                sitesofbonds.push_back(sitesatbond);
                bondsatsite[i] = bondcounter;
                //cout << bondcounter << endl;
                bondcounter++;
            }
            else // Well, what to do?
            {
                // If
                // Could periodic BC can mess this up? Probably not, since we start at the low indices...
                for(int j=0; j<bondcounter; j++)  // looping over our bonds so far
                {   // This is probably a bit inefficient...
                    sitesatbond = sitesofbonds[j];
                    // If we have found the bond, store index j in bondsatsite:
                    if(sitesatbond[0]==neighbours[i] && sitesatbond[1]==n)    bondsatsite[i] = j;
                }
            }
        }
        bondsofsites.push_back(bondsatsite);
    }
    setbondinteractions_test();
}

void Lattice::cubic_helical_initialize()
{
    //N = L*(L+1)*(L+1); // Look this up!
    N = L*L*L;
    no_of_neighbours = 6;
    //cout << "Do not use cubib_helical_initialize() yet. Not implemented" << endl;
    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    // Interactions. Should have a way of choosing which terms we look at. Maybe different initialization
    // functions in site? Send in a char for that instead of having all these ones. Quickly get a lot of
    // unneccessary calculations.
    // have some function for doing this:
    double hx = 1;
    double hy = 1;
    double hz = 1;
    double Dix = 1;
    double Diy = 1;
    double Diz = 1;

    // And, in the future, have it in the loop.

    //double J = 1;
    //double Dx = 1;
    //double Dy = 1;
    //double Dz = 1;

    // Move these when neccessary
    std::vector<double> siteint = givethesiteints(Dix, Diy, Diz, hx, hy, hz, sianisotropy, magfield);
    //std::vector<double> bondints = givethebondints(J, Dx, Dy, Dz, isotropic, dm);

    // Alternative implementation: Have Lattice tell which bond belongs to which sites. Initialization
    int bondcounter = 0;
    vecint bondsatsite = vecint(no_of_neighbours);
    vecint sitesatbond = vecint(2);


    // Could have these inside the loop and add randomness.

    for(int n=0; n<N; n++)
    {
        // Finding the neighbours to n
        // NB!: So far, I have only added neighbours that are a distance 1 apart, 2 being the length of the cell.
        // So no linking two cell corner atoms together as of now.
        // This should only be done once. And that is exactly what we are doing.
        // Doing modulo operations, as suggested in Newman & Barkema
        // These seems correct. However, I should probably find a second source to support this.
        int np1 = (n+1)%N;                   // Or should I call a function neighbours_xxxtype()?
        int npL = (n+L)%N;                   // Need to change these here, at least.
        int npL2 = (n+L*L)%N;
        int nm1 = (n+N-1)%N;
        int nmL = (n+N-L)%N;
        int nmL2 = (n+N-L*L)%N;

        //std::vector<Bond> bonds;

        // Making a lot of bond classes to be added to bonds.
        //bonds.push_back(Bond(n, np1, isotropic, dm, bondints));  // Do I really need to send in n?
        //bonds.push_back(Bond(n, nm1, isotropic, dm, bondints));
        //bonds.push_back(Bond(n, npL, isotropic, dm, bondints));
        //bonds.push_back(Bond(n, nmL, isotropic, dm, bondints));
        //bonds.push_back(Bond(n, npL2, isotropic, dm, bondints));
        //bonds.push_back(Bond(n, nmL2, isotropic, dm, bondints));

        // or
        //bonds.push_back(Bond(n, np1, bondints));

        // Is it too nested to make Site inherit Bond? ... Seems fair?
        // Send in bools
        sites.push_back(Site(n, sianisotropy, magfield, spinx, spiny, spinz, siteint)); //, bonds)
        // or
        //sites.push_back(Site(n, spinx, spiny, spinz, hx, hy, hz, Dix, Diy, Diz, bonds));
        vecint neighbours;
        neighbours.push_back(np1);
        neighbours.push_back(npL);
        neighbours.push_back(npL2);
        neighbours.push_back(nm1);
        neighbours.push_back(nmL);
        neighbours.push_back(nmL2);
        for(int i = 0; i<no_of_neighbours; i++)                // Traversing over all the neighbours
        {
            //cout << "in loop over neighbours" << endl;
            //cout << "n = " << n << endl;
            //cout << "neighbours[i]" << neighbours[i] << endl;

            if(n<neighbours[i])                  // If the bond is not set before
            {
                sitesatbond[0] = n;              // Store the lowest number first
                sitesatbond[1] = neighbours[i];
                sitesofbonds.push_back(sitesatbond);
                bondsatsite[i] = bondcounter;
                //cout << bondcounter << endl;
                bondcounter++;
            }
            else // Well, what to do?
            {
                // If
                // Could periodic BC can mess this up? Probably not, since we start at the low indices...
                for(int j=0; j<bondcounter; j++)  // looping over our bonds so far
                {   // This is probably a bit inefficient...
                    sitesatbond = sitesofbonds[j];
                    // If we have found the bond, store index j in bondsatsite:
                    if(sitesatbond[0]==neighbours[i] && sitesatbond[1]==n)    bondsatsite[i] = j;
                }
            }
        }
        bondsofsites.push_back(bondsatsite);
    }
    setbondinteractions_test();
}

void Lattice::fcc_helical_initialize()
{   // Should include something saying how the parameters are set.
    //N = L*(L+1)*(L+1); // Look this up!
    N = L*L*L;
    no_of_neighbours = 12;
    cout << "In fcc_helical_initialize" << endl;
    // Setting up the sites
    // We set up the matrix by having all spins in the same direction. Or maybe draw at random?
    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    // Interactions. Should have a way of choosing which terms we look at. Maybe different initialization
    // functions in site? Send in a char for that instead of having all these ones. Quickly get a lot of
    // unneccessary calculations.
    // have some function for doing this:
    double hx = 1;
    double hy = 1;
    double hz = 1;
    double Dix = 1;
    double Diy = 1;
    double Diz = 1;

    // And, in the future, have it in the loop.

    //double J = 1;
    //double Dx = 1;
    //double Dy = 1;
    //double Dz = 1;

    // Move these when neccessary
    std::vector<double> siteint = givethesiteints(Dix, Diy, Diz, hx, hy, hz, sianisotropy, magfield);
    //std::vector<double> bondints = givethebondints(J, Dx, Dy, Dz, isotropic, dm);

    // Alternative implementation: Have Lattice tell which bond belongs to which sites. Initialization
    int bondcounter = 0;
    vecint bondsatsite = vecint(12);
    vecint sitesatbond = vecint(2);


    // Could have these inside the loop and add randomness.

    for(int n=0; n<N; n++)
    {
        // Finding the neighbours to n
        // NB!: So far, I have only added neighbours that are a distance 1 apart, 2 being the length of the cell.
        // So no linking two cell corner atoms together as of now.
        // This should only be done once. And that is exactly what we are doing.
        // Doing modulo operations, as suggested in Newman & Barkema
        // Should probably verify this bit
        int np1 = (n+1)%N;
        int npL = (n+L)%N;
        int npL2 = (n+L*L)%N;
        int npLm1 = (n+L-1)%N;
        int npL2m1 = (n+L*L-1)%N;
        int npL2mL = (n+L*L-L)%N;
        int nm1 = (n+N-1)%N; // Seems to work for the 5x5 case, at least
        int nmL = (n+N-L)%N;
        int nmL2 = (n+N-L*L)%N;
        int nmLm1 = (n+N-L+1)%N;
        int nmL2m1 = (n+N-L*L+1)%N;
        int nmL2mL = (n+N-L*L+L)%N;

        //std::vector<Bond> bonds;

        // Making a lot of bond classes to be added to bonds.
        /*
        bonds.push_back(Bond(n, np1, isotropic, dm, bondints));  // Do I really need to send in n?
        bonds.push_back(Bond(n, nm1, isotropic, dm, bondints));
        bonds.push_back(Bond(n, npL, isotropic, dm, bondints));
        bonds.push_back(Bond(n, nmL, isotropic, dm, bondints));
        bonds.push_back(Bond(n, npL2, isotropic, dm, bondints));
        bonds.push_back(Bond(n, nmL2, isotropic, dm, bondints));
        bonds.push_back(Bond(n, npLm1, isotropic, dm, bondints));
        bonds.push_back(Bond(n, nmLm1, isotropic, dm, bondints));
        bonds.push_back(Bond(n, npL2m1, isotropic, dm, bondints));
        bonds.push_back(Bond(n, nmL2m1, isotropic, dm, bondints));
        bonds.push_back(Bond(n, npL2mL, isotropic, dm, bondints));
        bonds.push_back(Bond(n, nmL2mL, isotropic, dm, bondints));
        */
        // or
        //bonds.push_back(Bond(n, np1, bondints));

        // Is it too nested to make Site inherit Bond? ... Seems fair?
        // Send in bools
        sites.push_back(Site(n, sianisotropy, magfield, spinx, spiny, spinz, siteint)); // , bonds))
        // or
        //sites.push_back(Site(n, spinx, spiny, spinz, hx, hy, hz, Dix, Diy, Diz, bonds));
        vecint neighbours;
        neighbours.push_back(np1);
        neighbours.push_back(npL);
        neighbours.push_back(npL2);
        neighbours.push_back(npLm1);
        neighbours.push_back(npL2m1);
        neighbours.push_back(npL2mL);
        neighbours.push_back(nm1);
        neighbours.push_back(nmL);
        neighbours.push_back(nmL2);
        neighbours.push_back(nmLm1);
        neighbours.push_back(nmL2m1);
        neighbours.push_back(nmL2mL);
        for(int i = 0; i<no_of_neighbours; i++)                // Traversing over all the neighbours
        {
            //cout << "in loop over neighbours" << endl;
            //cout << "n = " << n << endl;
            //cout << "neighbours[i]" << neighbours[i] << endl;

            if(n<neighbours[i])                  // If the bond is not set before
            {
                sitesatbond[0] = n;              // Store the lowest number first
                sitesatbond[1] = neighbours[i];
                sitesofbonds.push_back(sitesatbond);
                bondsatsite[i] = bondcounter;
                //cout << bondcounter << endl;
                bondcounter++;
            }
            else // Well, what to do?
            {
                // If
                // Could periodic BC can mess this up? Probably not, since we start at the low indices...
                for(int j=0; j<bondcounter; j++)  // looping over our bonds so far
                {   // This is probably a bit inefficient...
                    sitesatbond = sitesofbonds[j];
                    // If we have found the bond, store index j in bondsatsite:
                    if(sitesatbond[0]==neighbours[i] && sitesatbond[1]==n)    bondsatsite[i] = j;
                }
            }
        }
        bondsofsites.push_back(bondsatsite);
    }
    setbondinteractions_test();
}

void Lattice::setbondinteractions_test()
{   // Should make this more flexible. Feeding in some chars or something.
    // Alternatively, have that as class variables set by the constructor.
    int limit = 0.5*no_of_neighbours*N; // Floor or something?
    vector<double> Js = vector<double>(limit);
    vector<double> Dxes = vector<double>(limit);
    vector<double> Dys = vector<double>(limit);
    vector<double> Dzs = vector<double>(limit);

    for(int i=0; i<limit; i++)
    {
        Js[i]   = 1;
        Dxes[i] = 1;
        Dys[i]  = 1;
        Dzs[i]  = 1;
    }
    bonds = Bond(isotropic, dm, Js, Dxes, Dys, Dzs);
}

std::vector<double> Lattice::givethesiteints(double Dix, double Diy, double Diz, double hx, double hy, double hz, bool sianisotropy, bool magfield)
{
    if(sianisotropy && magfield)
    {
        std::vector<double> siteint = std::vector<double>(6);
        siteint[0] = Dix;
        siteint[1] = Diy;
        siteint[2] = Diz;
        siteint[3] = hx;
        siteint[4] = hy;
        siteint[5] = hz;
        return siteint;
    }
    else if(sianisotropy && (!magfield))
    {
        std::vector<double> siteint = std::vector<double>(3);
        siteint[0] = Dix;
        siteint[1] = Diy;
        siteint[2] = Diz;
        return siteint;
    }
    else if((!sianisotropy) && magfield)
    {
        std::vector<double> siteint = std::vector<double>(3);
        siteint[0] = hx;
        siteint[1] = hy;
        siteint[2] = hz;
        return siteint;
    }
    else
    {
        std::vector<double> siteint = std::vector<double>(1);
        siteint[0] = 0;
        return siteint;
    }
}

std::vector<double> Lattice::givethebondints(double J, double Dx, double Dy, double Dz, bool isotropic, bool dm)
{
    if(isotropic && dm)
    {
        std::vector<double> bondints = std::vector<double>(4);
        bondints[0] = J;
        bondints[1] = Dx;
        bondints[2] = Dy;
        bondints[3] = Dz;
        return bondints;
    }
    else if(isotropic && (!dm))
    {
        std::vector<double> bondints = std::vector<double>(1);
        bondints[0] = J;
        return bondints;
    }
    else if((!isotropic) && dm)
    {
        std::vector<double> bondints = std::vector<double>(3);
        bondints[0] = Dx;
        bondints[1] = Dy;
        bondints[2] = Dz;
        return bondints;
    }
    else
    {
        std::vector<double> bondints = std::vector<double>(1);
        bondints[0] = 0;
        return bondints;
    }
}
