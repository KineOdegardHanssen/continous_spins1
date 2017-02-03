#include "lattice.h"

Lattice::Lattice()
{

}

Lattice::Lattice(int L, bool isotropic, bool sianisotropy, bool magfield, bool dm)
{
    this->L = L;
    this->isotropic = isotropic;
    this->sianisotropy = sianisotropy;
    this->magfield = magfield;
    this->dm = dm;

    notperiodic = false; // Default value. Changes if we choose a lattice with closed boundary conditions
}

void Lattice::chain_open_initialize()
{
    notperiodic = true;
    bool randomspins = false;  // Adding the option to change the spins
    cout << "NB! Bool notperiodic changed to true. Now operating with closed BCs!" << endl;
    dim = 1;
    N = L;
    no_of_neighbours = 2; // For the most part

    dimlengths = vector<int>(1);
    dimlengths[0] = L;
    a1 = vector<double>(1);
    a1[0] = 1;
    b1 = vector<double>(1);
    b1[0] = 2*M_PI;

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

    double J =  1; // As in Ising model work
    double Dx = 1;
    double Dy = 1;
    double Dz = 1;

    // Move these when neccessary
    std::vector<double> siteint = givethesiteints(Dix, Diy, Diz, hx, hy, hz, sianisotropy, magfield);
    std::vector<double> bondints = givethebondints(J, Dx, Dy, Dz, isotropic, dm);

    for(int n=0; n<N; n++)
    {
        // Finding the neighbours to n
        int no_of_neighbours_site = 0;
        int np1 = (n+1)%N;
        int nm1 = (n+N-1)%N;

        std::vector<Bond> bonds;

        if(randomspins)
        {
            double u = ran2(&seed);
            double v = ran2(&seed);

            double theta = acos(1.0-2.0*u);
            double phi = 2.0*M_PI*v;

            double sintheta = sin(theta);
            spinx = sintheta*cos(phi);
            spiny = sintheta*sin(phi);
            spinz = cos(theta);
        }


        // Excluding periodic neighbours
        if(np1>n)
        {
            bonds.push_back(Bond(n, np1, isotropic, dm, bondints));
            no_of_neighbours_site++;
        }
        if(nm1<n)
        {
            bonds.push_back(Bond(n, nm1, isotropic, dm, bondints));
            no_of_neighbours_site++;
        }


        //cout << "Spin " << n << ", neighbours: " << np1 << ", " << nm1 << endl;

        // Is it too nested to make Site inherit Bond? ... Seems fair?
        // Send in bools
        sites.push_back(Site(n, no_of_neighbours_site, sianisotropy, magfield, spinx, spiny, spinz, siteint, bonds));
        // or
        //sites.push_back(Site(n, spinx, spiny, spinz, hx, hy, hz, Dix, Diy, Diz, bonds));
    }
}


void Lattice::chain_periodic_initialize()
{
    dim = 1;
    N = L;
    no_of_neighbours = 2;

    dimlengths = vector<int>(1);
    dimlengths[0] = L;
    a1 = vector<double>(1);
    a1[0] = 1;
    b1 = vector<double>(1);
    b1[0] = 2*M_PI;

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

    double J =  1;
    double Dx = 1;
    double Dy = 1;
    double Dz = 1;

    // Move these when neccessary
    std::vector<double> siteint = givethesiteints(Dix, Diy, Diz, hx, hy, hz, sianisotropy, magfield);
    std::vector<double> bondints = givethebondints(J, Dx, Dy, Dz, isotropic, dm);

    for(int n=0; n<N; n++)
    {
        // Finding the neighbours to n
        // This should only be done once. And that is exactly what we are doing.
        // Doing modulo operations, as suggested in Newman & Barkema
        // These neighbours are consistent with the sketch in Newman & Barkema
        int np1 = (n+1)%N;
        int nm1 = (n+N-1)%N;

        cout << "Spin " << n << ", neighbours: " << np1 << ", " << nm1 << endl;

        std::vector<Bond> bonds;

        // Making a lot of bond classes to be added to bonds.
        bonds.push_back(Bond(n, np1, isotropic, dm, bondints));  // Do I really need to send in n?
        bonds.push_back(Bond(n, nm1, isotropic, dm, bondints));

        // or
        //bonds.push_back(Bond(n, np1, bondints));

        // Is it too nested to make Site inherit Bond? ... Seems fair?
        // Send in bools
        sites.push_back(Site(n, sianisotropy, magfield, spinx, spiny, spinz, siteint, bonds));
        // or
        //sites.push_back(Site(n, spinx, spiny, spinz, hx, hy, hz, Dix, Diy, Diz, bonds));
    }
}

void Lattice::quadratic_helical_initialize()
{   // This one is primarily for testing.
    //N = L*(L+1); // Look this up!
    dim = 2;
    N = L*L;
    no_of_neighbours = 4;

    // No of particles in each direction
    dimlengths = vector<int>(dim);
    dimlengths[0] = L;
    dimlengths[0] = L;

    // Lattice vectors
    a1 = vector<double>(2);
    a2 = vector<double>(2);
    a1[0] = 1;
    a1[1] = 0;
    a2[0] = 0;
    a2[1] = 1;
    // Reciprocal lattice vectors
    b1 = vector<double>(2);
    b2 = vector<double>(2);
    b1[0] = 2*M_PI;
    b1[1] = 0;
    b2[0] = 0;
    b2[1] = 2*M_PI;


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

    double J =  1; // As in Ising model work
    double Dx = 1;
    double Dy = 1;
    double Dz = 1;

    // Move these when neccessary
    std::vector<double> siteint = givethesiteints(Dix, Diy, Diz, hx, hy, hz, sianisotropy, magfield);
    std::vector<double> bondints = givethebondints(J, Dx, Dy, Dz, isotropic, dm);

    std::vector<double> position_n = std::vector<double>(2);
    std::vector<int> coord_n = std::vector<int>(2);

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

        std::vector<Bond> bonds;

        // Making a lot of bond classes to be added to bonds.
        bonds.push_back(Bond(n, np1, isotropic, dm, bondints));  // Do I really need to send in n?
        bonds.push_back(Bond(n, nm1, isotropic, dm, bondints));
        bonds.push_back(Bond(n, npL, isotropic, dm, bondints));
        bonds.push_back(Bond(n, nmL, isotropic, dm, bondints));
        // or
        //bonds.push_back(Bond(n, np1, bondints));

        // Is it too nested to make Site inherit Bond? ... Seems fair?
        // Send in bools
        sites.push_back(Site(n, sianisotropy, magfield, spinx, spiny, spinz, siteint, bonds));
        // or
        //sites.push_back(Site(n, spinx, spiny, spinz, hx, hy, hz, Dix, Diy, Diz, bonds));

        position_n[0] = 1.0*((int)n%L); // n1
        position_n[1] = 1.0*((int)n/L); // n2. Should allow for grid length a.

        coord_n[0] = ((int)n%L); // n1
        coord_n[1] = ((int)n/L); // n2

        sitepositions.push_back(position_n);
        sitecoordinates.push_back(coord_n);

    } // End of loop over all sites ( = loop over n)
}

void Lattice::cubic_helical_initialize()
{
    //N = L*(L+1)*(L+1); // Look this up!
    dim = 3;
    N = L*L*L;
    no_of_neighbours = 6;

    // No of particles in each direction
    dimlengths    = vector<int>(dim);
    dimlengths[0] = L;
    dimlengths[1] = L;
    dimlengths[2] = L;

    // Lattice vectors
    a1 = vector<double>(3); // Should they be double?
    a2 = vector<double>(3);
    a3 = vector<double>(3);
    a1[0] = 1;
    a1[1] = 0;
    a1[2] = 0;
    a2[0] = 0;
    a2[1] = 1;
    a2[2] = 0;
    a3[0] = 0;
    a3[1] = 0;
    a3[2] = 1;
    // Reciprocal lattice vectors
    b1 = vector<double>(3); // Should they be double?
    b2 = vector<double>(3);
    b3 = vector<double>(3);
    b1[0] = 2*M_PI;
    b1[1] = 0;
    b1[2] = 0;
    b2[0] = 0;
    b2[1] = 2*M_PI;
    b2[2] = 0;
    b3[0] = 0;
    b3[1] = 0;
    b3[2] = 2*M_PI;

    // Temporary fix. May want to send in L1, L2 and L3 to Lattice and deal with a 'rectangular' crystal
    int L1 = L;
    int L2 = L;

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

    double J = 1;
    double Dx = 1;
    double Dy = 1;
    double Dz = 1;

    // Move these when neccessary
    std::vector<double> siteint = givethesiteints(Dix, Diy, Diz, hx, hy, hz, sianisotropy, magfield);
    std::vector<double> bondints = givethebondints(J, Dx, Dy, Dz, isotropic, dm);

    std::vector<double> position_n = std::vector<double>(3);
    std::vector<int> coord_n = std::vector<int>(3);

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

        std::vector<Bond> bonds;

        // Making a lot of bond classes to be added to bonds.
        bonds.push_back(Bond(n, np1, isotropic, dm, bondints));  // Do I really need to send in n?
        bonds.push_back(Bond(n, nm1, isotropic, dm, bondints));
        bonds.push_back(Bond(n, npL, isotropic, dm, bondints));
        bonds.push_back(Bond(n, nmL, isotropic, dm, bondints));
        bonds.push_back(Bond(n, npL2, isotropic, dm, bondints));
        bonds.push_back(Bond(n, nmL2, isotropic, dm, bondints));

        // or
        //bonds.push_back(Bond(n, np1, bondints));

        // Is it too nested to make Site inherit Bond? ... Seems fair?
        // Send in bools
        sites.push_back(Site(n, sianisotropy, magfield, spinx, spiny, spinz, siteint, bonds));
        // or
        //sites.push_back(Site(n, spinx, spiny, spinz, hx, hy, hz, Dix, Diy, Diz, bonds));

        // Giving the position
        int n1 = n%L1;
        int n2 = n/L1 - n/(L1*L2)*L2;
        int n3 = n/(L1*L2);

        position_n[0] = 1.0*n1;     // Could possibly multiply by grid length a
        position_n[1] = 1.0*n2;
        position_n[2] = 1.0*n3;

        coord_n[0] = n1;     // Could possibly multiply by grid length a
        coord_n[1] = n2;
        coord_n[2] = n3;

        sitepositions.push_back(position_n);
        sitecoordinates.push_back(coord_n);
    } // End loop over n (all sites)

}

void Lattice::fcc_helical_initialize()
{
    bool DEBUG = true;
    // Should include something saying how the parameters are set.
    //N = L*(L+1)*(L+1); // Look this up!
    dim = 3;
    N = L*L*L;
    no_of_neighbours = 12;

    // No of particles in each direction
    dimlengths    = vector<int>(dim);
    dimlengths[0] = L;
    dimlengths[1] = L;
    dimlengths[2] = L;


    // Lattice vectors
    a1 = vector<double>(3); // Should they be double?
    a2 = vector<double>(3);
    a3 = vector<double>(3);
    a1[0] = 0.5;
    a1[1] = 0.5;
    a1[2] = 0;
    a2[0] = 0;
    a2[1] = 0.5;
    a2[2] = 0.5;
    a3[0] = 0.5;
    a3[1] = 0;
    a3[2] = 0.5;
    // Reciprocal lattice vectors
    b1 = vector<double>(3); // Should they be double?
    b2 = vector<double>(3);
    b3 = vector<double>(3);
    b1[0] = M_PI;
    b1[1] = M_PI;
    b1[2] = 0;
    b2[0] = 0;
    b2[1] = M_PI;
    b2[2] = M_PI;
    b3[0] = M_PI;
    b3[1] = 0;
    b3[2] = M_PI;

    cout << "In fcc_helical_initialize" << endl;

    // Temporary fix. May want to send L1, L2 and L3 into Lattice and deal with a 'rectangular' lattice
    int L1 = L;
    int L2 = L;

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
    // magfield
    double hx = 1;
    double hy = 1;
    double hz = 1;
    //sianisotropy
    double Dix = 1;
    double Diy = 1;
    double Diz = 1;
    // And, in the future, have it in the loop.
    double J =  1;
    double Dx = 1;
    double Dy = 1;
    double Dz = 1;

    // Move these when neccessary
    std::vector<double> siteint = givethesiteints(Dix, Diy, Diz, hx, hy, hz, sianisotropy, magfield);
    std::vector<double> bondints = givethebondints(J, Dx, Dy, Dz, isotropic, dm);

    std::vector<double> position_n = std::vector<double>(3);
    std::vector<int> coord_n = std::vector<int>(3);

    // Could have these inside the loop and add randomness.

    if(DEBUG)    cout << "Now entering the loop" << endl;
    for(int n=0; n<N; n++)
    {
        if(DEBUG)    cout << "In loop, n = " << n << endl;
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

        if(DEBUG)
        {
            cout << "np1" << " " << np1 << endl;
            cout << "npL" << " " << npL << endl;
            cout << "npL2" << " " << npL2 << endl;
            cout << "npLm1" << " " << npLm1 << endl;
            cout << "npL2m1" << " " << npL2m1 << endl;
            cout << "npL2mL" << " " << npL2mL << endl;
            cout << "nm1" << " " << nm1 << endl;
            cout << "nmL" << " " << nmL << endl;
            cout << "nmL2" << " " << nmL2 << endl;
            cout << "nmLm1" << " " << nmLm1 << endl;
            cout << "nmL2m1" << " " << nmL2m1 << endl;
            cout << "nmL2mL" << " " << nmL2mL << endl;
        }

        std::vector<Bond> bonds;

        if(DEBUG)    cout << "Setting the bonds" << endl;
        // Making a lot of bond classes to be added to bonds.
        bonds.push_back(Bond(n, np1, isotropic, dm, bondints));  // Do I really need to send in n?
        if(DEBUG)    cout << "Bond 1 done" << endl;
        bonds.push_back(Bond(n, nm1, isotropic, dm, bondints));
        if(DEBUG)    cout << "Bond 2 done" << endl;
        bonds.push_back(Bond(n, npL, isotropic, dm, bondints));
        if(DEBUG)    cout << "Bond 3 done" << endl;
        bonds.push_back(Bond(n, nmL, isotropic, dm, bondints));
        if(DEBUG)    cout << "Bond 4 done" << endl;
        bonds.push_back(Bond(n, npL2, isotropic, dm, bondints));
        if(DEBUG)    cout << "Bond 5 done" << endl;
        bonds.push_back(Bond(n, nmL2, isotropic, dm, bondints));
        if(DEBUG)    cout << "Bond 6 done" << endl;
        bonds.push_back(Bond(n, npLm1, isotropic, dm, bondints));
        if(DEBUG)    cout << "Bond 7 done" << endl;
        bonds.push_back(Bond(n, nmLm1, isotropic, dm, bondints));
        if(DEBUG)    cout << "Bond 8 done" << endl;
        bonds.push_back(Bond(n, npL2m1, isotropic, dm, bondints));
        if(DEBUG)    cout << "Bond 9 done" << endl;
        bonds.push_back(Bond(n, nmL2m1, isotropic, dm, bondints));
        if(DEBUG)    cout << "Bond 10 done" << endl;
        bonds.push_back(Bond(n, npL2mL, isotropic, dm, bondints));
        if(DEBUG)    cout << "Bond 11 done" << endl;
        bonds.push_back(Bond(n, nmL2mL, isotropic, dm, bondints));
        if(DEBUG)    cout << "Bond 12 done" << endl;

        // or
        //bonds.push_back(Bond(n, np1, bondints));

        // Is it too nested to make Site inherit Bond? ... Seems fair?
        // Send in bools
        cout << "Done setting the bonds. Setting the sites" << endl;
        sites.push_back(Site(n, sianisotropy, magfield, spinx, spiny, spinz, siteint, bonds));
        // or
        //sites.push_back(Site(n, spinx, spiny, spinz, hx, hy, hz, Dix, Diy, Diz, bonds));

        // Listing the positions of the spins:

        //  /*
        if(DEBUG)    cout << "Giving the position of the site in the fcc" << endl;
        // Giving the position of the fcc
        int n1 = n%L1;
        int n2 = (int)n/L1 - (int)n/(L1*L2)*L2;
        int n3 = n/(L1*L2);

        double xpos = 0.5*(n1+n3);  // Could possibly include the grid length a
        double ypos = 0.5*(n1+n2);
        double zpos = 0.5*(n2+n3);

        cout << "[" << xpos << "," << ypos << "," << zpos << "]" << endl;

        position_n[0] = xpos;
        position_n[1] = ypos;
        position_n[2] = zpos;

        coord_n[0] = n1;     // Have these here in case I want to check
        coord_n[1] = n2;
        coord_n[2] = n3;

        sitepositions.push_back(position_n);
        sitecoordinates.push_back(coord_n);
        // */

    }
    cout << "Done with fcc_helical_initialize" << endl;

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
