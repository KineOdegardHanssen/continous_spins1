#include "lattice.h"

Lattice::Lattice()
{

}

Lattice::Lattice(int L, bool isotropic, bool sianisotropy, bool magfield, bool dm)
{
    this->L = L;
    this->L1 = L;
    this->L2 = L;
    this->L3 = L;
    this->isotropic = isotropic;
    this->sianisotropy = sianisotropy;
    this->magfield = magfield;
    this->dm = dm;

    notperiodic = false; // Default value. Changes if we choose a lattice with closed boundary conditions
    systemstrengthsgiven = false;
    extended = false;
    seed = 23;
    cout << "In L Lattice constructor. L1 = " << L1 << ", L2 = " << L2 << ", L3 = " << L3 << endl;
}


Lattice::Lattice(int L1, int L2, int L3, bool isotropic, bool sianisotropy, bool magfield, bool dm)
{
    //this->L = L;
    this->L1 = L1;
    this->L2 = L2;
    this->L3 = L3;
    this->L = L1;
    this->isotropic = isotropic;
    this->sianisotropy = sianisotropy;
    this->magfield = magfield;
    this->dm = dm;

    notperiodic = false; // Default value. Changes if we choose a lattice with closed boundary conditions
    systemstrengthsgiven = false;
    extended = false;
    seed = 23;
    cout << "In L1, L2, L3 Lattice constructor. L1 = " << L1 << ", L2 = " << L2 << ", L3 = " << L3 << endl;
}

// Function for feeding in system parameters
void Lattice::setstrengths(vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    /* The in put arrays are ordered like this:
    hx, hy, hz, Dix, Diy, Diz; // sitestrengthsin
    J, Jy, Jz, Jxy, Jxz, Jyz;  // heisenbergin
    Dx, Dy, Dz;                // dm_in
    */
    cout << "In Lattice::setstrengths. Your chosen parameters will now be used" << endl;

    this->hx  = sitestrengthsin[0]; // I already have something similar to this...
    this->hy  = sitestrengthsin[1];
    this->hz  = sitestrengthsin[2];
    this->Dix = sitestrengthsin[3];
    this->Diy = sitestrengthsin[4];
    this->Diz = sitestrengthsin[5];

    this->J   = heisenbergin[0];
    this->Jx  = heisenbergin[1];
    this->Jy  = heisenbergin[2];
    this->Jz  = heisenbergin[3];
    this->Jxy = heisenbergin[4];
    this->Jxz = heisenbergin[5];
    this->Jyz = heisenbergin[6];

    this->Dx = dm_in[0];
    this->Dy = dm_in[1];
    this->Dz = dm_in[2];

    // Just in case we forget to call this one:
    systemstrengthsgiven = true;
}

/*
void Lattice::majordebug()
{
    MBUG = true;
}
*/

// Default function setting system strengths if the first one is forgotten
void Lattice::givestrengths_automatic()
{
    hx  = 1; // I already have something similar to this...
    hy  = 1;
    hz  = 1;
    Dix = 0;
    Diy = 1;
    Diz = 0;

    J   = 1;
    Jx  = 1;
    Jy  = 1;
    Jz  = 1;
    Jxy = 1;
    Jxz = 1;
    Jyz = 1;

    Dx = 1;
    Dy = 1;
    Dz = 1;

    cout << "WARNING!!! System parameters were not given. Parameters set to default value." << endl;
}

int Lattice::findneighbour2D(int n, int toi, int toj)
{   // Not tested, should work
    int i = n/L2;
    int j = n%L2;

    return ((L1+i+toi)%L1)*L2+((L2+j+toj)%L2);
}

int Lattice::findneighbour(int n, int toi, int toj, int tok)
{   // Gives correct output (at least for 2x2x2 fcc)
    int i = n/(L2*L3);
    int j = (int)n/L3 - (int)n/(L2*L3)*L2;
    int k = n%L3;

    return ((L1+i+toi)%L1)*(L2*L3)+((L2+j+toj)%L2)*L3+(L3+k+tok)%L3;
}


void Lattice::chain_open_initialize()
{
    notperiodic = true;
    bool randomspins = false;  // Adding the option to change the spins
    cout << "NB! Bool notperiodic changed to true. Now operating with open BCs!" << endl;
    dim = 1;
    N = L1;
    no_of_neighbours = 2; // For the most part

    dimlengths = vector<int>(1);
    dimlengths[0] = L1;
    a1 = vector<double>(1);
    a1[0] = 1;
    b1 = vector<double>(1);
    b1[0] = 2*M_PI;

    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    if(!systemstrengthsgiven)    givestrengths_automatic();

    for(int n=0; n<N; n++)
    {
        // Finding the neighbours to n
        int no_of_neighbours_site  = 0;
        int no_of_nneighbours_site = 0;
        int np1 = (n+1)%N;
        int nm1 = (n+N-1)%N;
        int np2 = (n+2)%N;
        int nm2 = (n+N-2)%N;

        std::vector<Bond> bonds;
        std::vector<Bond> nextnearestbonds;

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

        bool increasing;
        // Excluding periodic neighbours
        if(np1>n)
        {
            increasing = true;
            bonds.push_back(Bond(n, np1, J, increasing));
            no_of_neighbours_site++;
        }
        if(nm1<n)
        {
            increasing = false;
            bonds.push_back(Bond(n, nm1, J, increasing));
            no_of_neighbours_site++;
        }

        // Next nearest neighbours
        if(np2>n)
        {
            cout << "n = " << n << "; np2 = " << np2 << endl;
            increasing = true;
            nextnearestbonds.push_back(Bond(n, np2, Jx, increasing));
            no_of_nneighbours_site++;
        }
        if(nm2<n)
        {
            cout << "n = " << n << "; nm2 = " << nm2 << endl;
            increasing = false;
            nextnearestbonds.push_back(Bond(n, nm2, Jx, increasing));
            no_of_nneighbours_site++;
        }



        //cout << "Spin " << n << ", neighbours: " << np1 << ", " << nm1 << endl;

        // Is it too nested to make Site inherit Bond? ... Seems fair?
        // Send in bools
        sites.push_back(Site(n, no_of_neighbours_site, no_of_nneighbours_site, spinx, spiny, spinz, bonds, nextnearestbonds));
        // or
        //sites.push_back(Site(n, spinx, spiny, spinz, hx, hy, hz, Dix, Diy, Diz, bonds));
    }
}


void Lattice::chain_periodic_initialize()
{
    cout << "In chain_periodic_initialize" << endl;
    dim = 1;
    N = L1;
    no_of_neighbours = 2;

    dimlengths = vector<int>(1);
    dimlengths[0] = L1;
    a1 = vector<double>(1);
    a1[0] = 1;
    b1 = vector<double>(1);
    b1[0] = 2*M_PI;

    double spinx, spiny, spinz;
    bool zdird = false;
    if(zdird)
    {
        spinx = 0;
        spiny = 0;
        spinz = 1;
    }
    else
    {
        double a = 1/sqrt(3);
        spinx = a;
        spiny = a;
        spinz = a;
    }

    if(!systemstrengthsgiven)    givestrengths_automatic();

    bool randomspins = true;
    for(int n=0; n<N; n++)
    {

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

        // Finding the neighbours to n
        // This should only be done once. And that is exactly what we are doing.
        // Doing modulo operations, as suggested in Newman & Barkema
        // These neighbours are consistent with the sketch in Newman & Barkema
        int np1 = (n+1)%N;
        int nm1 = (n+N-1)%N;

        // Next nearest neighbours. Do not make much sense for small chains:
        int nnp = (n+2)%N;
        int nnm = (n+N-2)%N;

        //cout << "Spin " << n << ", neighbours: " << np1 << ", " << nm1 << endl;


        std::vector<Bond> bonds;

        /*
        bool increasingp;
        bool increasingm;
        if(np1>n)    increasingp = true;
        else         increasingp = false;

        if(nm1>n)    increasingm = true;
        else         increasingm = false;

        bonds.push_back(Bond(n, np1, isotropic, dm, increasingp, bondints));  // Do I really need to send in n?
        bonds.push_back(Bond(n, nm1, isotropic, dm, increasingm, bondints));
        */

        // Making a lot of bond classes to be added to bonds.
        bonds.push_back(Bond(n, np1, J, true));  // Do I really need to send in n?
        bonds.push_back(Bond(n, nm1, J, false));

        // I use Jx here because I don't bother feeding in a designated variable
        std::vector<Bond> nextnearesty;
        std::vector<Bond> nextnearestz;

        nextnearesty.push_back(Bond(n, nnm, Jx, false));
        nextnearesty.push_back(Bond(n, nnp, Jx, true));
        nextnearestz.push_back(Bond(n, nnm, Jx, false));
        nextnearestz.push_back(Bond(n, nnp, Jx, true));

        sites.push_back(Site(n, spinx, spiny, spinz, bonds, nextnearesty, nextnearestz));
    }
}

void Lattice::quadratic_helical_initialize()
{   // This one is primarily for testing.
    //N = L*(L+1); // Look this up!
    dim = 2;
    N = L1*L2;
    no_of_neighbours = 4;

    // No of particles in each direction
    dimlengths = vector<int>(dim);
    dimlengths[0] = L1;
    dimlengths[1] = L2;

    // Lattice vectors
    a1 = vector<double>(2);
    a2 = vector<double>(2);
    a1[0] = 1;    a1[1] = 0;
    a2[0] = 0;    a2[1] = 1;

    // Reciprocal lattice vectors
    b1 = vector<double>(2);
    b2 = vector<double>(2);
    b1[0] = 2*M_PI;    b1[1] = 0;
    b2[0] = 0;         b2[1] = 2*M_PI;


    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    if(!systemstrengthsgiven)    givestrengths_automatic();

    std::vector<double> position_n = std::vector<double>(2);
    std::vector<int> coord_n = std::vector<int>(2);
    std::vector<int> neighbours_n = std::vector<int>(4);

    // Could have these inside the loop and add randomness.

    bool randomspins = true;
    for(int n=0; n<N; n++)
    {

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
        // Finding the neighbours
        int np1, nm1, npL, nmL;
        np1 = findneighbour2D(n,0,1);
        nm1 = findneighbour2D(n,0,-1);
        npL = findneighbour2D(n,1,0);
        nmL = findneighbour2D(n,-1,0);

        neighbours_n[0] = np1;
        neighbours_n[1] = nm1;
        neighbours_n[2] = npL;
        neighbours_n[3] = nmL;

        siteneighbours.push_back(neighbours_n);


        std::vector<Bond> bonds;

        /*
        bool increasingp1, increasingm1, increasingpL, increasingmL;
        if(np1>n)    increasingp1 = true;
        else         increasingp1 = false;

        if(nm1>n)    increasingm1 = true;
        else         increasingm1 = false;

        if(npL>n)    increasingpL = true;
        else         increasingpL = false;

        if(nmL>n)    increasingmL = true;
        else         increasingmL = false;
        */

        // Making a lot of Bond classes to be added to vector of bonds.
        bonds.push_back(Bond(n, np1, J, true));  // Do I really need to send in n?
        bonds.push_back(Bond(n, nm1, J, false));
        bonds.push_back(Bond(n, npL, J, true));
        bonds.push_back(Bond(n, nmL, J, false));

        // Send in bools
        sites.push_back(Site(n, spinx, spiny, spinz, bonds));

        // Positions (row-major order)
        position_n[0] = 1.0*((int)n/L2); // n1
        position_n[1] = 1.0*((int)n%L2); // n2. Grid length a set to one, regulate by strength of interactions

        coord_n[0] = ((int)n/L2); // n1
        coord_n[1] = ((int)n%L2); // n2

        sitepositions.push_back(position_n);
        sitecoordinates.push_back(coord_n);

    } // End of loop over all sites ( = loop over n)
}

void Lattice::quadratic_helical_initialize_extended()
{
    dim = 2;
    N = L1*L2;
    no_of_neighbours = 4;
    extended = true;

    // No of particles in each direction
    dimlengths = vector<int>(dim);
    dimlengths[0] = L1;
    dimlengths[1] = L2;

    // Lattice vectors
    a1 = vector<double>(2);
    a2 = vector<double>(2);
    a1[0] = 1;    a1[1] = 0;
    a2[0] = 0;    a2[1] = 1;

    // Reciprocal lattice vectors
    b1 = vector<double>(2);
    b2 = vector<double>(2);
    b1[0] = 2*M_PI;    b1[1] = 0;
    b2[0] = 0;         b2[1] = 2*M_PI;

    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    if(!systemstrengthsgiven)    givestrengths_automatic();

    std::vector<double> position_n = std::vector<double>(2);
    std::vector<int> coord_n = std::vector<int>(2);
    std::vector<int> neighbours_n = std::vector<int>(4);

    // Extended version: Can have different interactions in different directions
    string x = "x";
    string y = "y";

    bool randomspins = true;
    for(int n=0; n<N; n++)
    {

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
        // Finding the neighbours
        int np1, nm1, npL, nmL;
        np1 = findneighbour2D(n,0,1);
        nm1 = findneighbour2D(n,0,-1);
        npL = findneighbour2D(n,1,0);
        nmL = findneighbour2D(n,-1,0);

        neighbours_n[0] = np1;
        neighbours_n[1] = nm1;
        neighbours_n[2] = npL;
        neighbours_n[3] = nmL;

        siteneighbours.push_back(neighbours_n);

        /*
        bool increasingp1, increasingm1, increasingpL, increasingmL;
        if(np1>n)    increasingp1 = true;
        else         increasingp1 = false;

        if(nm1>n)    increasingm1 = true;
        else         increasingm1 = false;

        if(npL>n)    increasingpL = true;
        else         increasingpL = false;

        if(nmL>n)    increasingmL = true;
        else         increasingmL = false;
        */


        std::vector<Bond> bonds;

        // Making a lot of Bond classes to be added to vector of bonds.
        bonds.push_back(Bond(n, np1, Jy, true, y)); // Not sure there is much point of this...
        bonds.push_back(Bond(n, nm1, Jy, false, y));
        bonds.push_back(Bond(n, npL, Jx, true, x));
        bonds.push_back(Bond(n, nmL, Jx, false, x));

        // Send in bools
        sites.push_back(Site(n, spinx, spiny, spinz, bonds));

        // Positions (row-major order)
        position_n[0] = 1.0*((int)n/L2); // n1
        position_n[1] = 1.0*((int)n%L2); // n2. Grid length a set to one, regulate by strength of interactions

        coord_n[0] = ((int)n/L2); // n1
        coord_n[1] = ((int)n%L2); // n2

        sitepositions.push_back(position_n);
        sitecoordinates.push_back(coord_n);

    } // End of loop over all sites ( = loop over n)
}

void Lattice::cubic_helical_initialize()
{
    //N = L*(L+1)*(L+1); // Look this up!
    dim = 3;
    N = L1*L2*L3;
    no_of_neighbours = 6;

    // No of particles in each direction
    dimlengths    = vector<int>(dim);    
    dimlengths[0] = L1;
    dimlengths[1] = L2;
    dimlengths[2] = L3;


    // Lattice vectors
    a1 = vector<double>(3); // Should they be double?
    a2 = vector<double>(3);
    a3 = vector<double>(3);
    a1[0] = 1;    a1[1] = 0;    a1[2] = 0;
    a2[0] = 0;    a2[1] = 1;    a2[2] = 0;
    a3[0] = 0;    a3[1] = 0;    a3[2] = 1;

    // Reciprocal lattice vectors
    b1 = vector<double>(3); // Should they be double?
    b2 = vector<double>(3);
    b3 = vector<double>(3);
    b1[0] = 2*M_PI;    b1[1] = 0;         b1[2] = 0;
    b2[0] = 0;         b2[1] = 2*M_PI;    b2[2] = 0;
    b3[0] = 0;         b3[1] = 0;         b3[2] = 2*M_PI;

    // Temporary fix. May want to send in L1, L2 and L3 to Lattice and deal with a 'rectangular' crystal
    //int L1 = L;
    //int L2 = L;

    //cout << "Do not use cubib_helical_initialize() yet. Not implemented" << endl;
    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    if(!systemstrengthsgiven)    givestrengths_automatic();

    std::vector<double> position_n = std::vector<double>(3);
    std::vector<int> coord_n = std::vector<int>(3);
    std::vector<int> neighbours_n = std::vector<int>(6);

    // Could have these inside the loop and add randomness.

    bool randomspins = true;
    for(int n=0; n<N; n++)
    {

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
        int np1, nm1, npL, nmL, npL2, nmL2;
        np1  = findneighbour(n,0,0,1);
        nm1  = findneighbour(n,0,0,-1);
        npL  = findneighbour(n,0,1,0);
        nmL  = findneighbour(n,0,-1,0);
        npL2 = findneighbour(n,1,0,0);
        nmL2 = findneighbour(n,-1,0,0);

        neighbours_n[0] = np1;
        neighbours_n[1] = nm1;
        neighbours_n[2] = npL;
        neighbours_n[3] = nmL;
        neighbours_n[4] = npL2;
        neighbours_n[5] = nmL2;

        siteneighbours.push_back(neighbours_n);

        std::vector<Bond> bonds;

        /*
        bool increasingp1, increasingm1, increasingpL, increasingmL, increasingpL2, increasingmL2;
        if(np1>n)     increasingp1 = true;
        else          increasingp1 = false;
        if(nm1>n)     increasingm1 = true;
        else          increasingm1 = false;
        if(npL>n)     increasingpL = true;
        else          increasingpL = false;
        if(nmL>n)     increasingmL = true;
        else          increasingmL = false;
        if(npL2>n)    increasingpL2 = true;
        else          increasingpL2 = false;
        if(nmL2>n)    increasingmL2 = true;
        else          increasingmL2 = false;
        */

        // Making a lot of bond classes to be added to bonds.
        bonds.push_back(Bond(n, np1, J, true));  // Do I really need to send in n?
        bonds.push_back(Bond(n, nm1, J, false));
        bonds.push_back(Bond(n, npL, J, true));
        bonds.push_back(Bond(n, nmL, J, false));
        bonds.push_back(Bond(n, npL2, J, true));
        bonds.push_back(Bond(n, nmL2, J, false));

        sites.push_back(Site(n, spinx, spiny, spinz, bonds));

        // Giving the position (row-major order)
        int n1, n2, n3;
        n1 = n/(L2*L3);
        n2 = n/L3 - n/(L2*L3)*L2;
        n3 = n%L3;

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

void Lattice::cubic_helical_initialize_extended()
{
    dim = 3;
    N = L1*L2*L3;
    no_of_neighbours = 6;
    extended = true;

    // No of particles in each direction
    dimlengths    = vector<int>(dim);
    dimlengths[0] = L1;
    dimlengths[1] = L2;
    dimlengths[2] = L3;

    // Lattice vectors
    a1 = vector<double>(3); // Should they be double?
    a2 = vector<double>(3);
    a3 = vector<double>(3);
    a1[0] = 1;    a1[1] = 0;    a1[2] = 0;
    a2[0] = 0;    a2[1] = 1;    a2[2] = 0;
    a3[0] = 0;    a3[1] = 0;    a3[2] = 1;

    // Reciprocal lattice vectors
    b1 = vector<double>(3); // Should they be double?
    b2 = vector<double>(3);
    b3 = vector<double>(3);
    b1[0] = 2*M_PI;    b1[1] = 0;         b1[2] = 0;
    b2[0] = 0;         b2[1] = 2*M_PI;    b2[2] = 0;
    b3[0] = 0;         b3[1] = 0;         b3[2] = 2*M_PI;

    // Temporary fix. May want to send in L1, L2 and L3 to Lattice and deal with a 'rectangular' crystal
    //int L1 = L;
    //int L2 = L;

    //cout << "Do not use cubib_helical_initialize() yet. Not implemented" << endl;
    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    if(!systemstrengthsgiven)    givestrengths_automatic();

    std::vector<double> position_n = std::vector<double>(3);
    std::vector<int> coord_n = std::vector<int>(3);
    std::vector<int> neighbours_n = std::vector<int>(6);

    // Extended version: Can have different interactions in different directions
    string x = "x";
    string y = "y";
    string z = "z";

    for(int n=0; n<N; n++)
    {
        int np1, nm1, npL, nmL, npL2, nmL2;
        np1  = findneighbour(n,0,0,1);
        nm1  = findneighbour(n,0,0,-1);
        npL  = findneighbour(n,0,1,0);
        nmL  = findneighbour(n,0,-1,0);
        npL2 = findneighbour(n,1,0,0);
        nmL2 = findneighbour(n,-1,0,0);

        neighbours_n[0] = np1;
        neighbours_n[1] = nm1;
        neighbours_n[2] = npL;
        neighbours_n[3] = nmL;
        neighbours_n[4] = npL2;
        neighbours_n[5] = nmL2;

        siteneighbours.push_back(neighbours_n);

        /*
        bool increasingp1, increasingm1, increasingpL, increasingmL, increasingpL2, increasingmL2;
        if(np1>n)     increasingp1 = true;
        else          increasingp1 = false;
        if(nm1>n)     increasingm1 = true;
        else          increasingm1 = false;
        if(npL>n)     increasingpL = true;
        else          increasingpL = false;
        if(nmL>n)     increasingmL = true;
        else          increasingmL = false;
        if(npL2>n)    increasingpL2 = true;
        else          increasingpL2 = false;
        if(nmL2>n)    increasingmL2 = true;
        else          increasingmL2 = false;
        */

        std::vector<Bond> bonds;

        // Making a lot of bond classes to be added to bonds.
        bonds.push_back(Bond(n, np1, Jz, true, z));
        bonds.push_back(Bond(n, nm1, Jz, true, z));
        bonds.push_back(Bond(n, npL, Jy, true, y));
        bonds.push_back(Bond(n, nmL, Jy, true, y));
        bonds.push_back(Bond(n, npL2, Jx, true, x));
        bonds.push_back(Bond(n, nmL2, Jx, true, x));

        sites.push_back(Site(n, spinx, spiny, spinz, bonds));

        // Giving the position (row-major order)
        int n1, n2, n3;
        n1 = n/(L2*L3);
        n2 = n/L3 - n/(L2*L3)*L2;
        n3 = n%L3;

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


void Lattice::fcc_helical_initialize_extended()
{
    cout << "In fcc_helical_initialize_extended" << endl;
    bool DEBUG = false;
    dim = 3;
    N = L1*L2*L3;
    cout << "N: " << N << endl;
    no_of_neighbours = 12;
    extended = true;

    // No of particles in each direction
    dimlengths    = vector<int>(dim);

    dimlengths[0] = L1;
    dimlengths[1] = L2;
    dimlengths[2] = L3;


    // Lattice vectors
    a1 = vector<double>(3); // Should they be double?
    a2 = vector<double>(3);
    a3 = vector<double>(3);
    a1[0] = 0.5;    a1[1] = 0.5;    a1[2] = 0;
    a2[0] = 0;      a2[1] = 0.5;    a2[2] = 0.5;
    a3[0] = 0.5;    a3[1] = 0;      a3[2] = 0.5;

    // Reciprocal lattice vectors
    b1 = vector<double>(3); // Should they be double?
    b2 = vector<double>(3);
    b3 = vector<double>(3);
    b1[0] =  2*M_PI;     b1[1] =  2*M_PI;     b1[2] = -2*M_PI;
    b2[0] = -2*M_PI;     b2[1] =  2*M_PI;     b2[2] =  2*M_PI;
    b3[0] =  2*M_PI;     b3[1] = -2*M_PI;     b3[2] =  2*M_PI;

    //cout << "In fcc_helical_initialize_extended!!!!!!!" << endl;

    // Setting up the sites
    // We set up the matrix by having all spins in the same direction. Or maybe draw at random?
    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    if(!systemstrengthsgiven)    givestrengths_automatic();

    if(DEBUG)
    {
        cout << "hx : " << hx << "; hy : " << hy << "; hz : " << hz << "; Dix : " << Dix << "; Diy : " << Diy << "; Diz : " << Diz << endl;
        cout << "J : " << J << "; Jy : " << Jy << "; Jz : " << Jz << "; Jxy : " << Jxy << "; Jxz : " << Jxz << "; Jyz : " << Jyz << endl;
        cout << " Dx : " << Dx << "; Dy : " << Dy << "; Dz : " << Dz << endl;
    }

    std::vector<double> position_n = std::vector<double>(3);
    std::vector<int> coord_n       = std::vector<int>(3);
    std::vector<int> neighbours_n  = std::vector<int>(12);

    string xy = "xy";
    string xz = "xz";
    string yz = "yz";

    // Could have these inside the loop and add randomness.
    bool randomspins = false;
    if(DEBUG)    cout << "Now entering the loop" << endl;
    for(int n=0; n<N; n++)
    {
        if(DEBUG)    cout << "In loop, n = " << n << endl;

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

        // Finding the neighbours to n
        int np1, nm1, npL, npL2, npLm1, npL2m1, npL2mL, nmL, nmL2, nmLm1, nmL2m1, nmL2mL;
        np1    = findneighbour(n,0,0,1);  // neighbour in positive a3-direction (row-major ordering). Jxz
        nm1    = findneighbour(n,0,0,-1); // neighbour in negative a3-direction (row-major ordering). Jxz
        npL    = findneighbour(n,0,1,0);  // neighbour in positive a2-direction (row-major ordering). Jyz
        nmL    = findneighbour(n,0,-1,0); // neighbour in negative a2-direction (row-major ordering). Jyz
        npL2   = findneighbour(n,1,0,0);  // neighbour in positive a1-direction (row-major ordering). Jxy
        nmL2   = findneighbour(n,-1,0,0); // neighbour in negative a1-direction (row-major ordering). Jxy
        npLm1  = findneighbour(n,0,1,-1); // Jxy
        nmLm1  = findneighbour(n,0,-1,1); // Jxy
        npL2m1 = findneighbour(n,1,0,-1); // Jyz
        nmL2m1 = findneighbour(n,-1,0,1); // Jyz
        npL2mL = findneighbour(n,1,-1,0); // Jxz
        nmL2mL = findneighbour(n,-1,1,0); // Jxz

        neighbours_n[0] = np1;
        neighbours_n[1] = nm1;
        neighbours_n[2] = npL;
        neighbours_n[3] = nmL;
        neighbours_n[4] = npL2;
        neighbours_n[5] = nmL2;
        neighbours_n[6] = npLm1;
        neighbours_n[7] = nmLm1;
        neighbours_n[8] = npL2m1;
        neighbours_n[9] = nmL2m1;
        neighbours_n[10] = npL2mL;
        neighbours_n[11] = nmL2mL;

        siteneighbours.push_back(neighbours_n);

        /*
        bool increasingp1, increasingm1, increasingpL, increasingmL, increasingpL2, increasingmL2;
        bool increasingpLm1, increasingmLm1, increasingpL2m1, increasingmL2m1, increasingpL2mL, increasingmL2mL;
        if(np1>n)     increasingp1 = true;
        else          increasingp1 = false;
        if(nm1>n)     increasingm1 = true;
        else          increasingm1 = false;
        if(npL>n)     increasingpL = true;
        else          increasingpL = false;
        if(nmL>n)     increasingmL = true;
        else          increasingmL = false;
        if(npL2>n)    increasingpL2 = true;
        else          increasingpL2 = false;
        if(nmL2>n)    increasingmL2 = true;
        else          increasingmL2 = false;
        if(npLm1>n)     increasingpLm1  = true;
        else            increasingpLm1  = false;
        if(nmLm1>n)     increasingmLm1  = true;
        else            increasingmLm1  = false;
        if(npL2m1>n)    increasingpL2m1 = true;
        else            increasingpL2m1 = false;
        if(nmL2m1>n)    increasingmL2m1 = true;
        else            increasingmL2m1 = false;
        if(npL2mL>n)    increasingpL2mL = true;
        else            increasingpL2mL = false;
        if(nmL2mL>n)    increasingmL2mL = true;
        else            increasingmL2mL = false;
        */


        // Should I have some bool that determines whether or not I need these?
        // This part probably need fixing!!
        // Finding the next nearest neighbours in the y- and z- direction
        // (The distance between sites in the x-direction is larger, so this dir can be neglected.)
        int nnyp, nnym, nnzp, nnzm;

        nnyp = findneighbour(n, 1, 1,-1);
        nnym = findneighbour(n,-1,-1, 1);
        nnzp = findneighbour(n,-1, 1, 1);
        nnzm = findneighbour(n, 1,-1,-1);

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

        // Should I have some bool that determines whether or not I need these?
        std::vector<Bond> nextnearesty;
        std::vector<Bond> nextnearestz;

        nextnearesty.push_back(Bond(n, nnym, Jy, false));
        nextnearesty.push_back(Bond(n, nnyp, Jy, true));
        nextnearestz.push_back(Bond(n, nnzm, Jz, false));
        nextnearestz.push_back(Bond(n, nnzp, Jz, true));


        if(DEBUG)    cout << "Setting the bonds" << endl;

        //Bond(siteindex1, siteindex2, J, increasing, direction, bondints);

        // Making a lot of bond classes to be added to bonds.
        // Should I send in J separately (to easier allow for difference in strength in different directions)
        bonds.push_back(Bond(n, np1, Jxz,  true, xz));  // Do I really need to send in n?
        if(DEBUG)    cout << "Bond 1 done" << endl;
        bonds.push_back(Bond(n, nm1, Jxz, false, xz));
        if(DEBUG)    cout << "Bond 2 done" << endl;
        bonds.push_back(Bond(n, npL, Jyz, true, yz));
        //cout << "Jyz bond set" << endl;
        if(DEBUG)    cout << "Bond 3 done" << endl;
        bonds.push_back(Bond(n, nmL, Jyz, false, yz));
        //cout << "Jyz bond set" << endl;
        if(DEBUG)    cout << "Bond 4 done" << endl;
        bonds.push_back(Bond(n, npL2, Jxy, true, xy));
        if(DEBUG)    cout << "Bond 5 done" << endl;
        bonds.push_back(Bond(n, nmL2, Jxy, false, xy));
        if(DEBUG)    cout << "Bond 6 done" << endl;
        bonds.push_back(Bond(n, npLm1, Jxy, true, xy));
        if(DEBUG)    cout << "Bond 7 done" << endl;
        bonds.push_back(Bond(n, nmLm1, Jxy, false, xy));
        if(DEBUG)    cout << "Bond 8 done" << endl;
        bonds.push_back(Bond(n, npL2m1, Jyz, true, yz));
        //cout << "Jyz bond set" << endl;
        if(DEBUG)    cout << "Bond 9 done" << endl;
        bonds.push_back(Bond(n, nmL2m1, Jyz, false, yz));
        //cout << "Jyz bond set" << endl;
        if(DEBUG)    cout << "Bond 10 done" << endl;
        bonds.push_back(Bond(n, npL2mL, Jxz, true, xz));
        if(DEBUG)    cout << "Bond 11 done" << endl;
        bonds.push_back(Bond(n, nmL2mL, Jxz, false, xz));
        if(DEBUG)    cout << "Bond 12 done" << endl;

        if(DEBUG)    cout << "Done setting the bonds. Setting the sites" << endl;
        sites.push_back(Site(n, spinx, spiny, spinz, bonds, nextnearesty, nextnearestz));

        if(DEBUG)    cout << "Giving the position of the site in the fcc" << endl;
        // Giving the position of the fcc (when saved in row-major order)
        // Could change it back, probably. Ordering not important, I guess. The same as rotating
        // Should probably just choose something and stick with it.
        int n1, n2, n3;

        n1 = (int)n/(L2*L3);
        n2 = (int)n/L3 - (int)n/(L3*L2)*L2;
        n3 = (int)n%L3;

        //cout << "Site " << n << "; [i,j,k] = [" << n1 << "," << n2 << "," << n3 << "]" << endl;

        double xpos = 0.5*(n1+n3);  // Could possibly include the grid length a
        double ypos = 0.5*(n1+n2);
        double zpos = 0.5*(n2+n3);

        //cout << "n : " << n <<  " Position: [" << xpos << "," << ypos << "," << zpos << "]" << endl;

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



void Lattice::fcc_helical_initialize_extended_yopen()
{
    cout << "In fcc_helical_initialize_extended" << endl;
    bool DEBUG = false;
    dim = 3;
    N = L1*L2*L3;
    cout << "N: " << N << endl;
    no_of_neighbours = 12;
    extended = true;

    // No of particles in each direction
    dimlengths    = vector<int>(dim);
    dimlengths[0] = L1;
    dimlengths[1] = L2;
    dimlengths[2] = L3;

    // Lattice vectors
    a1 = vector<double>(3);
    a2 = vector<double>(3);
    a3 = vector<double>(3);
    a1[0] = 0.5;    a1[1] = 0.5;    a1[2] = 0;
    a2[0] = 0;      a2[1] = 0.5;    a2[2] = 0.5;
    a3[0] = 0.5;    a3[1] = 0;      a3[2] = 0.5;

    // Reciprocal lattice vectors
    b1 = vector<double>(3);
    b2 = vector<double>(3);
    b3 = vector<double>(3);
    b1[0] =  2*M_PI;     b1[1] =  2*M_PI;     b1[2] = -2*M_PI;
    b2[0] = -2*M_PI;     b2[1] =  2*M_PI;     b2[2] =  2*M_PI;
    b3[0] =  2*M_PI;     b3[1] = -2*M_PI;     b3[2] =  2*M_PI;

    // Setting up the sites
    // We set up the matrix by having all spins in the same direction. Or maybe draw at random?
    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    if(!systemstrengthsgiven)    givestrengths_automatic();

    if(DEBUG)
    {
        cout << "hx : " << hx << "; hy : " << hy << "; hz : " << hz << "; Dix : " << Dix << "; Diy : " << Diy << "; Diz : " << Diz << endl;
        cout << "J : " << J << "; Jy : " << Jy << "; Jz : " << Jz << "; Jxy : " << Jxy << "; Jxz : " << Jxz << "; Jyz : " << Jyz << endl;
        cout << " Dx : " << Dx << "; Dy : " << Dy << "; Dz : " << Dz << endl;
    }

    std::vector<double> position_n = std::vector<double>(3);
    std::vector<int> coord_n       = std::vector<int>(3);
    std::vector<int> neighbours_n;

    string xy = "xy";
    string xz = "xz";
    string yz = "yz";

    // Giving the position of the fcc (when saved in row-major order)
    // Could change it back, probably. Ordering not important, I guess. The same as rotating
    // Should probably just choose something and stick with it.
    int n1, n2, n3;

    n1 = (int)n/(L2*L3);
    n2 = (int)n/L3 - (int)n/(L3*L2)*L2;
    n3 = (int)n%L3;

    //cout << "Site " << n << "; [i,j,k] = [" << n1 << "," << n2 << "," << n3 << "]" << endl;

    double xpos = 0.5*(n1+n3);  // Could possibly include the grid length a
    double ypos = 0.5*(n1+n2);
    double zpos = 0.5*(n2+n3);

    //cout << "n : " << n <<  " Position: [" << xpos << "," << ypos << "," << zpos << "]" << endl;

    position_n[0] = xpos;
    position_n[1] = ypos;
    position_n[2] = zpos;

    coord_n[0] = n1;     // Have these here in case I want to check
    coord_n[1] = n2;
    coord_n[2] = n3;

    sitepositions.push_back(position_n);
    sitecoordinates.push_back(coord_n);

    // Could have these inside the loop and add randomness.
    bool randomspins = false;
    if(DEBUG)    cout << "Now entering the loop" << endl;
    for(int n=0; n<N; n++)
    {
        if(DEBUG)    cout << "In loop, n = " << n << endl;

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

        // Finding the neighbours to n
        int no_of_neighbours_site = 12;
        if(ypos==(L2-1)) no_of_neighbours_site = 8; // If we are at the end
        else if(ypos==0) no_of_neighbours_site = 8; // If we are at the start

        int np1, nm1, npL, npL2, npLm1, npL2m1, npL2mL, nmL, nmL2, nmLm1, nmL2m1, nmL2mL;
        np1    = findneighbour(n,0,0,1);  // neighbour in positive a3-direction (row-major ordering). Jxz
        nm1    = findneighbour(n,0,0,-1); // neighbour in negative a3-direction (row-major ordering). Jxz
        npL    = findneighbour(n,0,1,0);  // neighbour in positive a2-direction (row-major ordering). Jyz
        nmL    = findneighbour(n,0,-1,0); // neighbour in negative a2-direction (row-major ordering). Jyz
        npL2   = findneighbour(n,1,0,0);  // neighbour in positive a1-direction (row-major ordering). Jxy
        nmL2   = findneighbour(n,-1,0,0); // neighbour in negative a1-direction (row-major ordering). Jxy
        npLm1  = findneighbour(n,0,1,-1); // Jxy
        nmLm1  = findneighbour(n,0,-1,1); // Jxy
        npL2m1 = findneighbour(n,1,0,-1); // Jyz
        nmL2m1 = findneighbour(n,-1,0,1); // Jyz
        npL2mL = findneighbour(n,1,-1,0); // Jxz
        nmL2mL = findneighbour(n,-1,1,0); // Jxz

        neighbours_n.push_back(np1);
        neighbours_n.push_back(nm1);
        // Bonds in the y-direction
        if(ypos==(L2-1))
        {   // At the end of the lattice in the y-dir
            // We have open BCs in that direction
            neighbours_n.push_back(nmL);
            neighbours_n.push_back(nmL2);
            neighbours_n.push_back(nmLm1);
            neighbours_n.push_back(nmL2m1);
        }
        else if(ypos==0)
        {   // At the beginning of the lattice in the y-dir
            // We have open BCs in that direction
            neighbours_n.push_back(npL);
            neighbours_n.push_back(npL2);
            neighbours_n.push_back(npLm1);
            neighbours_n.push_back(npL2m1);
        }
        else
        {   // In the middle of the lattice,
            // as y-coords are concerned.
            // Here we have all bonds.
            neighbours_n.push_back(npL);
            neighbours_n.push_back(nmL);
            neighbours_n.push_back(npL2);
            neighbours_n.push_back(nmL2);
            neighbours_n.push_back(npLm1);
            neighbours_n.push_back(nmLm1);
            neighbours_n.push_back(npL2m1);
            neighbours_n.push_back(nmL2m1);
        }
        neighbours_n.push_back(npL2mL);
        neighbours_n.push_back(nmL2mL);

        siteneighbours.push_back(neighbours_n);

        // Should I have some bool that determines whether or not I need these?
        // This part probably need fixing!!
        // Finding the next nearest neighbours in the y- and z- direction
        // (The distance between sites in the x-direction is larger, so this dir can be neglected.)
        int nnyp, nnym, nnzp, nnzm;

        nnyp = findneighbour(n, 1, 1,-1);
        nnym = findneighbour(n,-1,-1, 1);
        nnzp = findneighbour(n,-1, 1, 1);
        nnzm = findneighbour(n, 1,-1,-1);

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

        // Should I have some bool that determines whether or not I need these?
        std::vector<Bond> nextnearesty;
        std::vector<Bond> nextnearestz;

        if(ypos!=0)         nextnearesty.push_back(Bond(n, nnym, Jy, false));
        if(ypos!=(L2-1))    nextnearesty.push_back(Bond(n, nnyp, Jy, true));
        nextnearestz.push_back(Bond(n, nnzm, Jz, false));
        nextnearestz.push_back(Bond(n, nnzp, Jz, true));

        int no_of_nneighbours_site = 2;
        if(ypos==0 || ypos==(L2-1))    no_of_nneighbours_site = 1;

        if(DEBUG)    cout << "Setting the bonds" << endl;

        //Bond(siteindex1, siteindex2, J, increasing, direction, bondints);

        // Making a lot of bond classes to be added to bonds.
        // Should I send in J separately (to easier allow for difference in strength in different directions)
        bonds.push_back(Bond(n, np1, Jxz,  true, xz));  // Do I really need to send in n?
        if(DEBUG)    cout << "Bond 1 done" << endl;
        bonds.push_back(Bond(n, nm1, Jxz, false, xz));
        if(DEBUG)    cout << "Bond 2 done" << endl;
        if(ypos!=(L2-1)) bonds.push_back(Bond(n, npL, Jyz, true, yz));
        //cout << "Jyz bond set" << endl;
        if(DEBUG)    cout << "Bond 3 done" << endl;
        if(ypos!=0)     bonds.push_back(Bond(n, nmL, Jyz, false, yz));
        //cout << "Jyz bond set" << endl;
        if(DEBUG)    cout << "Bond 4 done" << endl;
        if(ypos!=(L2-1)) bonds.push_back(Bond(n, npL2, Jxy, true, xy));
        if(DEBUG)    cout << "Bond 5 done" << endl;
        if(ypos!=0)     bonds.push_back(Bond(n, nmL2, Jxy, false, xy));
        if(DEBUG)    cout << "Bond 6 done" << endl;
        if(ypos!=(L2-1)) bonds.push_back(Bond(n, npLm1, Jxy, true, xy));
        if(DEBUG)    cout << "Bond 7 done" << endl;
        if(ypos!=0)     bonds.push_back(Bond(n, nmLm1, Jxy, false, xy));
        if(DEBUG)    cout << "Bond 8 done" << endl;
        if(ypos!=(L2-1)) bonds.push_back(Bond(n, npL2m1, Jyz, true, yz));
        //cout << "Jyz bond set" << endl;
        if(DEBUG)    cout << "Bond 9 done" << endl;
        if(ypos!=0)     bonds.push_back(Bond(n, nmL2m1, Jyz, false, yz));
        //cout << "Jyz bond set" << endl;
        if(DEBUG)    cout << "Bond 10 done" << endl;
        bonds.push_back(Bond(n, npL2mL, Jxz, true, xz));
        if(DEBUG)    cout << "Bond 11 done" << endl;
        bonds.push_back(Bond(n, nmL2mL, Jxz, false, xz));
        if(DEBUG)    cout << "Bond 12 done" << endl;

        if(DEBUG)    cout << "Done setting the bonds. Setting the sites" << endl;
        sites.push_back(Site(n, no_of_neighbours_site, no_of_nneighbours_site, spinx, spiny, spinz, bonds, nextnearestbonds));


        if(DEBUG)    cout << "Giving the position of the site in the fcc" << endl;
    }
    cout << "Done with fcc_helical_initialize" << endl;
}

void Lattice::fcc_helical_initialize()
{
    bool DEBUG = false;
    // Should include something saying how the parameters are set.
    //N = L*(L+1)*(L+1); // Look this up!
    dim = 3;
    N = L1*L2*L3;
    no_of_neighbours = 12;

    // No of particles in each direction
    dimlengths    = vector<int>(dim);       
    dimlengths[0] = L1;
    dimlengths[1] = L2;
    dimlengths[2] = L3;


    // Lattice vectors
    a1 = vector<double>(3); // Should they be double?
    a2 = vector<double>(3);
    a3 = vector<double>(3);
    a1[0] = 0.5;    a1[1] = 0.5;    a1[2] = 0;
    a2[0] = 0;      a2[1] = 0.5;    a2[2] = 0.5;
    a3[0] = 0.5;    a3[1] = 0;      a3[2] = 0.5;

    // Reciprocal lattice vectors
    b1 = vector<double>(3); // Should they be double?
    b2 = vector<double>(3);
    b3 = vector<double>(3);
    b1[0] =  2*M_PI;     b1[1] =  2*M_PI;     b1[2] = -2*M_PI;
    b2[0] = -2*M_PI;     b2[1] =  2*M_PI;     b2[2] =  2*M_PI;
    b3[0] =  2*M_PI;     b3[1] = -2*M_PI;     b3[2] =  2*M_PI;

    cout << "In fcc_helical_initialize" << endl;

    // Setting up the sites
    // We set up the matrix by having all spins in the same direction. Or maybe draw at random?
    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    if(!systemstrengthsgiven)    givestrengths_automatic();

    std::vector<double> position_n = std::vector<double>(3);
    std::vector<int> coord_n = std::vector<int>(3);
    std::vector<int> neighbours_n = std::vector<int>(12);

    // Could have these inside the loop and add randomness.

    if(DEBUG)    cout << "Now entering the loop" << endl;
    bool randomspins = true;
    for(int n=0; n<N; n++)
    {

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
        if(DEBUG)    cout << "In loop, n = " << n << endl;
        // Finding the neighbours to n
        // NB!: So far, I have only added neighbours that are a distance 1 apart, 2 being the length of the cell.
        // So no linking two cell corner atoms together as of now.
        // This should only be done once. And that is exactly what we are doing.
        // Doing modulo operations, as suggested in Newman & Barkema
        // Verified this bit (against Newman & Barkema [i,j,k]).

        int np1, nm1, npL, npL2, npLm1, npL2m1, npL2mL, nmL, nmL2, nmLm1, nmL2m1, nmL2mL;
        np1    = findneighbour(n,0,0,1);  // neighbour in positive a3-direction (row-major ordering)
        nm1    = findneighbour(n,0,0,-1); // neighbour in negative a3-direction (row-major ordering)
        npL    = findneighbour(n,0,1,0);  // neighbour in positive a2-direction (row-major ordering)
        nmL    = findneighbour(n,0,-1,0); // neighbour in negative a2-direction (row-major ordering)
        npL2   = findneighbour(n,1,0,0);  // neighbour in positive a1-direction (row-major ordering)
        nmL2   = findneighbour(n,-1,0,0); // neighbour in negative a3-direction (row-major ordering)
        npLm1  = findneighbour(n,0,1,-1);
        nmLm1  = findneighbour(n,0,-1,1);
        npL2m1 = findneighbour(n,1,0,-1);
        nmL2m1 = findneighbour(n,-1,0,1);
        npL2mL = findneighbour(n,1,-1,0);
        nmL2mL = findneighbour(n,-1,1,0);

        neighbours_n[0] = np1;
        neighbours_n[1] = nm1;
        neighbours_n[2] = npL;
        neighbours_n[3] = nmL;
        neighbours_n[4] = npL2;
        neighbours_n[5] = nmL2;
        neighbours_n[6] = npLm1;
        neighbours_n[7] = nmLm1;
        neighbours_n[8] = npL2m1;
        neighbours_n[9] = nmL2m1;
        neighbours_n[10] = npL2mL;
        neighbours_n[11] = nmL2mL;

        siteneighbours.push_back(neighbours_n);

        /*
        bool increasingp1, increasingm1, increasingpL, increasingmL, increasingpL2, increasingmL2;
        bool increasingpLm1, increasingmLm1, increasingpL2m1, increasingmL2m1, increasingpL2mL, increasingmL2mL;

        if(np1>n)     increasingp1 = true;
        else          increasingp1 = false;
        if(nm1>n)     increasingm1 = true;
        else          increasingm1 = false;
        if(npL>n)     increasingpL = true;
        else          increasingpL = false;
        if(nmL>n)     increasingmL = true;
        else          increasingmL = false;
        if(npL2>n)    increasingpL2 = true;
        else          increasingpL2 = false;
        if(nmL2>n)    increasingmL2 = true;
        else          increasingmL2 = false;

        if(npLm1>n)     increasingpLm1  = true;
        else            increasingpLm1  = false;
        if(nmLm1>n)     increasingmLm1  = true;
        else            increasingmLm1  = false;
        if(npL2m1>n)    increasingpL2m1 = true;
        else            increasingpL2m1 = false;
        if(nmL2m1>n)    increasingmL2m1 = true;
        else            increasingmL2m1 = false;
        if(npL2mL>n)    increasingpL2mL = true;
        else            increasingpL2mL = false;
        if(nmL2mL>n)    increasingmL2mL = true;
        else            increasingmL2mL = false;
        */

        /*
        // Test this in some way...
        npL = (n+L3)%N;        // neighbour in positive a2-direction (row-major ordering)
        npL2 = (n+L2*L3)%N;     // neighbour in positive a1-direction (row-major ordering)
        npLm1 = (n+L3-1)%N;
        npL2m1 = (n+L2*L3-1)%N;
        npL2mL = (n+L2*L3-L3)%N;
        nmL = (n+N-L3)%N;      // neighbour in negative a2-direction (row-major ordering)
        nmL2 = (n+N-L2*L3)%N;   // neighbour in negative a1-direction (row-major ordering)
        nmLm1 = (n+N-L3+1)%N;
        nmL2m1 = (n+N-L2*L3+1)%N;
        nmL2mL = (n+N-L2*L3+L3)%N;
        */

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

        // In order to have a stronger interaction in one direction, I only need to send in true before
        // isotropic in the Bond initialization. Must find out which direction, though...
        if(DEBUG)    cout << "Setting the bonds" << endl;
        // Making a lot of bond classes to be added to bonds.
        bonds.push_back(Bond(n, np1, J, true));  // Do I really need to send in n?
        if(DEBUG)    cout << "Bond 1 done" << endl;
        bonds.push_back(Bond(n, nm1, J, false));
        if(DEBUG)    cout << "Bond 2 done" << endl;
        bonds.push_back(Bond(n, npL, J, true));
        if(DEBUG)    cout << "Bond 3 done" << endl;
        bonds.push_back(Bond(n, nmL, J, false));
        if(DEBUG)    cout << "Bond 4 done" << endl;
        bonds.push_back(Bond(n, npL2, J, true));
        if(DEBUG)    cout << "Bond 5 done" << endl;
        bonds.push_back(Bond(n, nmL2, J, false));
        if(DEBUG)    cout << "Bond 6 done" << endl;
        bonds.push_back(Bond(n, npLm1, J, true));
        if(DEBUG)    cout << "Bond 7 done" << endl;
        bonds.push_back(Bond(n, nmLm1, J, false));
        if(DEBUG)    cout << "Bond 8 done" << endl;
        bonds.push_back(Bond(n, npL2m1, J, true));
        if(DEBUG)    cout << "Bond 9 done" << endl;
        bonds.push_back(Bond(n, nmL2m1, J, false));
        if(DEBUG)    cout << "Bond 10 done" << endl;
        bonds.push_back(Bond(n, npL2mL, J, true));
        if(DEBUG)    cout << "Bond 11 done" << endl;
        bonds.push_back(Bond(n, nmL2mL, J, false));
        if(DEBUG)    cout << "Bond 12 done" << endl;

        cout << "Done setting the bonds. Setting the sites" << endl;
        sites.push_back(Site(n, spinx, spiny, spinz, bonds));

        // Listing the positions of the spins:

        //  /*
        if(DEBUG)    cout << "Giving the position of the site in the fcc" << endl;
        // Giving the position of the fcc (when saved in row-major order)
        // Could change it back, probably. Ordering not important, I guess. The same as rotating
        // However, this way, it is easy to map to the fftw-output.
        // Should probably just choose something and stick with it.
        int n1, n2, n3;

        n1 = n/(L2*L3);
        n2 = (int)n/L3 - (int)n/(L2*L3)*L2;
        n3 = n%L3;

        double xpos = 0.5*(n1+n3);  // Could possibly include the grid length a
        double ypos = 0.5*(n1+n2);
        double zpos = 0.5*(n2+n3);

        //cout << "[" << xpos << "," << ypos << "," << zpos << "]" << endl;

        position_n[0] = xpos;
        position_n[1] = ypos;
        position_n[2] = zpos;

        coord_n[0] = n1;     // Have these here in case I want to check
        coord_n[1] = n2;
        coord_n[2] = n3;

        sitepositions.push_back(position_n);
        sitecoordinates.push_back(coord_n);
        // */ /*

    }
    cout << "Done with fcc_helical_initialize_extended" << endl;
}

std::vector<double> Lattice::giveposition_fcc_lines(int i, int j, int k, char letter)
{
    // Our permutations
    //cout << "in giveposition_fcc_lines" << endl;

    if(letter=='x')         if(j>0)    j-=L2;
    if(letter=='y')         if(k>0)    k-=L3;
    if(letter=='z')         if(i>0)    i-=L1;

    //cout << "Gonna give'm positions" << endl;

    vector<double> position(3);
    position[0] = 0.5*(i+k);
    position[1] = 0.5*(i+j);
    position[2] = 0.5*(j+k);

    //cout << "Have set the positions" << endl;

    // Do something similar for the reciprocal lattice?
    return position;
}

std::vector<double> Lattice::giveqvector_fcc_lines(int i, int j, int k, char letter)
{
    // This function is sort of caput.
    // So we don't need to call fcc_initalize(_extended)
    // Reciprocal lattice vectors
    b1 = vector<double>(3); // Should they be double?
    b2 = vector<double>(3);
    b3 = vector<double>(3);
    b1[0] =  2*M_PI;     b1[1] =  2*M_PI;     b1[2] = -2*M_PI;
    b2[0] = -2*M_PI;     b2[1] =  2*M_PI;     b2[2] =  2*M_PI;
    b3[0] =  2*M_PI;     b3[1] = -2*M_PI;     b3[2] =  2*M_PI;

    // Our permutations
    // Drop this? For now at least?
    if(letter=='x')         if(j>0)    j-=L2;
    if(letter=='y')         if(k>0)    k-=L3;
    if(letter=='z')         if(i>0)    i-=L1;

    vector<double> qvector(3);
    qvector[0] = i*b1[0]/L1 + k*b3[0]/L3;
    qvector[1] = i*b1[1]/L1 + j*b2[1]/L2;
    qvector[2] = j*b2[2]/L2 + k*b3[2]/L3;

    cout << "WARNING! WARNING!:" << endl;
    cout << "This part of the program is outdates and incorrect. Call giveqvector_fcc() instead of giveqvector_fcc_lines()" << endl;

    return qvector;
}

std::vector<double> Lattice::giveqvector_fcc(int i, int j, int k)
{
    // So we don't need to call fcc_initalize(_extended)
    // Reciprocal lattice vectors
    b1 = vector<double>(3); // Should they be double?
    b2 = vector<double>(3);
    b3 = vector<double>(3);
    b1[0] =  2*M_PI;     b1[1] =  2*M_PI;     b1[2] = -2*M_PI;
    b2[0] = -2*M_PI;     b2[1] =  2*M_PI;     b2[2] =  2*M_PI;
    b3[0] =  2*M_PI;     b3[1] = -2*M_PI;     b3[2] =  2*M_PI;

    vector<double> qvector(3);
    qvector[0] =  i*b1[0]/L1 + j*b2[0]/L2 + k*b3[0]/L3;
    qvector[1] =  i*b1[1]/L1 + j*b2[1]/L2 + k*b3[1]/L3;
    qvector[2] =  i*b1[2]/L1 + j*b2[2]/L2 + k*b3[2]/L3;

    return qvector;
}

std::vector<int> Lattice::fccqyline()
{
    vector<int> qyline;
    vector<int> index;
    int i, j, k;
    double qx, qy, qz;
    vector<double> qvec = vector<double>(3);
    for(int n=0; n<N; n++)
    {
        index = sitecoordinates[n];
        i = index[0]; j = index[1]; k = index[2];

        qvec = giveqvector_fcc(i,j,k);
        qx =  qvec[0];        qy =  qvec[1];        qz =  qvec[2];

        if(abs(qx)<1e-18 && abs(qz)<1e-18)    qyline.push_back(n); // Increase it a little? Test if <1e-16?

    }
    return qyline;
}

std::vector<int> Lattice::fccqxline()
{
    vector<int> qxline;
    vector<int> index;
    int i, j, k;
    double qx, qy, qz;
    vector<double> qvec = vector<double>(3);
    for(int n=0; n<N; n++)
    {
        index = sitecoordinates[n];
        i = index[0]; j = index[1]; k = index[2];

        qvec = giveqvector_fcc(i,j,k);
        qx =  qvec[0];        qy =  qvec[1];        qz =  qvec[2];

        if(abs(qz)<1e-18 && abs(qy)<1e-18)    qxline.push_back(n); // Increase it a little? Test if <1e-16?

    }
    return qxline;
}

std::vector<int> Lattice::fccqzline()
{
    vector<int> qyline;
    vector<int> index;
    int i, j, k;
    double qx, qy, qz;
    vector<double> qvec = vector<double>(3);
    for(int n=0; n<N; n++)
    {
        index = sitecoordinates[n];
        i = index[0]; j = index[1]; k = index[2];

        qvec = giveqvector_fcc(i,j,k);
        qx =  qvec[0];        qy =  qvec[1];        qz =  qvec[2];

        if(abs(qx)<1e-18 && abs(qy)<1e-18)    qyline.push_back(n); // Increase it a little? Test if <1e-16?

    }
    return qyline;
}

std::vector<int> Lattice::fccqdline()
{
    vector<int> qdline;
    vector<int> index;
    int i, j, k;
    double qx, qy, qz;
    vector<double> qvec = vector<double>(3);
    for(int n=0; n<N; n++)
    {
        index = sitecoordinates[n];
        i = index[0]; j = index[1]; k = index[2];

        qvec = giveqvector_fcc(i,j,k);
        qx =  qvec[0];        qy =  qvec[1];        qz =  qvec[2];

        if(abs(qx-qy)<1e-16 && abs(qy-qz)<1e-16)    qdline.push_back(n); // Increase it a little? Test if <1e-16?

    }
    return qdline;
}


std::vector<int> Lattice::fccyline()
{   // Spatial y-line for fcc
    // Should this be written to file instead?
    int n = 0;
    int nprev = -1; // To start the loop
    vector<int> yline;
    while(nprev<n)  // This should fix it, I think...
    {
        nprev = n;
        yline.push_back(n);
        n = findneighbour(n, 1, 1, -1);
        cout << n << " ";
    }
    cout << endl;
    // I hope this doesn't run forever.
    return yline;
}

std::vector<int> Lattice::fccxline()
{   // Spatial x-line for fcc
    // Should this be written to file instead?
    int n = 0;
    int nprev = -1; // To start the loop
    vector<int> xline;
    while(nprev<n)  // This should fix it, I think...
    {
        nprev = n;
        xline.push_back(n);
        n = findneighbour(n, 1, -1, 1);
        cout << n << " ";
    }
    cout << endl;
    // I hope this doesn't run forever.
    return xline;
}

std::vector<int> Lattice::fcczline()
{   // Spatial z-line for fcc
    // Should this be written to file instead?
    int n = 0;
    int nprev = -1; // To start the loop
    vector<int> zline;

    /*zline.push_back(n);
    // Something odd
    while(nprev<n)  // This should fix it, I think...
    {
        nprev = n;
        n = findneighbour(n, -1, 1, 1);
        zline.push_back(n);
        cout << n << " ";
    }
    */

    int counter = 0;
    while(counter<L)  // Kind of a hack-y solution
    {
        zline.push_back(n);
        n = findneighbour(n, -1, 1, 1);
        cout << n << " ";
        counter++;
    }
    cout << endl;
    // I hope this doesn't run forever.
    return zline;
}

std::vector<int> Lattice::fccyline_shifted(double xshift, double zshift)
{
    cout << "In fccyline_shifted" << endl;
    int n1, n2, n3;
    double xc, yc, zc;

    vector<int> yline;
    for(int n=0; n<N; n++)
    {
        n1 = n/(L2*L3);
        n2 = n/L3 - n/(L2*L3)*L2;
        n3 = n%L3;

        xc = 0.5*(n1+n3);
        yc = 0.5*(n1+n2);
        zc = 0.5*(n2+n3);
        cout << "x: " << xc << "; y: " << yc << "; z:" << zc << endl;
        if(xc==xshift && zc==zshift)        yline.push_back(n);
    }
    return yline;
}

std::vector<int> Lattice::cubicyline()
{
    cout << "In cubicyline" << endl;
    int n1, n2, n3;

    vector<int> yline;
    for(int n=0; n<N; n++)
    {
        n1 = n/(L2*L3);
        n2 = n/L3 - n/(L2*L3)*L2;
        n3 = n%L3;
        if(n1==0 && n3==0)        yline.push_back(n);
    }
    return yline;
}

std::vector<int> Lattice::cubicxline()
{
    cout << "In cubicxline" << endl;
    int n1, n2, n3;

    vector<int> xline;
    for(int n=0; n<N; n++)
    {
        n1 = n/(L2*L3);
        n2 = n/L3 - n/(L2*L3)*L2;
        n3 = n%L3;
        cout << "n1 = " << n1 << "; n2 = " << n2 << "; n3 = " << n3 << endl;
        if(n2==0 && n3==0)            xline.push_back(n);
    }
    return xline;
}

std::vector<int> Lattice::cubiczline()
{
    cout << "In cubiczline" << endl;
    int n1, n2, n3;

    vector<int> zline;
    for(int n=0; n<N; n++)
    {
        n1 = n/(L2*L3);
        n2 = n/L3 - n/(L2*L3)*L2;
        n3 = n%L3;
        if(n1==0 && n2==0)        zline.push_back(n);
    }
    return zline;
}

std::vector<int> Lattice::quadrxline()
{
    cout << "In quadrxline" << endl;
    int n1, n2;

    vector<int> xline;
    for(int n=0; n<N; n++)
    {
        n1 = ((int)n/L2);
        n2 = ((int)n%L2);
        if(n2==0)        xline.push_back(n);
    }
    return xline;
}

std::vector<int> Lattice::quadryline()
{
    cout << "In quadryline" << endl;
    int n1, n2;

    vector<int> yline;
    for(int n=0; n<N; n++)
    {
        n1 = ((int)n/L2);
        n2 = ((int)n%L2);
        if(n1==0)        yline.push_back(n);
    }
    return yline;
}

std::vector<int> Lattice::diagline_fcc()
{
    cout << "In diagline_fcc" << endl;
    int n1, n2, n3, xint, yint,zint;

    vector<int> dline;
    for(int n=0; n<N; n++)
    {
        n1 = n/(L2*L3);
        n2 = n/L3 - n/(L2*L3)*L2;
        n3 = n%L3;

        xint = n1+n3;
        yint = n1+n2;
        zint = n2+n3;
        if(xint==yint && yint==zint)        dline.push_back(n);
    }
    return dline;
}

std::vector<int> Lattice::diagline_cubic()
{
    cout << "In diagline_cubic" << endl;
    int n1, n2, n3;

    vector<int> dline;
    for(int n=0; n<N; n++)
    {
        n1 = n/(L2*L3);
        n2 = n/L3 - n/(L2*L3)*L2;
        n3 = n%L3;
        if(n1==n2 && n2==n3)        dline.push_back(n);
    }
    return dline;
}

std::vector<int> Lattice::diagline_quadr()
{
    cout << "In diagline_quadr" << endl;
    int n1, n2;

    vector<int> dline;
    for(int n=0; n<N; n++)
    {
        n1 = ((int)n/L2); // n1
        n2 = ((int)n%L2);
        if(n1==n2)        dline.push_back(n);
    }
    return dline;
}
