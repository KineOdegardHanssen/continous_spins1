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
    dimequal    = true;  // Default for this constructor.
    systemstrengthsgiven = false;
}


Lattice::Lattice(int L1, int L2, int L3, bool isotropic, bool sianisotropy, bool magfield, bool dm)
{
    //this->L = L;
    this->L1 = L1;
    this->L2 = L2;
    this->L3 = L3;
    this->isotropic = isotropic;
    this->sianisotropy = sianisotropy;
    this->magfield = magfield;
    this->dm = dm;

    notperiodic = false; // Default value. Changes if we choose a lattice with closed boundary conditions
    dimequal    = false;  // Default for this constructor.
    systemstrengthsgiven = false;
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
    this->Jy  = heisenbergin[1];
    this->Jz  = heisenbergin[2];
    this->Jxy = heisenbergin[3];
    this->Jxz = heisenbergin[4];
    this->Jyz = heisenbergin[5];

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

int Lattice::findneighbour(int n, int toi, int toj, int tok)
{
    int i = n/(L2*L3);
    int j = (int)n/L3 - (int)n/(L2*L3)*L2;
    int k = n%L3;

    return ((L1+i+toi)%L1)*(L2*L3)+((L2+j+toj)%L2)*L2+(L3+k+tok)%L3;
}


void Lattice::chain_open_initialize()
{
    if(!dimequal)    L = L1;   // In case I am really stupid
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

    if(!systemstrengthsgiven)    givestrengths_automatic();

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
    cout << "In chain_periodic_initialize" << endl;
    if(!dimequal)    L = L1;   // In case I am really stupid
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

    if(!systemstrengthsgiven)    givestrengths_automatic();

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
    if(dimequal)    N = L*L;
    else            N = L1*L2;
    no_of_neighbours = 4;

    // No of particles in each direction
    dimlengths = vector<int>(dim);
    if(dimequal)
    {
        dimlengths[0] = L;
        dimlengths[1] = L;
    }
    else
    {
        dimlengths[0] = L1;
        dimlengths[1] = L2;
    }

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

    if(!systemstrengthsgiven)    givestrengths_automatic();

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
        int nm1 = (n+N-1)%N;

        int npL, nmL;
        if(dimequal)
        {
            npL = (n+L)%N;
            nmL = (n+N-L)%N;
        }
        else
        {
            npL = (n+L2)%N;
            nmL = (n+N-L2)%N;
        }

        std::vector<Bond> bonds;

        // Making a lot of bond classes to be added to bonds.
        bonds.push_back(Bond(n, np1, isotropic, dm, bondints));  // Do I really need to send in n?
        bonds.push_back(Bond(n, nm1, isotropic, dm, bondints));
        bonds.push_back(Bond(n, npL, isotropic, dm, bondints));
        bonds.push_back(Bond(n, nmL, isotropic, dm, bondints));


        // Send in bools
        sites.push_back(Site(n, sianisotropy, magfield, spinx, spiny, spinz, siteint, bonds));

        // Positions (row-major order)
        if(dimequal)
        {
            position_n[0] = 1.0*((int)n/L); // n1
            position_n[1] = 1.0*((int)n%L); // n2. Should allow for grid length a.

            coord_n[0] = ((int)n/L); // n1
            coord_n[1] = ((int)n%L); // n2
        }
        else
        {
            position_n[0] = 1.0*((int)n/L2); // n1
            position_n[1] = 1.0*((int)n%L2); // n2. Should allow for grid length a.

            coord_n[0] = ((int)n/L2); // n1
            coord_n[1] = ((int)n%L2); // n2
        }

        sitepositions.push_back(position_n);
        sitecoordinates.push_back(coord_n);

    } // End of loop over all sites ( = loop over n)
}

void Lattice::cubic_helical_initialize()
{
    //N = L*(L+1)*(L+1); // Look this up!
    dim = 3;
    if(dimequal)    N = L*L*L;
    else            N = L1*L2*L3;
    no_of_neighbours = 6;

    // No of particles in each direction
    dimlengths    = vector<int>(dim);
    if(dimequal)
    {
        dimlengths[0] = L;
        dimlengths[1] = L;
        dimlengths[2] = L;
    }
    else
    {
        dimlengths[0] = L1;
        dimlengths[1] = L2;
        dimlengths[2] = L3;
    }


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
    //int L1 = L;
    //int L2 = L;

    //cout << "Do not use cubib_helical_initialize() yet. Not implemented" << endl;
    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    if(!systemstrengthsgiven)    givestrengths_automatic();

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
        int nm1 = (n+N-1)%N;

        int npL, npL2, nmL, nmL2;
        if(dimequal)
        {
            npL = (n+L)%N;
            npL2 = (n+L*L)%N;
            nmL = (n+N-L)%N;
            nmL2 = (n+N-L*L)%N;
        }
        else
        {   // Must test these, I guess... Seems correct, however.
            npL = (n+L3)%N;
            npL2 = (n+L2*L3)%N;
            nmL = (n+N-L3)%N;
            nmL2 = (n+N-L2*L3)%N;
        }

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

        // Giving the position (row-major order)
        int n1, n2, n3;
        if(dimequal)    // Could also set L1 = L, L2 = L, L3 = L at the beginning of this function...
        {
            n1 = n/(L*L);
            n2 = n/L - n/(L*L)*L;
            n3 = n%L;
        }
        else
        {
            n1 = n/(L1*L2);
            n2 = n/L1 - n/(L1*L2)*L2;
            n3 = n%L1;
        }


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
    //cout << "In fcc_helical_initialize_extended" << endl;
    bool DEBUG = false;
    // Should include something saying how the parameters are set.
    //N = L*(L+1)*(L+1); // Look this up!
    dim = 3;
    if(dimequal)
    {
        L1 = L;
        L2 = L;
        L3 = L;
    }
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

    // Move these when neccessary
    std::vector<double> siteint = givethesiteints(Dix, Diy, Diz, hx, hy, hz, sianisotropy, magfield);
    // This could be simpler, actually:
    std::vector<double> bondints = std::vector<double>(3);
    bondints[0] = Dx;
    bondints[1] = Dy;
    bondints[2] = Dz;

    std::vector<double> position_n = std::vector<double>(3);
    std::vector<int> coord_n       = std::vector<int>(3);

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
        // NB!: So far, I have only added neighbours that are a distance 1 apart, 2 being the length of the cell.
        // So no linking two cell corner atoms together as of now.
        // This should only be done once. And that is exactly what we are doing.
        // Doing modulo operations, as suggested in Newman & Barkema
        // Verified this bit (against Newman & Barkema [i,j,k]).
        int np1 = (n+1)%N;        // neighbour in positive a3-direction (row-major ordering)        
        int nm1 = (n+N-1)%N;     // neighbour in negative a3-direction (row-major ordering)

        int npL, npL2, npLm1, npL2m1, npL2mL, nmL, nmL2, nmLm1, nmL2m1, nmL2mL;
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

        /*
        if(dimequal)
        {
            npL = (n+L)%N;        // neighbour in positive a2-direction (row-major ordering)
            npL2 = (n+L*L)%N;     // neighbour in positive a1-direction (row-major ordering)
            npLm1 = (n+L-1)%N;
            npL2m1 = (n+L*L-1)%N;
            npL2mL = (n+L*L-L)%N;
            nmL = (n+N-L)%N;      // neighbour in negative a2-direction (row-major ordering)
            nmL2 = (n+N-L*L)%N;   // neighbour in negative a1-direction (row-major ordering)
            nmLm1 = (n+N-L+1)%N;
            nmL2m1 = (n+N-L*L+1)%N;
            nmL2mL = (n+N-L*L+L)%N;
        }
        else
        {   // Test this in some way...
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
        }
        */

        // Should I have some bool that determines whether or not I need these?
        int nnyp, nnym, nnzp, nnzm;
        int deltanny, deltannz;
        if(dimequal)
        {
            deltanny = L*L+L-1;
            deltannz = -L*L+L-1;
        }
        else
        {
            deltanny = L2*L3+L3-1;
            deltannz = -L2*L3+L3-1; // This is negative...
        }
        nnyp = (n+deltanny)%N;
        nnym = (n+N-deltanny)%N;  // +N added to ensure it is not negative.
        nnzp = (n+N+deltannz)%N;  //            --"--
        nnzm = (n-deltannz)%N;


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

        nextnearesty.push_back(Bond(n, nnym, Jy));
        nextnearesty.push_back(Bond(n, nnyp, Jy));
        nextnearestz.push_back(Bond(n, nnzm, Jz));
        nextnearestz.push_back(Bond(n, nnzp, Jz));


        if(DEBUG)    cout << "Setting the bonds" << endl;
        // Making a lot of bond classes to be added to bonds.
        // Should I send in J separately (to easier allow for difference in strength in different directions)
        bonds.push_back(Bond(n, np1, Jxz, bondints));  // Do I really need to send in n?
        if(DEBUG)    cout << "Bond 1 done" << endl;
        bonds.push_back(Bond(n, nm1, Jxz, bondints));
        if(DEBUG)    cout << "Bond 2 done" << endl;
        bonds.push_back(Bond(n, npL, Jyz, bondints));
        if(DEBUG)    cout << "Bond 3 done" << endl;
        bonds.push_back(Bond(n, nmL, Jyz, bondints));
        if(DEBUG)    cout << "Bond 4 done" << endl;
        bonds.push_back(Bond(n, npL2, Jxy, bondints));
        if(DEBUG)    cout << "Bond 5 done" << endl;
        bonds.push_back(Bond(n, nmL2, Jxy, bondints));
        if(DEBUG)    cout << "Bond 6 done" << endl;
        bonds.push_back(Bond(n, npLm1, Jxy, bondints));
        if(DEBUG)    cout << "Bond 7 done" << endl;
        bonds.push_back(Bond(n, nmLm1, Jxy, bondints));
        if(DEBUG)    cout << "Bond 8 done" << endl;
        bonds.push_back(Bond(n, npL2m1, Jyz, bondints));
        if(DEBUG)    cout << "Bond 9 done" << endl;
        bonds.push_back(Bond(n, nmL2m1, Jyz, bondints));
        if(DEBUG)    cout << "Bond 10 done" << endl;
        bonds.push_back(Bond(n, npL2mL, Jxz, bondints));
        if(DEBUG)    cout << "Bond 11 done" << endl;
        bonds.push_back(Bond(n, nmL2mL, Jxz, bondints));
        if(DEBUG)    cout << "Bond 12 done" << endl;

        if(DEBUG)    cout << "Done setting the bonds. Setting the sites" << endl;
        sites.push_back(Site(n, sianisotropy, magfield, spinx, spiny, spinz, siteint, bonds, nextnearesty, nextnearestz));

        if(DEBUG)    cout << "Giving the position of the site in the fcc" << endl;
        // Giving the position of the fcc (when saved in row-major order)
        // Could change it back, probably. Ordering not important, I guess. The same as rotating
        // Should probably just choose something and stick with it.
        int n1, n2, n3;

        n1 = n/(L1*L2);
        n2 = (int)n/L3 - (int)n/(L3*L2)*L2;
        n3 = n%L3;

        double xpos = 0.5*(n1+n3);  // Could possibly include the grid length a
        double ypos = 0.5*(n1+n2);
        double zpos = 0.5*(n2+n3);

       // cout << "[" << xpos << "," << ypos << "," << zpos << "]" << endl;

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


void Lattice::fcc_helical_initialize()
{
    bool DEBUG = false;
    // Should include something saying how the parameters are set.
    //N = L*(L+1)*(L+1); // Look this up!
    dim = 3;
    if(dimequal)
    {
        L1 = L;
        L2 = L;
        L3 = L;
    }
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

    // Setting up the sites
    // We set up the matrix by having all spins in the same direction. Or maybe draw at random?
    double a = 1/sqrt(3);
    double spinx = a;
    double spiny = a;
    double spinz = a;

    if(!systemstrengthsgiven)    givestrengths_automatic();

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
        // Verified this bit (against Newman & Barkema [i,j,k]).
        int np1 = (n+1)%N;
        int nm1 = (n+N-1)%N;


        int npL, npL2, npLm1, npL2m1, npL2mL, nmL, nmL2, nmLm1, nmL2m1, nmL2mL;
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

        /*
        if(dimequal)
        {
            npL = (n+L)%N;        // neighbour in positive a2-direction (row-major ordering)
            npL2 = (n+L*L)%N;     // neighbour in positive a1-direction (row-major ordering)
            npLm1 = (n+L-1)%N;
            npL2m1 = (n+L*L-1)%N;
            npL2mL = (n+L*L-L)%N;
            nmL = (n+N-L)%N;      // neighbour in negative a2-direction (row-major ordering)
            nmL2 = (n+N-L*L)%N;   // neighbour in negative a1-direction (row-major ordering)
            nmLm1 = (n+N-L+1)%N;
            nmL2m1 = (n+N-L*L+1)%N;
            nmL2mL = (n+N-L*L+L)%N;
        }
        else
        {   // Test this in some way...
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
        }
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

        cout << "Done setting the bonds. Setting the sites" << endl;
        sites.push_back(Site(n, sianisotropy, magfield, spinx, spiny, spinz, siteint, bonds));

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
