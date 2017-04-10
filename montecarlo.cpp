#include "montecarlo.h"

MonteCarlo::MonteCarlo()
{
}

MonteCarlo::MonteCarlo(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{   // Should probably just send this in to the other initializer as MonteCarlo(int L, int L, int L, ...)
    randomtest = false;
    // Handling runningints (wouldn't want to vary this in one class instance, I guess.)
    this->eqsteps = eqsteps;
    this->mcsteps_inbin = mcsteps_inbin;
    this->no_of_bins = no_of_bins;

    // Setting the bools
    this->isotropic = isotropic;
    this->sianisotropy = sianisotropy;
    this->magfield = magfield;
    this->dm = dm;

    if(isotropic)       cout << "bool isotropic true. Heisenberg terms will be considered." << endl;
    if(sianisotropy)    cout << "bool sianisotropy true" << endl;
    if(magfield)        cout << "bool magfield true" << endl;
    if(dm)              cout << "bool dm true" << endl;

    this->filenamePrefix = filenamePrefix; // Do I need filenamePrefix any other places than here?
    this->printeveryMCstep = printeveryMCstep;
    this->calculatespincorrelationfunction = calculatespincorrelationfunction;

    if(periodic)    notperiodic = false;
    else            notperiodic = true;

    //if(sianisotropy)    cout << "Siiiii are the people!" << endl;

    // Making the Lattice
    double starttime = clock();
    mylattice = Lattice(L,isotropic, sianisotropy, magfield, dm);
    mylattice.setstrengths(sitestrengthsin, heisenbergin, dm_in);
    cout << "Instance of class Lattice initialized" << endl;

    // Type of Lattice
    if(periodic)
    {
        if(type_lattice=='F')      mylattice.fcc_helical_initialize();                // F for fcc
        if(type_lattice=='E')      mylattice.fcc_helical_initialize_extended();       // E for fcc extended
        else if(type_lattice=='C') mylattice.cubic_helical_initialize();              // C for cubic
        else if(type_lattice=='D') mylattice.cubic_helical_initialize_extended();     // D for cubic extended
        else if(type_lattice=='Q') mylattice.quadratic_helical_initialize();          // Q for quadratic
        else if(type_lattice=='R') mylattice.quadratic_helical_initialize_extended(); // R for quadr. ext.
        else if(type_lattice=='O') mylattice.chain_periodic_initialize();             // O for one-dimensional
    }
    else
    {
        if(type_lattice=='O')      mylattice.chain_open_initialize();
        else                       cout << "WARNING! type_lattice " << type_lattice << " only periodic. You have asked for open BCs! Failure! Failure!" << endl;
    }
    double endtime = clock();
    double total_time = (endtime - starttime)/(double) CLOCKS_PER_SEC;
    cout << "Lattice set, time: " << total_time << endl;

    // The lattice  // Could just extract the bools from here...
    this->N = mylattice.N;
    this->no_of_neighbours = mylattice.no_of_neighbours;
    cout << "Number of sites and neighbours retrieved to MonteCarlo." << endl;

    // Random generators
    // Should clean up in these...
    distribution_u    = std::uniform_real_distribution<double>(0,1);  // Varies depending on distribution
    distribution_v    = std::uniform_real_distribution<double>(0,1);  // Varies depending on distribution
    distribution_prob = std::uniform_real_distribution<double>(0,1);
    // For index. This is given helical boundary conditions, then I only need one index
    distribution_n    = std::uniform_int_distribution<int>(0,N-1);

    cout << "Distributions set" << endl;

    // Initializing some other quantities
    acceptancerate = 0;
    seed1 = 59;  // Seed to start random number generator
    seed2 = 61;
    testseed = 29;
    DEBUG = false;
    MAJORDEBUG = false;

    // Setting up files to print to. Might want to allow more flexibility here
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_everybeta.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    allFile.open(filename);
    delete filename;

    cout << "allFile set" << endl;
    if(printeveryMCstep)
    {
        char *filename = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename, "%s_everyMCstep.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        bigFile.open(filename);
        delete filename;
    }
    cout << "Passed if-test printeveryMCstep" << endl;
    if(calculatespincorrelationfunction)
    {
        cout << "In if-test calculatespincorrelationfunction" << endl;
        // Setting up files to print to
        // Should probably add yet another switch...
        char *filenamex = new char[1000];                                // File name can have max 1000 characters
        sprintf(filenamex, "%s_spincorrelationfunctionx.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        spcorFilex.open(filenamex);
        delete filenamex;

        char *filenamey = new char[1000];                                // File name can have max 1000 characters
        sprintf(filenamey, "%s_spincorrelationfunctiony.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        spcorFiley.open(filenamey);
        delete filenamey;

        char *filenamez = new char[1000];                                // File name can have max 1000 characters
        sprintf(filenamez, "%s_spincorrelationfunctionz.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        spcorFilez.open(filenamez);
        delete filenamez;

        char *filenametot = new char[1000];                                // File name can have max 1000 characters
        sprintf(filenametot, "%s_spincorrelationfunctiontot.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        spcorFiletot.open(filenametot);
        delete filenametot;

        /*
        char *filename2 = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename2, "%s_spincorrelationfunction_transformed.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        ftspcorFile.open(filename2);
        delete filename2;
        */

        cout << "File set" << endl;
    }
    if(randomtest)
    {
        char *filename = new char[1000];                                     // File name can have max 1000 characters
        sprintf(filename, "%s_testofrandom.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        randomtestFile.open(filename);
        delete filename;
    }
}



MonteCarlo::MonteCarlo(int L1, int L2, int L3, int eqsteps, int mcsteps_inbin, int no_of_bins, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    randomtest = false;

    // Handling runningints (wouldn't want to vary this in one class instance, I guess.)
    this->eqsteps = eqsteps;
    this->mcsteps_inbin = mcsteps_inbin;
    this->no_of_bins = no_of_bins;

    // Setting the bools
    this->isotropic = isotropic;
    this->sianisotropy = sianisotropy;
    this->magfield = magfield;
    this->dm = dm;

    this->filenamePrefix = filenamePrefix; // Do I need filenamePrefix any other places than here?
    this->printeveryMCstep = printeveryMCstep;
    this->calculatespincorrelationfunction = calculatespincorrelationfunction;

    if(periodic)    notperiodic = false;
    else            notperiodic = true;

    //if(sianisotropy)    cout << "Siiiii are the people!" << endl;

    // Making the Lattice
    double starttime = clock();
    mylattice = Lattice(L1, L2, L3, isotropic, sianisotropy, magfield, dm);
    mylattice.setstrengths(sitestrengthsin, heisenbergin, dm_in);
    cout << "Instance of class Lattice initialized" << endl;

    // Type of Lattice
    if(periodic)
    {
        if(type_lattice=='F')      mylattice.fcc_helical_initialize();                // F for fcc
        if(type_lattice=='E')      mylattice.fcc_helical_initialize_extended();       // E for fcc extended
        else if(type_lattice=='C') mylattice.cubic_helical_initialize();              // C for cubic
        else if(type_lattice=='D') mylattice.cubic_helical_initialize_extended();     // D for cubic extended
        else if(type_lattice=='Q') mylattice.quadratic_helical_initialize();          // Q for quadratic
        else if(type_lattice=='R') mylattice.quadratic_helical_initialize_extended(); // R for quadr. ext.
        else if(type_lattice=='O') mylattice.chain_periodic_initialize();             // O for one-dimensional
    }
    else
    {
        if(type_lattice=='O')      mylattice.chain_open_initialize();
        else                       cout << "WARNING! type_lattice " << type_lattice << " only periodic. You have asked for open BCs! Failure! Failure!" << endl;
    }

    //else if(type_lattice=='T') mylattice.chain_2p_periodic_initialize(); // T for two particles
    double endtime = clock();
    double total_time = (endtime - starttime)/(double) CLOCKS_PER_SEC;
    cout << "Lattice set, time: " << total_time << endl;

    // The lattice  // Could just extract the bools from here...
    this->N = mylattice.N;
    this->no_of_neighbours = mylattice.no_of_neighbours;
    cout << "Number of sites and neighbours retrieved to MonteCarlo." << endl;

    // Random generators
    // Should clean up in these...
    distribution_u    = std::uniform_real_distribution<double>(0,1);  // Varies depending on distribution
    distribution_v    = std::uniform_real_distribution<double>(0,1);  // Varies depending on distribution
    distribution_prob = std::uniform_real_distribution<double>(0,1);
    // For index. This is given helical boundary conditions, then I only need one index
    distribution_n    = std::uniform_int_distribution<int>(0,N-1);

    cout << "Distributions set" << endl;

    // Initializing some other quantities
    acceptancerate = 0;
    seed1 = 59;  // Seed to start random number generator
    seed2 = 63;
    DEBUG = false;
    MAJORDEBUG = false;

    // Setting up files to print to. Might want to allow more flexibility here
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_everybeta.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    allFile.open(filename);
    delete filename;

    cout << "allFile set" << endl;
    if(printeveryMCstep)
    {
        char *filename = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename, "%s_everyMCstep.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        bigFile.open(filename);
        delete filename;
    }
    cout << "Passed if-test printeveryMCstep" << endl;
    if(calculatespincorrelationfunction)
    {
        cout << "In if-test calculatespincorrelationfunction" << endl;
        // Setting up files to print to
        // Should probably add yet another switch
        char *filenamex = new char[1000];                                // File name can have max 1000 characters
        sprintf(filenamex, "%s_spincorrelationfunctionx.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        spcorFilex.open(filenamex);
        delete filenamex;

        char *filenamey = new char[1000];                                // File name can have max 1000 characters
        sprintf(filenamey, "%s_spincorrelationfunctiony.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        spcorFiley.open(filenamey);
        delete filenamey;

        char *filenamez = new char[1000];                                // File name can have max 1000 characters
        sprintf(filenamez, "%s_spincorrelationfunctionz.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        spcorFilez.open(filenamez);
        delete filenamez;

        char *filenametot = new char[1000];                                // File name can have max 1000 characters
        sprintf(filenametot, "%s_spincorrelationfunctiontot.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        spcorFiletot.open(filenametot);
        delete filenametot;

        /*
        char *filename2 = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename2, "%s_spincorrelationfunction_transformed.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        ftspcorFile.open(filename2);
        delete filename2;
        */

        cout << "File set" << endl;  

    }
    if(randomtest)
    {
        char *filename = new char[1000];                                     // File name can have max 1000 characters
        sprintf(filename, "%s_testofrandom.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        randomtestFile.open(filename);
        delete filename;
    }
}

void MonteCarlo::writeallqstofile()
{
    // Printing information about the q-vectors straight away
    char *filenameq = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenameq, "%s_qs.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    qFile.open(filenameq);
    delete filenameq;

    cout << "File for printing q-vector to file is initiated" << endl;

    double qx;
    double qy;
    double qz;

    int L1 = mylattice.L1;
    int L2 = mylattice.L2;

    if(mylattice.dim==1) // Not really neccessary. Not double tested, either
    {
        for(int i=0; i<N; i++)    qFile << i << " " << 2*M_PI/L1*i << endl;
    }
    if(mylattice.dim==2)
    {

        cout << "For quadratic, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns = mylattice.sitecoordinates[i];
            //cout << "ns retrieved" << endl;
            qx = ns[0]*mylattice.b1[0]/L1; // Not so sure about this, just need it to compile
            qy = ns[1]*mylattice.b2[1]/L2; // Double check
            //cout << "qvec set" << endl;
            qFile << i << " " << qx << " " << qy << " " << endl;

        }
        cout << "Done printing to qFile" << endl;

    }

    else if(mylattice.dim==3)
    {
        int L3 = mylattice.L3;
        cout << "For cubic or fcc, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns = mylattice.sitecoordinates[i];
            //cout << "ns retrieved" << endl;
            // These must be changed if we change into possibly setting L1, L2, L3 different
            // Don't really need b1, b2, b3, could just use a1, a2, a3 multiplied by 2*M_PI...
            // Could be more general, but we don't need to play around with our lattices that much...
            qx = ns[0]*mylattice.b1[0]/L1 + ns[2]*mylattice.b3[0]/L3;
            qy = ns[0]*mylattice.b1[1]/L1 + ns[1]*mylattice.b2[1]/L2;
            qz = ns[1]*mylattice.b2[2]/L2 + ns[2]*mylattice.b3[2]/L3;
            // Print to file. Site number, qx, qy, qz.
            qFile << i << " " << qx << " " << qy << " " << qz << endl;
        }
        //cout << "Done printing to qFile" << endl;
    }
}



/*
void MonteCarlo::chooseprintfile(string filenamePrefix)
{
    // Initializing class for printing to file
    print(filenamePrefix);
    //print(filenamePrefix);
    //print.givePrefix(filenamePrefix);
    // For now: Just print to all files. Have other options, of course.
    print.open_all();
}
*/

void MonteCarlo::debugmode(bool on)
{
    if(on)    DEBUG = true;
    else      DEBUG = false;
}

void MonteCarlo::majordebugtrue()
{
    MAJORDEBUG = true;

    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_comparetheory_results.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    compareFile.open(filename);
    delete filename;

}

void MonteCarlo::initialize_energy()
{
    if(DEBUG)    cout << "In MonteCarlo_initialize_energy" << endl;
    energy_old = 0;
    double energy_contribution_sites = 0;
    double energy_contribution_bonds = 0;

    bool BEDBUG = false;
    for(int i=0; i<N; i++)
    {
        if(BEDBUG) cout << "In loop to set initial energy" << endl;
        // Spin i
        double sx = mylattice.sites[i].spinx;
        double sy = mylattice.sites[i].spiny;
        double sz = mylattice.sites[i].spinz;
        // Contribution from sites
        if(sianisotropy)
        {
            //cout << "in sianisotropy, init" << endl;
            if(BEDBUG)    cout << "In sianisotropy" << endl;
            double Dix = mylattice.sites[i].Dix;   // Should these lie in Lattice instead of in Lattice::site?
            double Diy = mylattice.sites[i].Diy;
            double Diz = mylattice.sites[i].Diz;           
            energy_contribution_sites += (Dix*sx*sx + Diy*sy*sy+ Diz*sz*sz);
            //cout << "Energy contribution from sites: " << energy_contribution_sites << endl;
        }
        if(magfield)
        {
            if(BEDBUG)    cout << "In magfield" << endl;
            double hx = mylattice.sites[i].hx;
            double hy = mylattice.sites[i].hy;
            double hz = mylattice.sites[i].hz;
            energy_contribution_sites -= hx*sx + hy*sy + hz*sz;
        }
        // Contribution from bonds
        if(isotropic)
        {
            if(BEDBUG)   cout << "In isotropic" << endl;
            double nneighbours;
            double partnerspinx = 0;
            double partnerspiny = 0;
            double partnerspinz = 0;
            // Determining the number of neighbours of the site
            if(notperiodic)    nneighbours = mylattice.sites[i].no_of_neighbours_site;
            else               nneighbours = no_of_neighbours;
            //cout << "Looping over " << nneighbours << " neighbours" << endl;
            for(int j=0; j<nneighbours; j++)
            {
                if(BEDBUG)    cout << "in loop in isotropic, j = " << j << endl;
                int k = mylattice.sites[i].bonds[j].siteindex2;
                if(i<k)    // To avoid double counting. More computationally efficient than halving the energy.
                {
                    if(BEDBUG)    cout << "Have set k in loop, k = " << k << "; N = " << N << endl;
                    double J = mylattice.sites[i].bonds[j].J;
                    if(BEDBUG)    cout << "Have accessed J in bond between spins" << endl;
                    if(J!=0)
                    {
                        double sxk = mylattice.sites[k].spinx;
                        double syk = mylattice.sites[k].spiny;
                        double szk = mylattice.sites[k].spinz;
                        if(BEDBUG)    cout << "Have accessed the components of the spin on the other end" << endl;
                        partnerspinx += J*sxk;
                        partnerspiny += J*syk;
                        partnerspinz += J*szk;
                        if(BEDBUG)    cout << "Have gathered this contribution into partnerspin" << endl;
                    }
                }
            }
            if(BEDBUG)   cout << "Done with the loop in isotropic" << endl;
            if(BEDBUG)   cout << "partnerspinx = " << partnerspinx << "; partnerspiny = " << partnerspiny << "; partnerspinz = " << partnerspinz << endl;
            energy_contribution_bonds += (partnerspinx*sx + partnerspiny*sy + partnerspinz*sz);
            //cout << "Energy contribution from bonds: " << energy_contribution_bonds << endl;
            if(BEDBUG)     cout << "Done with isotropic" << endl;
        }
        if(dm)
        {
            if(BEDBUG)    cout << "In dm" << endl;
            // Double loops and stuff. Could maybe make this more efficient
            int nneighbours;
            // Determining the number of neighbours of the site
            if(notperiodic)    nneighbours = mylattice.sites[i].no_of_neighbours_site;
            else               nneighbours = no_of_neighbours;
            for(int j=0; j<nneighbours; j++)
            {
                int k = mylattice.sites[i].bonds[j].siteindex2; // Hope I can actually get to this value.
                if(i<k)    // To avoid double counting. Hopefully saves time for large systems
                {
                    double Dx = mylattice.sites[i].bonds[j].Dx;
                    double Dy = mylattice.sites[i].bonds[j].Dy;
                    double Dz = mylattice.sites[i].bonds[j].Dz;

                    double sxk = mylattice.sites[k].spinx;
                    double syk = mylattice.sites[k].spiny;
                    double szk = mylattice.sites[k].spinz;

                    energy_contribution_bonds += Dx*(sy*szk-syk*sz)+Dy*(sz*sxk-szk*sx)+Dz*(sx*syk-sy*sxk);
                }
            }
        }
        if(BEDBUG) cout << "Done with one, onto the others" << endl;
    }
    energy_old = energy_contribution_sites + energy_contribution_bonds;
    if(BEDBUG)    cout << "Energy initialized, energy = " << energy_old << endl;
}

void MonteCarlo::reset_energy()
{  // In case it gets stuck in some region, I guess...
    //cout << "In reset_energy()" << endl;
    bool randomspins = true;
    double originalvalue = 1/sqrt(3);
    for(int i=0; i<N; i++)
    {   // Reset all spins
        mylattice.sites[i].spinx = originalvalue;
        mylattice.sites[i].spiny = originalvalue;
        mylattice.sites[i].spinz = originalvalue;
        if(randomspins)
        {
            double u = ran2(&seed1);
            double v = ran2(&seed1);

            double theta = acos(1.0-2.0*u);
            double phi = 2.0*M_PI*v;

            double sintheta = sin(theta);
            double spinx = sintheta*cos(phi);
            double spiny = sintheta*sin(phi);
            double spinz = cos(theta);

            mylattice.sites[i].spinx = spinx;
            mylattice.sites[i].spiny = spiny;
            mylattice.sites[i].spinz = spinz;
        }
    }

    initialize_energy();
    //setrandomgenerators();
}

void MonteCarlo::givexplanforFFT(vector<double>& r, vector<complex<double> >& q)  // Return p? Or have p as a class variable?
{
    int rank = mylattice.dim;               // Dimension of lattice
    vector<int> Ls = mylattice.dimlengths;  // List containing [L], [L1,L2], [L1,L2,L3],
                                            // depending on the lattice
    // p declared as a class variable
    px = fftw_plan_dft_r2c(rank,
                          &Ls[0],
                          &r[0],
                          reinterpret_cast<fftw_complex*>(&q[0]),
                          FFTW_ESTIMATE);
}

void MonteCarlo::giveyplanforFFT(vector<double>& r, vector<complex<double> >& q)  // Return p? Or have p as a class variable?
{
    int rank = mylattice.dim;               // Dimension of lattice
    vector<int> Ls = mylattice.dimlengths;  // List containing [L], [L1,L2], [L1,L2,L3],
                                            // depending on the lattice
    // p declared as a class variable
    py = fftw_plan_dft_r2c(rank,
                          &Ls[0],
                          &r[0],
                          reinterpret_cast<fftw_complex*>(&q[0]),
                          FFTW_ESTIMATE);
}

void MonteCarlo::givezplanforFFT(vector<double>& r, vector<complex<double> >& q)  // Return p? Or have p as a class variable?
{
    int rank = mylattice.dim;               // Dimension of lattice
    vector<int> Ls = mylattice.dimlengths;  // List containing [L], [L1,L2], [L1,L2,L3],
                                            // depending on the lattice
    // p declared as a class variable
    pz = fftw_plan_dft_r2c(rank,
                          &Ls[0],
                          &r[0],
                          reinterpret_cast<fftw_complex*>(&q[0]),
                          FFTW_ESTIMATE);
}


void MonteCarlo::giveplanforFFT_inverse(vector<double>& rout, vector<complex<double> >& q)  // Return p? Or have p as a class variable?
{
    int rank = mylattice.dim;               // Dimension of lattice
    vector<int> Ls = mylattice.dimlengths;  // List containing [L], [L1,L2], [L1,L2,L3],
                                            // depending on the lattice
    // p declared as a class variable
    pinv = fftw_plan_dft_c2r(rank,
                          &Ls[0],
                          reinterpret_cast<fftw_complex*>(&q[0]),
                          &rout[0],
                          FFTW_ESTIMATE);
}


void MonteCarlo::runmetropolis(double beta)
{
    if(DEBUG)    cout << "In runmetropolis in MonteCarlo" << endl;
    bool HUMBUG  = false;
    bool LADYBUG = false;
    bool SCBUG   = false; // For finding out what's wrong with the spin correlation function
    bool DMEBUG  = false;
    bool ENGYBUG = true; // For checking that our two ways of computing the energy agrees (somewhat)
    bool bincout = false;

    // Header for spcorFiles
    if(calculatespincorrelationfunction)
    {
        spcorFilex   << beta << " " << N << endl;
        spcorFiley   << beta << " " << N << endl;
        spcorFilez   << beta << " " << N << endl;
        spcorFiletot << beta << " " << N << endl;
        ftspcorFile  << beta << " " << N << endl;
    }
    //allFile << "N: " << N << "; Jxz: " << mylattice.sites[0].bonds[0].J << "; Jyz: " << mylattice.sites[0].bonds[2].J << "; Jxy: " << mylattice.sites[0].bonds[4].J << "; Dix: " << mylattice.sites[0].Dix << "; Diy: " << mylattice.sites[0].Diy << "; Diz: " << mylattice.sites[0].Diz << endl;
    //allFile << "eqsteps: " << eqsteps << "; mcsteps_inbin: " << mcsteps_inbin << "; no_of_bins: " << no_of_bins << endl;

    // Initializing the energy
    initialize_energy();

    if(DEBUG)    cout << "Done with initialize_energy()" << endl << "Starting equilibration steps" << endl;
    cout << "Initial energy: " << energy_old << endl;
    cout << "Initial energy, check: " << check_the_energy() << endl;

    // Equilibration steps
    double starttime = clock();
    for(int i=0; i<eqsteps; i++)
    {
        if(HUMBUG)    cout << "In equilibration steps loop, i = " << i << endl;
        mcstepf_metropolis(beta); //, generator_u, generator_v, generator_n, generator_prob, distribution_prob, distribution_u, distribution_v, distribution_n);
        if(i<11)      cout << "i = " << i << ", energy: " << energy_old << endl;
        if(i<11)      cout << "Hardcoded energy, step " << i << ": " << check_the_energy() << endl;
    }
    double endtime = clock();
    double total_time = (endtime - starttime)/(double) CLOCKS_PER_SEC;
    cout << "Time equilibration steps: " << total_time << endl;


    if(DEBUG)    cout << "Done with equilibration steps" << endl << "Starting the Monte Carlo steps and measurements" << endl;

    // Monte Carlo steps and measurements
    starttime = clock();
    // Measurable quantities
    std::vector<double> acceptancerates    = std::vector<double>(no_of_bins);
    std::vector<double> energies           = std::vector<double>(no_of_bins);
    std::vector<double> energies_sq        = std::vector<double>(no_of_bins);
    std::vector<double> cvs                = std::vector<double>(no_of_bins);
    std::vector<double> mxs                = std::vector<double>(no_of_bins);
    std::vector<double> mys                = std::vector<double>(no_of_bins);
    std::vector<double> mzs                = std::vector<double>(no_of_bins);
    std::vector<double> mxs_abs            = std::vector<double>(no_of_bins);
    std::vector<double> mys_abs            = std::vector<double>(no_of_bins);
    std::vector<double> mzs_abs            = std::vector<double>(no_of_bins);
    std::vector<double> mxssq              = std::vector<double>(no_of_bins);
    std::vector<double> myssq              = std::vector<double>(no_of_bins);
    std::vector<double> mzssq              = std::vector<double>(no_of_bins);
    std::vector<double> mxsquad            = std::vector<double>(no_of_bins);
    std::vector<double> mysquad            = std::vector<double>(no_of_bins);
    std::vector<double> mzsquad            = std::vector<double>(no_of_bins);
    // For the correlation function
    // The spins
    std::vector<double> spins_in_x = std::vector<double>(N);
    std::vector<double> spins_in_y = std::vector<double>(N);
    std::vector<double> spins_in_z = std::vector<double>(N);
    // Declare qconf and set the plan
    vector< complex<double> > qconfx(N);  // Output array
    vector< complex<double> > qconfy(N);  // Output array
    vector< complex<double> > qconfz(N);  // Output array
    givexplanforFFT(spins_in_x, qconfx);
    giveyplanforFFT(spins_in_y, qconfy);
    givezplanforFFT(spins_in_z, qconfz);

    // Array for the results
    // Determining the length of the array
    int dim = mylattice.dim;
    int qlimit = 1; // To be multiplied;
    for(int l=0; l<(dim-1); l++)
    {    // Looping over all dimensions but the last
        qlimit *= mylattice.dimlengths[l];
    }
    qlimit *= mylattice.dimlengths[dim-1]/2+1;

    // Do the inverse?
    std::vector<double> rout = std::vector<double>(N); // Output array
    vector< complex<double> > qconfstore(qlimit);  // Array to store spin correlation function for input
    giveplanforFFT_inverse(rout, qconfstore);

    vector<double> correlation_functionx_av_bin    = vector<double>(N);
    vector<double> correlation_functiony_av_bin    = vector<double>(N);
    vector<double> correlation_functionz_av_bin    = vector<double>(N);
    vector<double> correlation_functiontot_av_bin  = vector<double>(N);
    vector<double> correlation_functionx_av        = vector<double>(N);
    vector<double> correlation_functiony_av        = vector<double>(N);
    vector<double> correlation_functionz_av        = vector<double>(N);
    vector<double> correlation_functiontot_av      = vector<double>(N);
    for(int i=0; i<N; i++)    correlation_functionx_av[i]   = 0;
    for(int i=0; i<N; i++)    correlation_functiony_av[i]   = 0;
    for(int i=0; i<N; i++)    correlation_functionz_av[i]   = 0;
    for(int i=0; i<N; i++)    correlation_functiontot_av[i] = 0;
    vector<vector<double> > correlation_functionx_store;
    vector<vector<double> > correlation_functiony_store;
    vector<vector<double> > correlation_functionz_store;
    vector<vector<double> > correlation_functiontot_store;

    // For the inverse Fourier transform (after we have taken fftw)
    // This should probably be done to the correlation function, not spins_in_z...
    // Because otherwise we would just get the same result back (though not normalized...)
    vector<double> ftcorrelation_function_av_bin    = vector<double>(N);
    vector<double> ftcorrelation_function_av        = vector<double>(N);
    for(int i=0; i<N; i++)    ftcorrelation_function_av[i] = 0;
    vector<vector<double> > ftcorrelation_function_store;

    // Resetting quantities
    double ar_av        = 0;
    double energy_av    = 0;
    double energy_sq_av = 0;
    double cv_average   = 0;
    double mx_av        = 0;
    double my_av        = 0;
    double mz_av        = 0;
    double mx_abs_av    = 0;
    double my_abs_av    = 0;
    double mz_abs_av    = 0;
    double mxsq_av      = 0;
    double mysq_av      = 0;
    double mzsq_av      = 0;
    double mxquad_av    = 0;
    double myquad_av    = 0;
    double mzquad_av    = 0;
    for(int i=0; i<no_of_bins; i++)  // Loop over the bins
    {   // For every bin
        if(SCBUG)    cout << "Starting bin " << i << endl;
        // Reset quantities
        energies[i]    = 0;
        energies_sq[i] = 0;
        cvs[i]         = 0;
        mxs[i]         = 0;
        mys[i]         = 0;
        mzs[i]         = 0;
        mxs_abs[i]     = 0;
        mys_abs[i]     = 0;
        mzs_abs[i]     = 0;
        mxssq[i]       = 0;
        myssq[i]       = 0;
        mzssq[i]       = 0;
        mxsquad[i]     = 0;
        mysquad[i]     = 0;
        mzsquad[i]     = 0;
        // Resetting the correlation function bin average for every bin
        for(int k=0; k<N; k++)    correlation_functionx_av_bin[k]     = 0;
        for(int k=0; k<N; k++)    correlation_functiony_av_bin[k]     = 0;
        for(int k=0; k<N; k++)    correlation_functionz_av_bin[k]     = 0;
        for(int k=0; k<N; k++)    correlation_functiontot_av_bin[k]   = 0;
        for(int k=0; k<N; k++)    ftcorrelation_function_av_bin[k]    = 0;
        if(LADYBUG)
        {
            if(i>0)
            {
                cout << "i = " << i << "; Average energy before loop: " << energy_av/(mcsteps_inbin*i) << endl;
                cout << "i = " << i << "; energies[i]: " << energies[i] << endl;
            }
        }

        // Setting vectors

        //cout << "Going to start the MCsteps" << endl;
        for(int j=0; j<mcsteps_inbin; j++)    // Loop over mcsteps in bin
        {   // For each mcstep
            if(SCBUG)    cout << "in mcstepsloop, j= " << j << endl;
            mcstepf_metropolis(beta);
            if(SCBUG)    cout << "mcstep done" << endl;
            if(randomtest)
            {
                randomtestFile << ran2(&testseed) << endl;
            }

            // acceptancerate
            acceptancerates[i] += acceptancerate;
            ar_av += acceptancerate;
            // energy
            //cout << "Current energy: " << energy_old << endl;
            energies[i]    += energy_old;    // Storing to get the standard deviation
            energies_sq[i] += energy_old*energy_old;
            if(DMEBUG) // Test to fix DM energy problem.
            {
                double energy_hardcoded = check_the_energy();
                cout << "Energy, MC routine: " << energy_old << "; energy, read off the spins: " << energy_hardcoded << endl;
            }
            if(ENGYBUG)
            {
                double energy_hardcoded = check_the_energy();
                double endiff = energy_hardcoded-energy_old;
                // Choose when I would like to be notified
                if(abs(endiff)>1e-9)    cout << "A (small?) difference between hardcoded and derived energy: endiff = " << endiff << endl;
                //cout << "A (small?) difference between hardcoded and derived energy: endiff = " << endiff << endl;

            }
            //energy_av      += energy_old;
            //energy_sq_av   += energy_old*energy_old;
            //cout << "Current energy: " << energy_old << "; Average energy so far: " << energies[i]/(j+1) << endl;
            // Magnetization
            double mx     = 0;
            double my     = 0;
            double mz     = 0;
            double mx_abs = 0;
            double my_abs = 0;
            double mz_abs = 0;
            double mxsq   = 0;
            double mysq   = 0;
            double mzsq   = 0;
            double mxquad = 0;
            double myquad = 0;
            double mzquad = 0;
            if(SCBUG)    cout << "Retrieving spins to feed to FFTW" << endl;
            for(int k=0; k<N; k++)
            {
                mx+= mylattice.sites[k].spinx;
                my+= mylattice.sites[k].spiny;
                mz+= mylattice.sites[k].spinz;
                // For the correlation function
                if(calculatespincorrelationfunction)
                {
                    spins_in_x[k] = mylattice.sites[k].spinx;
                    spins_in_y[k] = mylattice.sites[k].spiny;
                    spins_in_z[k] = mylattice.sites[k].spinz;
                }
            }
            if(SCBUG)    cout << "Done with that, ascribing magnetizations and such" << endl;
            mx = mx/N;
            my = my/N;
            mz = mz/N;
            mx_abs = abs(mx);
            my_abs = abs(my);
            mz_abs = abs(mz);
            mxsq = mx*mx;
            mysq = my*my;
            mzsq = mz*mz;
            mxquad = mxsq*mxsq;
            myquad = mysq*mysq;
            mzquad = mzsq*mzsq;

            mx_av += mx;
            my_av += my;
            mz_av += mz;
            mx_abs_av += mx_abs;
            my_abs_av += my_abs;
            mz_abs_av += mz_abs;
            mxsq_av += mxsq;
            mysq_av += mysq;
            mzsq_av += mzsq;
            mxquad_av += mxquad;
            myquad_av += myquad;
            mzquad_av += mzquad;

            mxs[i] += mx;
            mys[i] += my;
            mzs[i] += mz;
            mxs_abs[i] += mx_abs;
            mys_abs[i] += my_abs;
            mzs_abs[i] += mz_abs;
            mxssq[i] += mxsq;
            myssq[i] += mysq;
            mzssq[i] += mzsq;
            mxsquad[i] += mxquad;
            mysquad[i] += myquad;
            mzsquad[i] += mzquad;

            // FFT steps
            // Should I do this here? The manual said that we can reuse the plan.
            if(calculatespincorrelationfunction)
            {
                if(SCBUG)    cout << "Calculating the spin correlation function in mcsteps" << endl;
                fftw_execute(px);
                fftw_execute(py);
                fftw_execute(pz);

                // Back again:
                //Finding <S^z_{q}S^z_{-q}> to look at the spin correlation as a function of distance
                for(int n=0; n<qlimit; n++)    qconfstore[n] = qconfz[n]*conj(qconfz[n]); // This is unnormalized
                fftw_execute(pinv);

                double cx, cy, cz;
                if(mylattice.dim==1) // Chain
                {
                    int index;
                    for(int n=0; n<N; n++)
                    {
                        if(n<=(int)N/2)    index = n;
                        else               index = N-n;
                        cx = (qconfx[index]*conj(qconfx[index])).real()/(N*N);
                        cy = (qconfy[index]*conj(qconfy[index])).real()/(N*N);
                        cz = (qconfz[index]*conj(qconfz[index])).real()/(N*N);

                        correlation_functionx_av_bin[n]   += cx;
                        correlation_functiony_av_bin[n]   += cy;
                        correlation_functionz_av_bin[n]   += cz;
                        correlation_functiontot_av_bin[n] += cx + cy + cz;

                        // Normalization constant?
                        double normconst = 1/(N*N); // Only a factor of 1/N one way
                        for(int n=0; n<N ;n++)   ftcorrelation_function_av_bin[n] += normconst*rout[n];
                    }
                }
                if(mylattice.dim==2)
                {
                    int L1 = mylattice.L1;
                    int L2 = mylattice.L2;
                    int elinar = 0; // For retrieving the elements residing in the output array
                    for(int n=0; n<N; n++)
                    {
                        vector<int> cord = mylattice.sitecoordinates[n];
                        int n1 = cord[0];
                        int n2 = cord[1];
                        int index;
                        if(n2<=(int)L2/2)
                        {   // If our element is stored in the output array, we retrieve it
                            index = elinar;
                            elinar++;        // We move one index forward in the output array
                            //cout << "I have met ms elinar" << endl;
                        }
                        else
                        {   // If our element is not stored in the output array, we retrieve its
                            // complex conjugate. We needn't do anything with it as the result is a
                            // complex number times its complex conjugate
                            index = ((int)L2/2+1)*((L1-n1)%L1)+(L2-n2)%L2;
                            //cout << "Retrieving the index of the cc" << endl;
                        } // End if-tests
                        if(SCBUG)    cout << "Adding to the spin correlation function" << endl;
                        cx = (qconfx[index]*conj(qconfx[index])).real()/(N*N);
                        cy = (qconfy[index]*conj(qconfy[index])).real()/(N*N);
                        cz = (qconfz[index]*conj(qconfz[index])).real()/(N*N);

                        correlation_functionx_av_bin[n]   += cx;
                        correlation_functiony_av_bin[n]   += cy;
                        correlation_functionz_av_bin[n]   += cz;
                        correlation_functiontot_av_bin[n] += cx + cy + cz;

                        // Change this:

                        // Normalization constant?
                        double normconst = 1/(N*N); // Only a factor of 1/N one way
                        for(int n=0; n<N ;n++)   ftcorrelation_function_av_bin[n] += normconst*rout[n];
                    }
                }
                if(mylattice.dim==3) // Simple cubic or fcc
                {
                    // Patching together the spin correlation function to look at <S^z_{q}S^z_{-q}>
                    if(SCBUG)    cout << "Patching the sp.corr. output" << endl;
                    //cout << "Our dimension is 3" << endl;
                    int L1 = mylattice.L1;
                    int L2 = mylattice.L2;
                    int L3 = mylattice.L3;
                    int elinar = 0; // For retrieving the elements residing in the output array
                    for(int n=0; n<N; n++)
                    {
                        vector<int> cord = mylattice.sitecoordinates[n];
                        int n1 = cord[0];
                        int n2 = cord[1];
                        int n3 = cord[2];
                        int index;
                        if(n3<=(int)L3/2)
                        {   // If our element is stored in the output array, we retrieve it
                            index = elinar;
                            elinar++;        // We move one index forward in the output array
                            //cout << "I have met ms elinar" << endl;
                        }
                        else
                        {   // If our element is not stored in the output array, we retrieve its
                            // complex conjugate. We needn't do anything with it as the result is a
                            // complex number times its complex conjugate
                            index = L2*((int)L3/2+1)*((L1-n1)%L1)+((int)L3/2+1)*((L2-n2)%L2)+((L3-n3)%L3);
                            //cout << "Retrieving the index of the cc" << endl;
                        } // End if-tests
                        if(SCBUG)    cout << "Adding to the spin correlation function" << endl;
                        cx = (qconfx[index]*conj(qconfx[index])).real()/(N*N);
                        cy = (qconfy[index]*conj(qconfy[index])).real()/(N*N);
                        cz = (qconfz[index]*conj(qconfz[index])).real()/(N*N);

                        correlation_functionx_av_bin[n]   += cx;
                        correlation_functiony_av_bin[n]   += cy;
                        correlation_functionz_av_bin[n]   += cz;
                        correlation_functiontot_av_bin[n] += cx + cy + cz;
                        if(SCBUG)    cout << "Done adding to the spin correlation function" << endl;
                        if(SCBUG)    cout << "Which is: " << correlation_functionz_av_bin[n] << endl;
                    } // End loop over n
                    if(SCBUG)    cout << "Done patching together the corr.func.array" << endl;

                    // Normalization constant?
                    double normconst = 1/(N*N); // Only a factor of 1/N one way
                    for(int n=0; n<N ;n++)   ftcorrelation_function_av_bin[n] += normconst*rout[n];
                }
            }

            //Print to bigFile
            if(printeveryMCstep)
            {
                bigFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << beta << " " << energy_old << " " << energy_old*energy_old << " " << mx << " " << my << " " << mz << endl;
            }

            // Some sort of measurement of the magnetization... How to do this when we have a continuous spin?
        }  // End loop over mcsteps
        if(SCBUG)    cout << "Done with mcsteps, bin " << i << endl;

        //cout << "Sum energies, bin " << i << ": " << energies[i] << endl;

        // For every bin, we find the following quantities:
        mxs[i]         = mxs[i]/mcsteps_inbin;
        mys[i]         = mys[i]/mcsteps_inbin;
        mzs[i]         = mzs[i]/mcsteps_inbin;
        mxs_abs[i]     = mxs_abs[i]/mcsteps_inbin;
        mys_abs[i]     = mys_abs[i]/mcsteps_inbin;
        mzs_abs[i]     = mzs_abs[i]/mcsteps_inbin;
        mxssq[i]       = mxssq[i]/mcsteps_inbin;
        myssq[i]       = myssq[i]/mcsteps_inbin;
        mzssq[i]       = mzssq[i]/mcsteps_inbin;
        mxsquad[i]     = mxsquad[i]/mcsteps_inbin;
        mysquad[i]     = mysquad[i]/mcsteps_inbin;
        mzsquad[i]     = mzsquad[i]/mcsteps_inbin;
        energies[i]    = energies[i]/mcsteps_inbin;
        energies_sq[i] = energies_sq[i]/mcsteps_inbin;
        double cv_bin  = beta*beta*(energies_sq[i]-energies[i]*energies[i]);
        cvs[i]         = cv_bin;
        cv_average    += cv_bin;


        //cout << "Average energy, bin " << i << ": " << energies[i] << endl;
        if(SCBUG)    cout << "Have calculated all bin quantities, now doing the spin correlation function" << endl;

        if(bincout)    allFile << "Average energy, bin " << i << ": " << energies[i] << endl;

        if(calculatespincorrelationfunction)
        {   // Take the average and print to file
            // Make the averages (is this more efficient than just calculating the average?)
            //cout << "Looping over all particles to get their spin correlation function" << endl;
            for(int l = 0; l<N; l++)
            {
                //cout << "calculating the average of the spin correlation function in this bin" << endl;
                correlation_functionx_av_bin[l]   /=(mcsteps_inbin);
                correlation_functiony_av_bin[l]   /=(mcsteps_inbin);
                correlation_functionz_av_bin[l]   /=(mcsteps_inbin);
                correlation_functiontot_av_bin[l] /=(mcsteps_inbin);
                //cout << "Average for this bin, particle " << l << ": " << correlation_function_av_bin[l] << endl;
                correlation_functionx_av[l]    += correlation_functionx_av_bin[l]; // I divide by the number of bins later
                correlation_functiony_av[l]    += correlation_functiony_av_bin[l];
                correlation_functionz_av[l]    += correlation_functionz_av_bin[l];
                correlation_functiontot_av[l]  += correlation_functiontot_av_bin[l];
                //cout << "Done calculating average of the spin correlation function over every bin" << endl;
                // Could print for every bin, but is that really neccessary?
                //spcorFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << correlation_function_av[l] << " ";  // Should I include a beta, just in case?

                ftcorrelation_function_av_bin[l] /= mcsteps_inbin;
                ftcorrelation_function_av[l]  += ftcorrelation_function_av_bin[l];
            }
            // WHERE TO PUT THE TOTAL CORRELATION FUNCTION?
            //spcorFile << endl; // End line to get ready for new result
            correlation_functionx_store.push_back(correlation_functionx_av_bin);
            correlation_functiony_store.push_back(correlation_functiony_av_bin);
            correlation_functionz_store.push_back(correlation_functionz_av_bin);

            correlation_functiontot_store.push_back(correlation_functiontot_av_bin);

            ftcorrelation_function_store.push_back(ftcorrelation_function_av_bin);
        } // End if-test calculatecorrelationfunction
        if(SCBUG)    cout << "Have made correlation_function_store, bin " << i << endl;
        endtime = clock();
        total_time = (endtime - starttime)/(double) CLOCKS_PER_SEC;
        if((double)(i+1)/no_of_bins==0.10)    cout << "10% done. Time elapsed: " << total_time << endl;
        if((double)(i+1)/no_of_bins==0.25)    cout << "25% done. Time elapsed: " << total_time << endl;
        if((double)(i+1)/no_of_bins==0.50)    cout << "50% done. time elapsed: " << total_time << endl;
        if((double)(i+1)/no_of_bins==0.75)    cout << "75% done. Time elapsed: " << total_time << endl;
        if((double)(i+1)/no_of_bins==1.00)    cout << "100% done. Time elapsed: " << total_time<< endl;


    }  // End loops over bins
    //cout << "Done with the bins" << endl;
    //----------------------// Acceptance rate //-----------------------//
    ar_av = ar_av/(mcsteps_inbin*no_of_bins);

    //cout << "Acceptance rate std" << endl;
    // Standard deviation //
    double ar_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    ar_stdv +=(acceptancerates[l]-ar_av)*(acceptancerates[l]-ar_av);
    ar_stdv = sqrt(ar_stdv/(no_of_bins*(no_of_bins-1)));

    //----------------------// Energy //-----------------------//
    //energy_av = energy_av/(mcsteps_inbin*no_of_bins);
    //energy_sq_av = energy_sq_av/(mcsteps_inbin*no_of_bins);

    energy_av = 0; // May choose to set this here. Then I need to prefix with double

    //cout << "Average energy" << endl;
    // Average energy
    for(int l=0; l<no_of_bins; l++)    energy_av += energies[l];
    energy_av = energy_av/no_of_bins;

    // Average squared energy
    for(int l=0; l<no_of_bins; l++)    energy_sq_av += energies_sq[l];
    energy_sq_av = energy_sq_av/no_of_bins;

    // May do the same with the spins.

    // Error in the energy //
    //cout << "Energy std" << endl;
    double E_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    E_stdv += (energies[l]-energy_av)*(energies[l]-energy_av);
    E_stdv = sqrt(E_stdv/(no_of_bins*(no_of_bins-1)));

    double Esq_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    Esq_stdv += (energies_sq[l]-energy_sq_av)*(energies_sq[l]-energy_sq_av);
    Esq_stdv = sqrt(Esq_stdv/(no_of_bins*(no_of_bins-1)));

    //---------------------//Heat capacity//---------------------//
    //cout << "Heat capacity" << endl;
    double cv = beta*beta*(energy_sq_av-energy_av*energy_av);

    // Approximate error in the heat capacity:

    cv_average = cv_average/(mcsteps_inbin*no_of_bins);
    //cout << "cv_average: " << cv_average << endl;

    //cout << "Heat capacity std" << endl;
    double cv_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    cv_stdv += (cvs[l]-cv_average)*(cvs[l]-cv_average);
    cv_stdv = sqrt(cv_stdv/(no_of_bins*(no_of_bins-1)));

    //-----------------------//Magnetization//----------------------//
    //cout << "Magnetization" << endl;
    mx_av     = mx_av/(mcsteps_inbin*no_of_bins);
    my_av     = my_av/(mcsteps_inbin*no_of_bins);
    mz_av     = mz_av/(mcsteps_inbin*no_of_bins);
    mx_abs_av = mx_abs_av/(mcsteps_inbin*no_of_bins);
    my_abs_av = my_abs_av/(mcsteps_inbin*no_of_bins);
    mz_abs_av = mz_abs_av/(mcsteps_inbin*no_of_bins);
    mxsq_av   = mxsq_av/(mcsteps_inbin*no_of_bins);
    mysq_av   = mysq_av/(mcsteps_inbin*no_of_bins);
    mzsq_av   = mzsq_av/(mcsteps_inbin*no_of_bins);
    mxquad_av = mxquad_av/(mcsteps_inbin*no_of_bins);
    myquad_av = myquad_av/(mcsteps_inbin*no_of_bins);
    mzquad_av = mzquad_av/(mcsteps_inbin*no_of_bins);

    // Error in the magnetization //
    //cout << "Magnetization std" << endl;
    // First power
    double mx_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    mx_stdv += (mxs[l]-mx_av)*(mxs[l]-mx_av);
    mx_stdv = sqrt(mx_stdv/(no_of_bins*(no_of_bins-1)));

    double my_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    my_stdv += (mys[l]-my_av)*(mys[l]-my_av);
    my_stdv = sqrt(my_stdv/(no_of_bins*(no_of_bins-1)));

    double mz_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    mz_stdv += (mzs[l]-mz_av)*(mzs[l]-mz_av);
    mz_stdv = sqrt(mz_stdv/(no_of_bins*(no_of_bins-1)));

    // Absoulte value
    double mx_abs_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    mx_abs_stdv += (mxs_abs[l]-mx_abs_av)*(mxs_abs[l]-mx_abs_av);
    mx_abs_stdv = sqrt(mx_abs_stdv/(no_of_bins*(no_of_bins-1)));

    double my_abs_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    my_abs_stdv += (mys_abs[l]-my_abs_av)*(mys_abs[l]-my_abs_av);
    my_abs_stdv = sqrt(my_abs_stdv/(no_of_bins*(no_of_bins-1)));

    double mz_abs_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    mz_abs_stdv += (mzs_abs[l]-mz_abs_av)*(mzs_abs[l]-mz_abs_av);
    mz_abs_stdv = sqrt(mz_abs_stdv/(no_of_bins*(no_of_bins-1)));

    // Squared
    double mxsq_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    mxsq_stdv += (mxssq[l]-mxsq_av)*(mxssq[l]-mxsq_av);
    mxsq_stdv = sqrt(mxsq_stdv/(no_of_bins*(no_of_bins-1)));

    double mysq_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    mysq_stdv += (myssq[l]-mysq_av)*(myssq[l]-mysq_av);
    mysq_stdv = sqrt(mysq_stdv/(no_of_bins*(no_of_bins-1)));

    double mzsq_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    mzsq_stdv += (mzssq[l]-mzsq_av)*(mzssq[l]-mzsq_av);
    mzsq_stdv = sqrt(mzsq_stdv/(no_of_bins*(no_of_bins-1)));

    // Quadrupled
    double mxquad_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    mxquad_stdv += (mxsquad[l]-mxquad_av)*(mxsquad[l]-mxquad_av);
    mxquad_stdv = sqrt(mxquad_stdv/(no_of_bins*(no_of_bins-1)));

    double myquad_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    myquad_stdv += (mysquad[l]-myquad_av)*(mysquad[l]-myquad_av);
    myquad_stdv = sqrt(myquad_stdv/(no_of_bins*(no_of_bins-1)));

    double mzquad_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    mzquad_stdv += (mzsquad[l]-mzquad_av)*(mzsquad[l]-mzquad_av);
    mzquad_stdv = sqrt(mzquad_stdv/(no_of_bins*(no_of_bins-1)));


    //----------------The correlation function-------------------//
    //cout << "Should be entering the printing phase for spin correlation function" << endl;
    if(calculatespincorrelationfunction)
    {
        cout << "Going to print spin correlation to file" << endl;
        // Completing the average over the spin correlation function
        for(int k=0; k<N; k++)
        {
            correlation_functionx_av[k]   /= no_of_bins;
            correlation_functiony_av[k]   /= no_of_bins;
            correlation_functionz_av[k]   /= no_of_bins;
            correlation_functiontot_av[k] /= no_of_bins;
        }

        vector<double> spinncorrxstdv   = vector<double>(N);
        vector<double> spinncorrystdv   = vector<double>(N);
        vector<double> spinncorrzstdv   = vector<double>(N);
        vector<double> spinncorrtotstdv = vector<double>(N);
        //energy += (energies[l]-energy_av)*(energies[l]-energy_av);
        // Finding the standard deviation
        double littlex, littley, littlez, littletot;
        for(int l=0; l<no_of_bins; l++)
        {   // For every bin
            vector<double> corrfuncx_thisbin   = correlation_functionx_store[l];
            vector<double> corrfuncy_thisbin   = correlation_functiony_store[l];
            vector<double> corrfuncz_thisbin   = correlation_functionz_store[l];
            vector<double> corrfunctot_thisbin = correlation_functiontot_store[l];
            for(int k=0; k<N; k++)
            {   // For every particle
                littlex   = corrfuncx_thisbin[k]-correlation_functionx_av[k];
                littley   = corrfuncy_thisbin[k]-correlation_functiony_av[k];
                littlez   = corrfuncz_thisbin[k]-correlation_functionz_av[k];
                littletot = corrfunctot_thisbin[k]-correlation_functiontot_av[k];
                spinncorrxstdv[k]   += littlex*littlex;
                spinncorrystdv[k]   += littley*littley;
                spinncorrzstdv[k]   += littlez*littlez;
                spinncorrtotstdv[k] += littletot*littletot;
            }
        }
        for(int k=0; k<N; k++)
        {
            spinncorrxstdv[k]   = sqrt(spinncorrxstdv[k]/(no_of_bins*(no_of_bins-1)));
            spinncorrystdv[k]   = sqrt(spinncorrystdv[k]/(no_of_bins*(no_of_bins-1)));
            spinncorrzstdv[k]   = sqrt(spinncorrzstdv[k]/(no_of_bins*(no_of_bins-1)));
            spinncorrtotstdv[k] = sqrt(spinncorrtotstdv[k]/(no_of_bins*(no_of_bins-1)));
        }

        // Print it to file
        for(int k=0; k<N; k++)
        {
            //cout << "Correlation function, site " << k << ": " << correlation_function_av[k] << "; stdv: " << spinncorrstdv[k] << endl;
            spcorFilex << std::setprecision(std::numeric_limits<double>::digits10 + 1) << correlation_functionx_av[k] << " " << spinncorrxstdv[k] << endl; // Easier to get the desired output this way
            spcorFiley << std::setprecision(std::numeric_limits<double>::digits10 + 1) << correlation_functiony_av[k] << " " << spinncorrystdv[k] << endl; // Easier to get the desired output this way
            spcorFilez << std::setprecision(std::numeric_limits<double>::digits10 + 1) << correlation_functionz_av[k] << " " << spinncorrzstdv[k] << endl; // Easier to get the desired output this way
            spcorFiletot << std::setprecision(std::numeric_limits<double>::digits10 + 1) << correlation_functiontot_av[k] << " " << spinncorrtotstdv[k] << endl; // Easier to get the desired output this way
        }
        // Fourier transformed back:
        cout << "Going to print the s.c.f. f.t.ed back" << endl;
        for(int k=0; k<N; k++)    ftcorrelation_function_av[k] = ftcorrelation_function_av[k]/no_of_bins;

        vector<double> ftspinncorrstdv = vector<double>(N);
        //energy += (energies[l]-energy_av)*(energies[l]-energy_av);
        // Finding the standard deviation
        for(int l=0; l<no_of_bins; l++)
        {   // For every bin
            vector<double> ftcorrfunc_thisbin = ftcorrelation_function_store[l];
            //cout << "ftcorrfunc_thisbin retrieved" << endl;
            for(int k=0; k<N; k++)
            {   // For every particle
                double little = ftcorrfunc_thisbin[k]-ftcorrelation_function_av[k];
                ftspinncorrstdv[k] += little*little;
            }
        }
        for(int k=0; k<N; k++)    ftspinncorrstdv[k] = sqrt(ftspinncorrstdv[k]/(no_of_bins*(no_of_bins-1)));

        // Write to file?
        for(int k=0; k<N; k++)
        {
            //cout << "FT Correlation function, site " << k << ": " << ftcorrelation_function_av[k] << "; ftstdv: " << ftspinncorrstdv[k] << endl;
            ftspcorFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << ftcorrelation_function_av[k] << " " << ftspinncorrstdv[k] << endl; // Easier to get the desired output this way
        }
    }
    if(SCBUG)    cout << "DONE! with writing the correlation functions to file";
    //cout << "Now going to print to allFile" << endl;

    //-----------------------Printing----------------------------//
    allFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << beta << " " << energy_av << " " << E_stdv << " " << energy_sq_av << " " << Esq_stdv << " " << cv << " " << cv_stdv << " " <<  mx_av ;
    allFile << " " << mx_stdv << " " << my_av << " " << my_stdv << " " << mz_av << " " << mz_stdv << " " << ar_av << " " << ar_stdv;
    allFile << " " << mxsq_av << " " << mxsq_stdv << " " << mysq_av << " " << mysq_stdv << " " << mzsq_av << " " << mzsq_stdv;
    allFile << " " << mxquad_av << " " << mxquad_stdv << " " << myquad_av << " " << myquad_stdv << " " << mzquad_av << " " << mzquad_stdv;
    allFile << " " << mx_abs_av << " " << mx_abs_stdv << " " << my_abs_av << " " << my_abs_stdv << " " << mz_abs_av << " " << mz_abs_stdv << endl;

    endtime = clock();
    total_time = (endtime - starttime)/(double) CLOCKS_PER_SEC;
    cout << "Time MC steps and measurements: "  << total_time << endl;
    //cout << "Done with the Monte Carlo procedure this time around" << endl;

    // Destroy plan
    fftw_destroy_plan(px);
    fftw_destroy_plan(py);
    fftw_destroy_plan(pz);
    fftw_destroy_plan(pinv);

    // Print to terminal if desired
    bool printenergytoterminal = false;
    if(printenergytoterminal)    cout << "Energy for beta = " << beta << ": " << energy_av;
    if(printenergytoterminal)    cout << ";  Standard deviation: " << E_stdv << endl;

    bool CONFWR  = true;
    if(CONFWR)
    {
        double sx, sy, sz;
        // Make files to print to
        ofstream configFile;
        char *filename = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename, "%s_spinconfig.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        configFile.open(filename);
        delete filename;

        ofstream neighbourFile;
        char *filename2 = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename2, "%s_neighbours.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        neighbourFile.open(filename2);
        delete filename2;

        int dim = mylattice.dim; // Header with dimension and no of elements in each direction

        configFile << beta << " " << dim << " " << mylattice.L1 << " " << mylattice.L2 << " " << mylattice.L3 << endl;
        vector<double> position;
        vector<int> neighbours;
        int no_of_neighbours = mylattice.no_of_neighbours;
        for(int n=0; n<N; n++)
        {
            sx = mylattice.sites[n].spinx;
            sy = mylattice.sites[n].spiny;
            sz = mylattice.sites[n].spinz;

            if(dim==1)    configFile << n << " " << sx << " " << sy << " " << sz << endl;
            if(dim==2)
            {
                // Printing the config
                position = mylattice.sitepositions[n];
                configFile << position[0] << " " << position[1] << " " << sx << " " << sy << " " << sz << endl;

                // Printing the list of neighbours
                neighbours = mylattice.siteneighbours[n];
                for(int j; j<no_of_neighbours; j++)    neighbourFile << neighbours[j] << " ";
                neighbourFile << endl;
            }
            if(dim==3)
            {
                // Printing the config
                position = mylattice.sitepositions[n];
                configFile << position[0] << " " << position[1] << " " << position[2] << " " << sx << " " << sy << " " << sz << endl;

                // Printing the list of neighbours
                neighbours = mylattice.siteneighbours[n];
                for(int j=0; j<no_of_neighbours; j++)
                {
                    //cout << "In loop, j = " << j << "; Neighbour = " << neighbours[j] << endl;
                    neighbourFile << neighbours[j] << " ";
                }
                neighbourFile << endl;

            } // End if-test over dimensions

        } // End loop over particles
        configFile.close();
        neighbourFile.close();
    } // End if-test CONFWR
}


void MonteCarlo::mcstepf_metropolis(double beta) //, std::default_random_engine generator_u, std::default_random_engine generator_v, std::default_random_engine generator_n, std::default_random_engine generator_prob,  std::uniform_real_distribution<double> distribution_prob, std::uniform_real_distribution<double> distribution_u, std::uniform_real_distribution<double> distribution_v, std::uniform_int_distribution<int> distribution_n)
{   // Include a counter that measures how many 'flips' are accepted. But what to do with it? Write to file?

    bool HUMBUG   = false;  // The humbug is defeated.
    bool MINIBUG  = false;  // For very little output. Couts the probs and/or the energies at each step
    bool CONTRBUG = false;
    if(HUMBUG)    cout << "In mcstepf_metropolis" << endl;
    double changes = 0;
    if(HUMBUG)   cout << "In mcstepf. Looping over spins now" << endl;
    for(int n=0; n<N; n++)
    {
        if(HUMBUG)    cout << "Inside loop in mcstepf. n = " << n << endl;

        // Is this random enough?
        int k = distribution_n(generator_n);   // Draw a random site
        // Or, in case random is inaccurate:
        //(Tested against 2x2x2 fcc. Yields the same results within the error bars.)
        //int k = (int)floor(N*ran2(&seed));
        if(HUMBUG)    cout << "Random site k drawn. k = " << k << endl;

        double sx = mylattice.sites[k].spinx;
        double sy = mylattice.sites[k].spiny;
        double sz = mylattice.sites[k].spinz;

        if(HUMBUG)    cout << "Components of spin " << k << " accessed" << endl;
        //cout << "k :" << k << endl;

        // Changing the spin (tentatively):

        //double u = distribution_u(generator_u);
        //double v = distribution_v(generator_v);

        double u = ran2(&seed1);
        double v = ran2(&seed1);
        if(HUMBUG)    cout << "Have drawn random numbers in mcstepf" << endl;

        ///*
        double theta = acos(1.0-2.0*u);
        double phi = 2.0*M_PI*v;

        double sintheta = sin(theta);
        double sx_t = sintheta*cos(phi);
        double sy_t = sintheta*sin(phi);
        //double sz_t = sqrt(1-sintheta*sintheta); // This produces faulty output...
        double sz_t = cos(theta);
        /*
        double divisor = sqrt(sx_t*sx_t + sy_t*sy_t + sz_t*sz_t);
        sx_t /= divisor;
        sy_t /= divisor;
        sz_t /= divisor;
        */

        //cout << "Spin before trial: [" << sx << "," << sy << "," << sz << "]" << endl;
        //cout << "Spin after trial:  [" << sx_t << "," << sy_t << "," << sz_t << "]" << endl;

        //cout << "sx_t = " << sx_t << endl;
        //cout << "sy_t = " << sy_t << endl;
        //cout << "sz_t = " << sz_t << endl;
        //cout <<  "Checking if spin is normalized: " <<  (sx_t*sx_t+sy_t*sy_t+sz_t*sz_t) << endl;

        //*/

        /*
        double zetasq = 2.0;  // zeta squared. Value to initialize loop.
        double zeta1;
        double zeta2;
        while(zetasq>=1.0) // Is this correct? Could equal 1 too, right?
        {
            double r_1 = distribution_u(generator_u);
            double r_2 = distribution_v(generator_v);
            r_1 = 1.0-2.0*r_1;
            r_2 = 1.0 -2.0*r_2;

            zeta1 = 1.0-2.0*r_1;
            zeta2 = 1.0-2.0*r_2;

            zetasq = zeta1*zeta1 + zeta2*zeta2;
            //if(DEBUG)    cout << "Inside while loop, n = " << n << ": zeta1 = " << zeta1 << "; zeta2 = " << zeta2 << "; zeta^2 = " << zetasq << endl;

            //cout << "Inside loop, n = " << n << " " << zetasq << endl;
        }

        //cout << "Outside loop, n = " << n << " " << zetasq << endl;
        //if(DEBUG)    cout << "Outside while loop, n = " << n << ": zeta1 = " << zeta1 << "; zeta2 = " << zeta2 <<"; zeta^2 = " << zetasq << endl;
        double thesqrt = sqrt(1-zetasq);
        double sx_t = 2.0*zeta1*thesqrt;
        double sy_t = 2.0*zeta2*thesqrt;
        double sz_t = 1.0-2.0*zetasq;

        double stot = sqrt(sx_t*sx_t + sy_t*sy_t + sz_t*sz_t);
        sx_t = sx_t/stot;
        sy_t = sy_t/stot;
        sz_t = sz_t/stot;

        //cout << "sx_t = " << sx_t << endl;
        //cout << "sy_t = " << sy_t << endl;
        //cout << "sz_t = " << sz_t << endl;
        //cout <<  "Checking if spin is normalized: " <<  (sx_t*sx_t+sy_t*sy_t+sz_t*sz_t) << endl;
        // Check normalization
        if(HUMBUG)    cout << "Normalized? S^2 = " << (sx_t*sx_t + sy_t*sy_t + sz_t*sz_t) << endl;
        */

        if(HUMBUG)    cout << "Have made a uniform spherical distribution using them" << endl;
        // Energy contribution after spin change

        double energy_diff = 0; // Resetting the energy difference for every n
        //cout << "energy_diff reset: " << energy_diff << endl;
        if(sianisotropy)
        {
            if(HUMBUG)    cout << "In sianisotropy in mcstepf" << endl;
            if(CONTRBUG)  cout << "In sianisotropy in mcstepf" << endl;
            double Dix = mylattice.sites[k].Dix;
            double Diy = mylattice.sites[k].Diy;
            double Diz = mylattice.sites[k].Diz;
            if(CONTRBUG)  cout << "Dix: " << Dix << "; Diy: " << Diy << "; Diz = " << Diz << endl;
            //cout << "Dix : " << Dix << "; Diy : " << Diy << "; Diz : " << Diz << endl;
            energy_diff += (Dix*(sx_t*sx_t -sx*sx) + Diy*(sy_t*sy_t-sy*sy)+ Diz*(sz_t*sz_t -sz*sz));
            //cout << "Contr from sian: " << (Dix*(sx_t*sx_t -sx*sx) + Diy*(sy_t*sy_t-sy*sy)+ Diz*(sz_t*sz_t -sz*sz)) << endl;
        }
        if(magfield)
        {
            if(HUMBUG)    cout << "In magfield in mcstepf" << endl;
            if(CONTRBUG)  cout << "In magfield in mcstepf" << endl;
            double hx = mylattice.sites[k].hx;
            double hy = mylattice.sites[k].hy;
            double hz = mylattice.sites[k].hz;
            energy_diff += hx*(sx-sx_t) + hy*(sy-sy_t) + hz*(sz-sz_t);
        }
        if(isotropic)
        {
            if(HUMBUG)    cout << "In isotropic in mcstepf" << endl;
            if(CONTRBUG)  cout << "In isotropic in mcstepf" << endl;
            int nneighbours;
            double partnerspinx = 0;
            double partnerspiny = 0;
            double partnerspinz = 0;
            //cout << "Partnerspins before loop: partnerspinx = " << partnerspinx << " " << "; partnerspiny = " << partnerspiny << " " << "; partnerspinz = " << partnerspinz << endl;
            // Determining the number of neighbours
            if(notperiodic)    nneighbours = mylattice.sites[k].no_of_neighbours_site;
            else               nneighbours = no_of_neighbours;
            if(HUMBUG)    cout << "no_of_neighbours, spin " << k << ": " << nneighbours << endl;
            if(CONTRBUG)  cout << "no_of_neighbours, spin " << k << ": " << nneighbours << endl;
            //cout << nneighbours;
            for(int j=0; j<nneighbours; j++)
            {
                // Picking out the neighbour
                int l = mylattice.sites[k].bonds[j].siteindex2;
                if(HUMBUG)    cout << "Spin no. " << l << " chosen." << endl;

                // Picking out the J each time (may vary depending on bond type)
                double J = mylattice.sites[k].bonds[j].J;
                //cout << "Neighbour no. " << j << "; J = " << J << endl;
                if(HUMBUG)    cout << "Neighbour no. " << j << "; J = " << J << endl;
                if(CONTRBUG)
                {
                    string dir = mylattice.sites[k].bonds[j].direction;
                    if(!mylattice.extended)    cout << "Neighbour no. " << j << ", i.e. spin " << l << "; J = " << J << endl;
                    else                       cout << "Neighbour no. " << j << ", i.e. spin " << l << "; J = " << J << "; Direction: " << dir << endl;
                }
                if(J!=0)
                {
                    double sxk = mylattice.sites[l].spinx;  // The neighbours does not change
                    double syk = mylattice.sites[l].spiny;
                    double szk = mylattice.sites[l].spinz;
                    partnerspinx += J*sxk;
                    partnerspiny += J*syk;
                    partnerspinz += J*szk;
                    //cout << "j = " << j << "; partnerspinx = " << partnerspinx << " " << "; partnerspiny = " << partnerspiny << " " << "; partnerspinz = " << partnerspinz << endl;
                    //cout << "J!=0, get a contribution from bond " << j << ": " << J*(sxk*(sx_t-sx)+syk*(sy_t-sy)+szk*(sz_t-sz)) << endl;
                }
            }
            if(HUMBUG)    cout << "Out of that blasted loop!" << endl;
            energy_diff += partnerspinx*(sx_t-sx) + partnerspiny*(sy_t-sy) + partnerspinz*(sz_t-sz);
            if(CONTRBUG)  cout << "Contr from Heisenberg: " <<  partnerspinx*(sx_t-sx) + partnerspiny*(sy_t-sy) + partnerspinz*(sz_t-sz) << endl;
        }
        if(dm)
        {
            if(HUMBUG)    cout << "In dm in mcstepf" << endl;
            if(CONTRBUG)  cout << "In dm in mcstepf" << endl;
            // Determining the number of neighbours
            //double dmcontrib = 0;
            double detsign;
            int nneighbours;
            if(notperiodic)    nneighbours = mylattice.sites[k].no_of_neighbours_site;
            else               nneighbours = no_of_neighbours;
            for(int j=0; j<nneighbours; j++)
            {
                int l = mylattice.sites[k].bonds[j].siteindex2;
                bool increasing = mylattice.sites[k].bonds[j].increasing;
                if(increasing)           detsign = 1.0;
                else                     detsign = -1.0;

                //cout << "Sign: " << detsign << endl;
                if(HUMBUG)    cout << "Spin no. " << l << " chosen." << endl;

                double Dx = mylattice.sites[k].bonds[j].Dx;
                double Dy = mylattice.sites[k].bonds[j].Dy;
                double Dz = mylattice.sites[k].bonds[j].Dz;
                if(HUMBUG)    cout << "Bonds accessed" << endl;

                double sxk = mylattice.sites[l].spinx;
                double syk = mylattice.sites[l].spiny;
                double szk = mylattice.sites[l].spinz;
                if(HUMBUG)    cout << "Components of spin no. " << l << " accessed." << endl;

                //dmcontrib += detsign*(Dx*((sy_t-sy)*szk-syk*(sz_t-sz))+Dy*((sz_t-sz)*sxk-szk*(sx_t-sx))+Dz*((sx_t-sx)*syk-(sy_t-sy)*sxk));
                energy_diff += detsign*(Dx*((sy_t-sy)*szk-syk*(sz_t-sz))+Dy*((sz_t-sz)*sxk-szk*(sx_t-sx))+Dz*((sx_t-sx)*syk-(sy_t-sy)*sxk));
            }

            if(HUMBUG)    cout << "Finding the energy difference from dm" << endl;
            if(HUMBUG)    cout << "Done with dm in mcstepf" << endl;
            //cout << "Contribution from DM: " << dmcontrib << endl;
        }
        //cout << "energy_diff = " << energy_diff << endl;

        double energy_new = energy_old + energy_diff;
        //cout << "Spin difference in each direction:  [" << sx_t-sx << "," << sy_t-sy << "," << sz_t-sz << "]" << endl;
        //cout << "Energy_difference: " << energy_diff << endl;

        // Updating the energy and the state according to Metropolis
        if(energy_new <= energy_old)
        {
            // Updating the spin
            mylattice.sites[k].spinx = sx_t;
            mylattice.sites[k].spiny = sy_t;
            mylattice.sites[k].spinz = sz_t;

            if(MINIBUG) cout << "energy_old = " << energy_old << "; energy_new = " << energy_new << ". ACCEPTED!" << endl;

            // Updating the energy
            energy_old = energy_new;

            // Updating changes to get the acceptance rate
            changes+=1;

            // Misc
            if(HUMBUG)   cout << "ENERGY DECREASED!" << endl;            
            if(MAJORDEBUG)    debug1d2p();
        }
        else
        {
            //double prob = exp(-beta*(energy_new-energy_old));
            double prob = exp(-beta*(energy_diff));
            double drawn = distribution_prob(generator_prob);
            //double drawn = ran2(&seed2);

            if(MINIBUG)   cout << "MCProb: " << prob << "; Number drawn: " << drawn << endl;
            if(HUMBUG)    cout << "Suggesting energy increase. Probability of success: " << prob << "; number drawn: " << drawn << endl;
            if(drawn<prob)
            {
                if(MINIBUG) cout << "drawn<prob. MOVE ACCEPTED! Energy increased: energy_old = " << energy_old << "; energy_diff = " << energy_diff << "; energy_new = " << energy_new << endl;
                // Updating the spin
                mylattice.sites[k].spinx = sx_t;
                mylattice.sites[k].spiny = sy_t;
                mylattice.sites[k].spinz = sz_t;

                // Misc
                //cout << "ENERGY INCREASED!" << endl;
                if(HUMBUG)    cout << "ENERGY INCREASED! energy_old = " << energy_old << "; energy_diff = " << energy_diff << "; energy_new = " << energy_new << endl;
                if(HUMBUG)    cout << "Success" << endl;
                if(MAJORDEBUG)    debug1d2p();

                // Updating the energy
                energy_old = energy_new;

                // Updating changes to get the acceptance rate
                changes+=1;
            }
        }
        //if(n<10)    cout << "Energy, algorithm: " << energy_old << "; energy, hardcoded: " << check_the_energy() << endl;
    } // End loop over n. MCstep done
    acceptancerate = changes/N; // Write the percentage of hits to file.
}

void MonteCarlo::endsims()
{
    //print.closeAllFiles();
    if(allFile.is_open())      allFile.close();
    if(printeveryMCstep)       if(bigFile.is_open())      bigFile.close();
    if(calculatespincorrelationfunction)
    {
        if(spcorFilex.is_open())      spcorFilex.close();
        if(spcorFiley.is_open())      spcorFiley.close();
        if(spcorFilez.is_open())      spcorFilez.close();
        if(spcorFiletot.is_open())    spcorFiletot.close();
        if(ftspcorFile.is_open())     ftspcorFile.close();
    }
}


// Debugging/testing functions
double MonteCarlo::check_the_energy()
{ // Really similar to initialize_energy, but does not overwrite energy_old
    //if(DEBUG)    cout << "In MonteCarlo_check_the_energy" << endl;
    double energy_out = 0;
    double energy_contribution_sites = 0;
    double energy_contribution_bonds = 0;

    bool BEDBUG = false;
    for(int i=0; i<N; i++)
    {
        if(BEDBUG) cout << "In loop to set initial energy" << endl;
        // Spin i
        double sx = mylattice.sites[i].spinx;
        double sy = mylattice.sites[i].spiny;
        double sz = mylattice.sites[i].spinz;
        // Contribution from sites
        if(sianisotropy)
        {
            //cout << "in sianisotropy, init" << endl;
            if(BEDBUG)    cout << "In sianisotropy" << endl;
            double Dix = mylattice.sites[i].Dix;   // Should these lie in Lattice instead of in Lattice::site?
            double Diy = mylattice.sites[i].Diy;
            double Diz = mylattice.sites[i].Diz;
            energy_contribution_sites += (Dix*sx*sx + Diy*sy*sy+ Diz*sz*sz);
            //cout << "Energy contribution from sites: " << energy_contribution_sites << endl;
        }
        if(magfield)
        {
            if(BEDBUG)    cout << "In magfield" << endl;
            double hx = mylattice.sites[i].hx;
            double hy = mylattice.sites[i].hy;
            double hz = mylattice.sites[i].hz;
            energy_contribution_sites -= hx*sx + hy*sy + hz*sz;
        }
        // Contribution from bonds
        if(isotropic)
        {
            if(BEDBUG)   cout << "In isotropic" << endl;
            double nneighbours;
            double partnerspinx = 0;
            double partnerspiny = 0;
            double partnerspinz = 0;
            // Determining the number of neighbours of the site
            if(notperiodic)    nneighbours = mylattice.sites[i].no_of_neighbours_site;
            else               nneighbours = no_of_neighbours;
            //cout << "Looping over " << nneighbours << " neighbours" << endl;
            for(int j=0; j<nneighbours; j++)
            {
                if(BEDBUG)    cout << "in loop in isotropic, j = " << j << endl;
                int k = mylattice.sites[i].bonds[j].siteindex2;
                if(i<k)    // To avoid double counting. More computationally efficient than halving the energy.
                {
                    if(BEDBUG)    cout << "Have set k in loop, k = " << k << "; N = " << N << endl;
                    double J = mylattice.sites[i].bonds[j].J;
                    if(BEDBUG)    cout << "Have accessed J in bond between spins" << endl;
                    double sxk = mylattice.sites[k].spinx;
                    double syk = mylattice.sites[k].spiny;
                    double szk = mylattice.sites[k].spinz;
                    if(BEDBUG)    cout << "Have accessed the components of the spin on the other end" << endl;
                    partnerspinx += J*sxk;
                    partnerspiny += J*syk;
                    partnerspinz += J*szk;
                    if(BEDBUG)    cout << "Have gathered this contribution into partnerspin" << endl;
                }
            }
            if(BEDBUG)   cout << "Done with the loop in isotropic" << endl;
            if(BEDBUG)   cout << "partnerspinx = " << partnerspinx << "; partnerspiny = " << partnerspiny << "; partnerspinz = " << partnerspinz << endl;
            energy_contribution_bonds += (partnerspinx*sx + partnerspiny*sy + partnerspinz*sz);
            //cout << "Energy contribution from bonds: " << energy_contribution_bonds << endl;
            if(BEDBUG)     cout << "Done with isotropic" << endl;
        }
        if(dm)
        {
            if(BEDBUG)    cout << "In dm" << endl;
            // Double loops and stuff. Could maybe make this more efficient
            int nneighbours;
            // Determining the number of neighbours of the site
            if(notperiodic)    nneighbours = mylattice.sites[i].no_of_neighbours_site;
            else               nneighbours = no_of_neighbours;
            for(int j=0; j<nneighbours; j++)
            {
                int k = mylattice.sites[i].bonds[j].siteindex2; // Hope I can actually get to this value.
                if(i<k)    // To avoid double counting. Hopefully saves time for large systems
                {
                    double Dx = mylattice.sites[i].bonds[j].Dx;
                    double Dy = mylattice.sites[i].bonds[j].Dy;
                    double Dz = mylattice.sites[i].bonds[j].Dz;

                    double sxk = mylattice.sites[k].spinx;
                    double syk = mylattice.sites[k].spiny;
                    double szk = mylattice.sites[k].spinz;

                    energy_contribution_bonds += Dx*(sy*szk-syk*sz)+Dy*(sz*sxk-szk*sx)+Dz*(sx*syk-sy*sxk);
                }
            }
        }
        if(BEDBUG) cout << "Done with one, onto the others" << endl;
    }
    energy_out = energy_contribution_sites + energy_contribution_bonds;
    return energy_out;
}





void MonteCarlo::debug1d2p()
{
    // Only to work for chain with two particles, homogenous J
    double spin0x = mylattice.sites[0].spinx;
    double spin0y = mylattice.sites[0].spiny;
    double spin0z = mylattice.sites[0].spinz;
    double spin1x = mylattice.sites[1].spinx;
    double spin1y = mylattice.sites[1].spiny;
    double spin1z = mylattice.sites[1].spinz;
    double homJ   = mylattice.sites[0].bonds[0].J;

    double test_energy = 2*homJ*(spin0x*spin1x+spin0y*spin1y+spin0z*spin1z);

    compareFile  << test_energy << " " << energy_old << endl;
}

void MonteCarlo::testFFTW()
{
    vector<double> spins_in_z = vector<double>(N);
    double spinthing = 0;
    for(int i=0;i<N;i++)
    {
        spinthing += 0.1;
        spins_in_z[i] = spinthing;
        //spins_in_z[i] = 1.0; // or 1.0
        //spins_in_z[i] = pow(-1.0,i); // or 1.0
    }

    // To check special cases, we want to set the spins manually
    /*
    spins_in_z[0] = 1.0;
    spins_in_z[1] = 1.0;
    spins_in_z[2] = -1.0;
    spins_in_z[3] = -1.0;
    */

    // Giving qlimit
    int dim = mylattice.dim;
    cout << "Dimension of lattice = " << dim << endl;
    vector<int> dimlen = mylattice.dimlengths;
    int qlimit = 1; // To be multiplied;
    for(int l=0; l<(dim-1); l++)
    {    // Looping over all dimensions but the last
        qlimit *= dimlen[l];
    }
    qlimit *= (mylattice.dimlengths[dim-1]/2)+1;

    //cout << "qlimit: " << qlimit << endl;

    vector< complex<double> > qconf(N);  // Output array
    vector<double> correlation_function_av = vector<double >(qlimit); // Change this?
    //double time_start = clock();
    givezplanforFFT(spins_in_z, qconf);
    fftw_execute(pz);

    //double time_end = clock();
    //double time_realtocomplex = (time_end-time_start)/CLOCKS_PER_SEC;

    cout << "Printing the spin correlation function" << endl;
    for(int l=0; l<qlimit; l++)
    {   // Accumulating the average
        // Should consider whether I actually want an output array of half the length
        // of the input array.
        cout << "index l: " << l << endl;
        correlation_function_av[l] = (qconf[l]*conj(qconf[l])).real()/(N*N); //Should I divide by N? Or sqrt(N)?
        // Multiplying a complex number by its complex conjugate should yield a real number.
        // but I call .real() to get the right data type.
        cout << "Value of correlation function: " <<  correlation_function_av[l] << " " << endl;
    }

    cout << endl; // End line to get ready for new result

    for(int i=0; i<qlimit; i++)    cout << qconf[i] << endl;

    // Patching it:
    vector<double> correlation_function = vector<double>(N);
    if(mylattice.dim==1) // Chain
    {
        int L = mylattice.L;
        cout << "L = " << L << endl;
        for(int n=0; n<N; n++) // Loop over the q-values
        {   // FFTW procedure
            int index;
            if(n<=(int)N/2)    index = n;
            else               index = N-n;
            correlation_function[n] = (qconf[index]*conj(qconf[index])).real()/(N*N);
            cout << "By FFTW: n = " << n << "; correlation function = " << correlation_function[n] << endl;
        }
    }
    if(mylattice.dim==3) // Simple cubic or fcc
    {
        int L1 = mylattice.L1;
        int L2 = mylattice.L2;
        int L3 = mylattice.L3;
        //for(int k=0; k<N; k++)    correlation_function_av_bin[k] = 0;
        int elinar = 0; // For retrieving the elements residing in the output array
        for(int n=0; n<N; n++)
        {   // FFTW procedure
            vector<int> cord = mylattice.sitecoordinates[n];
            int n1 = cord[0];
            int n2 = cord[1];
            int n3 = cord[2];
            int index;
            if(n3<=(int)L3/2)
            {   // If our element is stored in the output array, we retrieve it
                index = elinar;
                elinar++;        // We move one index forward in the output array
            }
            else
            {   // If our element is not stored in the output array, we retrieve its
                // complex conjugate. We needn't do anything with it as the result is a
                // complex number times its complex conjugate
                index = L2*((int)L3/2+1)*(L1-n1)+(int)L3/2*(L2-n2)+L3-n3;
                //cout << "Retrieving the index of the cc" << endl;
            } // End if-tests
            correlation_function[n] = (qconf[index]*conj(qconf[index])).real()/(N*N);
            cout << "By FFTW: n = " << n << "; correlation function = " << correlation_function[n] << endl;
        }
    }


    //cout << "Time spent on FFT: " << time_realtocomplex << " s" << endl;

    // Doing it manually
    /*
    if(mylattice.dim==1)
    {
        cout << "Testing my 'manual' implementation of the correlation function: " << endl;
        double q, spincorr_manual_rel, spincorr_manual_cpl, relfac, cplfac, spinzk, spinzl;
        int L = mylattice.L;
        for(int n=0; n<N; n++)
        {
            spincorr_manual_rel = 0;
            spincorr_manual_cpl = 0;
            for(int k=0; k<N; k++) // Loop over r'
            {
                for(int l=0; l<N; l++) // Loop over r
                {
                    q = 2*M_PI*(k-l)*n/L;
                    //cout << "q = " << q << endl;
                    relfac = cos(q);
                    cplfac = sin(q);
                    spinzk = spins_in_z[k];
                    spinzl = spins_in_z[l];
                    spincorr_manual_rel += relfac*spinzk*spinzl;
                    spincorr_manual_cpl += cplfac*spinzk*spinzl;
                }

            }
            spincorr_manual_rel = spincorr_manual_rel/(N*N);
            spincorr_manual_cpl = spincorr_manual_cpl/(N*N);
            cout << "FFTW: n = " << n << "; correlation function = " << correlation_function[n] << endl;
            cout << "'Manually': n = " << n << "; correlation function = ( " << spincorr_manual_rel << " , " << spincorr_manual_cpl << " )" << endl;
            cout << "Difference, FFTW-manually (real part): " << correlation_function[n]-spincorr_manual_rel << endl << endl;
        } // End loop over n
    } // End if-test dim==1
    if(mylattice.dim==3)
    {
        double q, spincorr_manual_rel, spincorr_manual_cpl, relfac, cplfac, spinzk, spinzl;
        int L1 = mylattice.L1;
        int L2 = mylattice.L2;
        int L3 = mylattice.L3;
        cout << "L1 = " << L1 << "; L2 = " << L2 << "; L3 = " << L3 << endl;

        // Manual procedure
        for(int n=0; n<N; n++)
        {
            spincorr_manual_rel = 0;
            spincorr_manual_cpl = 0;
            vector<int> coord = mylattice.sitecoordinates[n];
            int n1 = coord[0];
            int n2 = coord[1];
            int n3 = coord[2];
            //cout << "n1 = " << n1 << "; n2 = " << n2 << "; n3 = " << n3 << endl << endl;
            for(int k=0; k<N; k++) // Loop over r'
            {
                vector<int> cord = mylattice.sitecoordinates[k];
                int k1 = cord[0]; // Number in the a1-direction
                int k2 = cord[1]; // Number in the a2-direction
                int k3 = cord[2]; // Number in the a3-direction
                //cout << "k1 = " << k1 << "; k2 = " << k2 << "; k3 = " << k3 << endl;
                for(int l=0; l<N; l++) // Loop over r
                {   // Something similar for cubic/fcc, but need to get indices.
                    vector<int> cord = mylattice.sitecoordinates[l];
                    int l1 = cord[0]; // Number in the a1-direction
                    int l2 = cord[1]; // Number in the a2-direction
                    int l3 = cord[2]; // Number in the a3-direction
                    //cout << "l1 = " << l1 << "; l2 = " << l2 << "; l3 = " << l3 << endl;

                    double q1 = 2*M_PI*(k1-l1)*n1/L1;
                    double q2 = 2*M_PI*(k2-l2)*n2/L2;
                    double q3 = 2*M_PI*(k3-l3)*n3/L3;
                    //cout << "k1-l1 = " <<  k1-l2 << endl;
                    q =  q1 + q2 + q3;
                    //cout << "q1 = " << q1 << "; q2 = " << q2 << "; q3 = " << q3 << endl;
                    //cout << "q = " << q << endl;
                    relfac = cos(q);
                    cplfac = sin(q);
                    spinzk = spins_in_z[k];
                    spinzl = spins_in_z[l];
                    spincorr_manual_rel += relfac*spinzk*spinzl;
                    spincorr_manual_cpl += cplfac*spinzk*spinzl;
                } // End loop over l (r)
            } // End loop over k (r')
            spincorr_manual_rel = spincorr_manual_rel/(N*N);
            spincorr_manual_cpl = spincorr_manual_cpl/(N*N);

            cout << "FFTW: n = " << n << "; correlation function = " << correlation_function[n] << endl;
            cout << "'Manually': n = " << n << "; correlation function = ( " << spincorr_manual_rel << " , " << spincorr_manual_cpl << " )" << endl;
            cout << "Difference, FFTW-manually (real part): " << correlation_function[n]-spincorr_manual_rel << endl << endl;
        } // End loop over n
    } // End loop over if-test dim==3
    */
}


void MonteCarlo::compareFFTW_withmanual(double beta)
{
    // For FFTW:
    // Making the plan
    std::vector<double> spins_in_z = std::vector<double>(N);
    vector< complex<double> > qconf(N);  // Output array
    givezplanforFFT(spins_in_z, qconf);

    // Making an array for the correlation function, one for each value of q
    vector<double> correlation_function        = vector<double>(N);
    for(int i=0; i<N; i++)    correlation_function[i] = 0;

    // For manual calculation:
    double q = 0;
    double q1 = 0;
    double q2 = 0;
    double q3 = 0;
    double relfac = 0;
    double cplfac = 0;
    double spincorr_manual_rel = 0;
    double spincorr_manual_cpl = 0;
    vector<double> spincorr_manual_rel_array(N);

    // Equilibration steps
    for(int i=0; i<eqsteps; i++)        mcstepf_metropolis(beta);

    for(int i=0; i<mcsteps_inbin; i++) // A small number of Monte Carlo steps
    {
        // Running the Monte Carlo procedure. One MCstep
        mcstepf_metropolis(beta);

        // Update the spin array
        //for(int j=0; j<N; j++)    spins_in_z[j] = mylattice.sites[j].spinz;
        for(int j=0; j<N; j++)
        {
            spins_in_z[j] = mylattice.sites[j].spinz;
            //cout << "Spin " << j << ": " << spins_in_z[j] << endl;
        }
        // Then execute the plan
        fftw_execute(pz);

        if(i==mcsteps_inbin-1)
        {// Print the result after some MCsteps, to see if something in the routine goes wrong
            if(mylattice.dim==1) // Chain
            {
                int L = mylattice.L;
                cout << "L = " << L << endl;
                for(int n=0; n<N; n++) // Loop over the q-values
                {   // FFTW procedure
                    // Resetting quantities
                    spincorr_manual_rel = 0;
                    spincorr_manual_cpl = 0;

                    int index;
                    if(n<=(int)N/2)    index = n;
                    else               index = N-n;
                    correlation_function[n] = (qconf[index]*conj(qconf[index])).real()/(N*N);
                    // Manual procedure
                    for(int k=0; k<N; k++) // Loop over r'
                    {
                        for(int l=0; l<N; l++) // Loop over r
                        {
                            q = 2*M_PI*(k-l)*n/L;
                            relfac = cos(q);
                            cplfac = sin(q);
                            //cout << "n = " << n << "; k = " << k << "; l =" << l << "; q = " << q << endl;
                            double spinzk = mylattice.sites[k].spinz;
                            double spinzl = mylattice.sites[l].spinz;
                            spincorr_manual_rel += relfac*spinzk*spinzl;
                            spincorr_manual_cpl += cplfac*spinzk*spinzl;
                        }
                    }
                    spincorr_manual_rel = spincorr_manual_rel/(N*N);
                    spincorr_manual_cpl = spincorr_manual_cpl/(N*N);
                    spincorr_manual_rel_array[n] = spincorr_manual_rel;

                    cout << "By FFTW: n = " << n << "; correlation function = " << correlation_function[n] << endl;
                    cout << "'Manually': n = " << n << "; correlation function = ( " << spincorr_manual_rel << " , " << spincorr_manual_cpl << " )" << endl;
                    cout << "Difference between FFTW and manual (real part): " << correlation_function[n]-spincorr_manual_rel << endl << endl;
                }
            }
            if(mylattice.dim==3) // Simple cubic or fcc
            {
                int L1 = mylattice.L1;
                int L2 = mylattice.L2;
                int L3 = mylattice.L3;
                //for(int k=0; k<N; k++)    correlation_function_av_bin[k] = 0;
                int elinar = 0; // For retrieving the elements residing in the output array
                for(int n=0; n<N; n++)
                {   // FFTW procedure
                    vector<int> cord = mylattice.sitecoordinates[n];
                    int n1 = cord[0];
                    int n2 = cord[1];
                    int n3 = cord[2];
                    int index;
                    if(n3<=(int)L3/2)
                    {   // If our element is stored in the output array, we retrieve it
                        index = elinar;
                        elinar++;        // We move one index forward in the output array
                    }
                    else
                    {   // If our element is not stored in the output array, we retrieve its
                        // complex conjugate. We needn't do anything with it as the result is a
                        // complex number times its complex conjugate
                        index = L2*((int)L3/2+1)*((L1-n1)%L1)+((int)L3/2+1)*((L2-n2)%L2)+((L3-n3)%L3);
                        //cout << "Retrieving the index of the cc" << endl;
                    } // End if-tests
                    correlation_function[n] = (qconf[index]*conj(qconf[index])).real()/(N*N);

                    // Manual procedure
                    // Resetting quantities
                    spincorr_manual_rel = 0;
                    spincorr_manual_cpl = 0;
                    for(int k=0; k<N; k++) // Loop over r'
                    {
                        vector<int> cord = mylattice.sitecoordinates[k];
                        int k1 = cord[0]; // Number in the a1-direction
                        int k2 = cord[1]; // Number in the a2-direction
                        int k3 = cord[2]; // Number in the a3-direction
                        for(int l=0; l<N; l++) // Loop over r
                        {   // Something similar for cubic/fcc, but need to get indices.
                            vector<int> cord = mylattice.sitecoordinates[l];
                            int l1 = cord[0]; // Number in the a1-direction
                            int l2 = cord[1]; // Number in the a2-direction
                            int l3 = cord[2]; // Number in the a3-direction

                            q1 = 2*M_PI*(k1-l1)*n1/L1;
                            q2 = 2*M_PI*(k2-l2)*n2/L2;
                            q3 = 2*M_PI*(k3-l3)*n3/L3;
                            q =  q1 + q2 + q3;
                            //cout << "q = " << q << endl;
                            relfac = cos(q);
                            cplfac = sin(q);
                            // Standard, should work:
                            double spinzk = mylattice.sites[k].spinz;
                            double spinzl = mylattice.sites[l].spinz;
                            // Safety route (doesn't work either):
                            //double spinzk = spins_in_z[k];
                            //double spinzl = spins_in_z[l];
                            spincorr_manual_rel += relfac*spinzk*spinzl;
                            spincorr_manual_cpl += cplfac*spinzk*spinzl;
                        }
                    }
                    spincorr_manual_rel = spincorr_manual_rel/(N*N);
                    spincorr_manual_cpl = spincorr_manual_cpl/(N*N);
                    spincorr_manual_rel_array[n] = spincorr_manual_rel;

                    /*
                    cout << "n = " << n << endl;
                    cout << "n1 = " << n1 << "; n2 = " << n2 << "; n3 = " << n3 << endl;
                    cout << "By FFTW: n = " << n << "; correlation function = " << correlation_function[n] << endl;
                    cout << "'Manually': n = " << n << "; correlation function = ( " << spincorr_manual_rel << " , " << spincorr_manual_cpl << " )" << endl;
                    cout << "Difference, FFTW-manually (real part): " << correlation_function[n]-spincorr_manual_rel << endl << endl;
                    */
                } // End loop over n

            } // End if-test dimension
        } // End if-test last step
    } // End Monte Carlo steps

    // For cleaner diagnostics output:
    cout << "DIAGNOSTICS" << endl;
    vector<int> cord;
    int n1, n2, n3;
    //int index;
    int errorcounter = 0;
    int signwarnings = 0;
    int L1 = mylattice.L1;
    int L2 = mylattice.L2;
    int L3 = mylattice.L3;

    for(int n=0; n<N; n++)
    {
        cord = mylattice.sitecoordinates[n];
        n1   = cord[0];
        n2   = cord[1];
        n3   = cord[2];
        if(abs(spincorr_manual_rel_array[n]-correlation_function[n])>1e-16)
        {

            cout << "n = " << n << endl;
            cout << "n1 = " << n1 << "; n2 = " << n2 << "; n3 = " << n3 << endl;
            cout << "By FFTW: n = " << n << "; correlation function = " << correlation_function[n] << endl;
            cout << "'Manually': n = " << n << "; correlation function = ( " << spincorr_manual_rel_array[n] << " , " << spincorr_manual_cpl << " )" << endl;
            cout << "Difference, FFTW-manually (real part): " << correlation_function[n]-spincorr_manual_rel_array[n] << endl << endl;
            errorcounter++;
        }
        if(L1-n1<0)
        {
            cout << "Warning! Wrong sign! n = " << n << ", L1-n1 = " << L1-n1 << endl;
            signwarnings++;
        }
        if(L2-n2<0)
        {
            cout << "Warning! Wrong sign! n = " << n << ", L2-n2 = " << L2-n2 << endl;
            signwarnings++;
        }
        if(L3-n3<0)
        {
            cout << "Warning! Wrong sign! n = " << n << ", L3-n3 = " << L3-n3 << endl;
            signwarnings++;
        }
    }

    if(!(errorcounter+signwarnings)==0)
    {
        cout << "Number of significant deviations in result: " << errorcounter << endl;
        cout << "Number of warnings about sign (Li-ni): " << signwarnings << endl;
        cout << "L1 = " << L1 << "; L2 = " << L2 << "; L3 = " << L3 << endl;
        cout << "(int)L3/2+1 = " << (int)L3/2+1 << endl;
    }
    else    cout << "No significant deviations. Program OK." << endl;


    /*
    if(mylattice.L==6) // This test is only for the LxLxL fcc lattice.
    { // Should have taken the average, though...
        vector<int> yline6(6);
        yline6[0] = 0; yline6[1] = 47; yline6[2] = 88; yline6[3] = 129; yline6[4] = 170; yline6[5] = 211;

        cout << "The yline:" << endl;
        cout << "n_y = " << 0 << endl;
        cout << "FFTW: " << correlation_function[0] << endl;
        cout << "Hard coded: " << spincorr_manual_rel_array[0] << endl;

        cout << "n_y = " << 1 << endl;
        cout << "FFTW: " << correlation_function[yline6[1]] << endl;
        cout << "Hard coded: " << spincorr_manual_rel_array[yline6[1]] << endl;

        cout << "n_y = " << 2 << endl;
        cout << "FFTW: " << correlation_function[yline6[2]] << endl;
        cout << "Hard coded: " << spincorr_manual_rel_array[yline6[2]] << endl;

        cout << "n_y = " << 3 << endl;
        cout << "FFTW: " << correlation_function[yline6[3]] << endl;
        cout << "Hard coded: " << spincorr_manual_rel_array[yline6[3]] << endl;

        cout << "n_y = " << 4 << endl;
        cout << "FFTW: " << correlation_function[yline6[4]] << endl;
        cout << "Hard coded: " << spincorr_manual_rel_array[yline6[4]] << endl;

        cout << "n_y = " << 5 << endl;
        cout << "FFTW: " << correlation_function[yline6[5]] << endl;
        cout << "Hard coded: " << spincorr_manual_rel_array[yline6[5]] << endl;
    }
    */
}

void MonteCarlo::compareFFTW_withmanual_av(double beta)
{
    // For FFTW:
    // Making the plan
    std::vector<double> spins_in_z = std::vector<double>(N);
    vector< complex<double> > qconf(N);  // Output array
    givezplanforFFT(spins_in_z, qconf);

    // Making an array for the correlation function, one for each value of q
    vector<double> correlation_function        = vector<double>(N);
    for(int i=0; i<N; i++)    correlation_function[i] = 0;
    vector<double> correlation_function_av     = vector<double>(N);
    for(int i=0; i<N; i++)    correlation_function_av[i] = 0;

    // For manual calculation:
    double q = 0;
    double q1 = 0;
    double q2 = 0;
    double q3 = 0;
    double relfac = 0;
    double cplfac = 0;
    double spincorr_manual_rel = 0;
    double spincorr_manual_cpl = 0;
    vector<double> spincorr_manual_rel_array(N);
    vector<double> spincorr_manual_rel_array_av(N);

    // Equilibration steps
    for(int i=0; i<eqsteps; i++)        mcstepf_metropolis(beta);

    for(int i=0; i<mcsteps_inbin; i++) // A small number of Monte Carlo steps
    {
        cout << "MCstep " << i << endl;
        // Running the Monte Carlo procedure. One MCstep
        mcstepf_metropolis(beta);

        // Update the spin array
        for(int j=0; j<N; j++)    spins_in_z[j] = mylattice.sites[j].spinz;

        // Then execute the plan
        fftw_execute(pz);
        // For different system types
        if(mylattice.dim==1) // Chain
        {
            int L = mylattice.L;
            cout << "L = " << L << endl;
            for(int n=0; n<N; n++) // Loop over the q-values
            {
                // Resetting quantities
                spincorr_manual_rel = 0;
                spincorr_manual_cpl = 0;

                // FFTW procedure
                int index;
                if(n<=(int)N/2)    index = n;
                else               index = N-n;
                correlation_function_av[n] += (qconf[index]*conj(qconf[index])).real()/(N*N);

                // Manual procedure
                for(int k=0; k<N; k++) // Loop over r'
                {
                    for(int l=0; l<N; l++) // Loop over r
                    {
                        q = 2*M_PI*(k-l)*n/L;
                        relfac = cos(q);
                        cplfac = sin(q);
                        double spinzk = mylattice.sites[k].spinz;
                        double spinzl = mylattice.sites[l].spinz;
                        spincorr_manual_rel += relfac*spinzk*spinzl;
                        spincorr_manual_cpl += cplfac*spinzk*spinzl;
                    }

                }
                    spincorr_manual_rel = spincorr_manual_rel/(N*N);
                    spincorr_manual_cpl = spincorr_manual_cpl/(N*N);
                    spincorr_manual_rel_array[n]    = spincorr_manual_rel;
                    spincorr_manual_rel_array_av[n] += spincorr_manual_rel;
                }
            }
        if(mylattice.dim==2)
        {
            int L1 = mylattice.L1;
            int L2 = mylattice.L2;
            int elinar = 0; // For retrieving the elements residing in the output array
            for(int n=0; n<N; n++)
            {
                vector<int> cord = mylattice.sitecoordinates[n];
                int n1 = cord[0];
                int n2 = cord[1];
                int index;
                if(n2<=(int)L2/2)
                {   // If our element is stored in the output array, we retrieve it
                    index = elinar;
                    elinar++;        // We move one index forward in the output array
                    //cout << "I have met ms elinar" << endl;
                }
                else
                {   // If our element is not stored in the output array, we retrieve its
                    // complex conjugate. We needn't do anything with it as the result is a
                    // complex number times its complex conjugate
                    index = ((int)L2/2+1)*((L1-n1)%L1)+(L2-n2)%L2;
                    //cout << "Retrieving the index of the cc" << endl;
                } // End if-tests
                correlation_function_av[n] += (qconf[index]*conj(qconf[index])).real()/(N*N);

                // Manual procedure
                // Resetting quantities
                spincorr_manual_rel = 0;
                spincorr_manual_cpl = 0;
                for(int k=0; k<N; k++) // Loop over r'
                {
                    vector<int> cord = mylattice.sitecoordinates[k];
                    int k1 = cord[0]; // Number in the x-direction
                    int k2 = cord[1]; // Number in the y-direction
                    for(int l=0; l<N; l++) // Loop over r
                    {   // Something similar for cubic/fcc, but need to get indices.
                        vector<int> cord = mylattice.sitecoordinates[l];
                        int l1 = cord[0]; // Number in the x-direction
                        int l2 = cord[1]; // Number in the y-direction

                        q1 = 2*M_PI*(k1-l1)*n1/L1;
                        q2 = 2*M_PI*(k2-l2)*n2/L2;

                        q =  q1 + q2;
                        //cout << "q = " << q << endl;
                        relfac = cos(q);
                        cplfac = sin(q);
                        // Standard, should work:
                        double spinzk = mylattice.sites[k].spinz;
                        double spinzl = mylattice.sites[l].spinz;
                        // Safety route (doesn't work either):
                        //double spinzk = spins_in_z[k];
                        //double spinzl = spins_in_z[l];
                        spincorr_manual_rel += relfac*spinzk*spinzl;
                        spincorr_manual_cpl += cplfac*spinzk*spinzl;
                    }
                }
                spincorr_manual_rel = spincorr_manual_rel/(N*N);
                spincorr_manual_cpl = spincorr_manual_cpl/(N*N);
                spincorr_manual_rel_array[n]     = spincorr_manual_rel;
                spincorr_manual_rel_array_av[n] += spincorr_manual_rel;
            }
        }
        if(mylattice.dim==3) // Simple cubic or fcc
        {
            int L1 = mylattice.L1;
            int L2 = mylattice.L2;
            int L3 = mylattice.L3;
            //for(int k=0; k<N; k++)    correlation_function_av_bin[k] = 0;
            int elinar = 0; // For retrieving the elements residing in the output array
            for(int n=0; n<N; n++)
            {   // FFTW procedure
                vector<int> cord = mylattice.sitecoordinates[n];
                int n1 = cord[0];
                int n2 = cord[1];
                int n3 = cord[2];
                int index;
                if(n3<=(int)L3/2)
                {   // If our element is stored in the output array, we retrieve it
                    index = elinar;
                    elinar++;        // We move one index forward in the output array
                }
                else
                {   // If our element is not stored in the output array, we retrieve its
                    // complex conjugate. We needn't do anything with it as the result is a
                    // complex number times its complex conjugate
                    index = L2*((int)L3/2+1)*((L1-n1)%L1)+((int)L3/2+1)*((L2-n2)%L2)+((L3-n3)%L3);
                    //cout << "Retrieving the index of the cc" << endl;
                } // End if-tests
                correlation_function_av[n] += (qconf[index]*conj(qconf[index])).real()/(N*N);

                // Manual procedure
                // Resetting quantities
                spincorr_manual_rel = 0;
                spincorr_manual_cpl = 0;
                for(int k=0; k<N; k++) // Loop over r'
                {
                    vector<int> cord = mylattice.sitecoordinates[k];
                    int k1 = cord[0]; // Number in the a1-direction
                    int k2 = cord[1]; // Number in the a2-direction
                    int k3 = cord[2]; // Number in the a3-direction
                    for(int l=0; l<N; l++) // Loop over r
                    {   // Something similar for cubic/fcc, but need to get indices.
                        vector<int> cord = mylattice.sitecoordinates[l];
                        int l1 = cord[0]; // Number in the a1-direction
                        int l2 = cord[1]; // Number in the a2-direction
                        int l3 = cord[2]; // Number in the a3-direction

                        q1 = 2*M_PI*(k1-l1)*n1/L1;
                        q2 = 2*M_PI*(k2-l2)*n2/L2;
                        q3 = 2*M_PI*(k3-l3)*n3/L3;
                        q =  q1 + q2 + q3;
                        //cout << "q = " << q << endl;
                        relfac = cos(q);
                        cplfac = sin(q);
                        // Standard, should work:
                        double spinzk = mylattice.sites[k].spinz;
                        double spinzl = mylattice.sites[l].spinz;
                        // Safety route (doesn't work either):
                        //double spinzk = spins_in_z[k];
                        //double spinzl = spins_in_z[l];
                        spincorr_manual_rel += relfac*spinzk*spinzl;
                        spincorr_manual_cpl += cplfac*spinzk*spinzl;
                    }
                }
                spincorr_manual_rel = spincorr_manual_rel/(N*N);
                spincorr_manual_cpl = spincorr_manual_cpl/(N*N);
                spincorr_manual_rel_array[n]     = spincorr_manual_rel;
                spincorr_manual_rel_array_av[n] += spincorr_manual_rel;
            } // End loop over n
        } // End if-test dimension
    } // End Monte Carlo steps

    for(int n=0; n<N; n++)
    {
        spincorr_manual_rel_array_av[n] /= mcsteps_inbin;
        correlation_function_av[n] /= mcsteps_inbin;
    }
    // For cleaner diagnostics output:
    cout << "DIAGNOSTICS" << endl;
    vector<int> cord;
    int n1, n2, n3;
    //int index;
    int errorcounter = 0;
    int signwarnings = 0;
    int L1 = mylattice.L1;
    int L2 = mylattice.L2;
    int L3 = mylattice.L3;

    for(int n=0; n<N; n++)
    {
        cord = mylattice.sitecoordinates[n];
        n1   = cord[0];
        n2   = cord[1];
        n3   = cord[2];
        if(abs(spincorr_manual_rel_array_av[n]-correlation_function_av[n])>1e-16)
        {

            cout << "n = " << n << endl;
            cout << "n1 = " << n1 << "; n2 = " << n2 << "; n3 = " << n3 << endl;
            cout << "By FFTW: n = " << n << "; correlation function = " << correlation_function_av[n] << endl;
            cout << "'Manually': n = " << n << "; correlation function = ( " << spincorr_manual_rel_array_av[n] << " , " << spincorr_manual_cpl << " )" << endl;
            cout << "Difference, FFTW-manually (real part): " << correlation_function_av[n]-spincorr_manual_rel_array_av[n] << endl << endl;
            errorcounter++;
        }
        if(L1-n1<0)
        {
            cout << "Warning! Wrong sign! n = " << n << ", L1-n1 = " << L1-n1 << endl;
            signwarnings++;
        }
        if(L2-n2<0)
        {
            cout << "Warning! Wrong sign! n = " << n << ", L2-n2 = " << L2-n2 << endl;
            signwarnings++;
        }
        if(L3-n3<0)
        {
            cout << "Warning! Wrong sign! n = " << n << ", L3-n3 = " << L3-n3 << endl;
            signwarnings++;
        }
    }

    if(!(errorcounter+signwarnings)==0)
    {
        cout << "Number of significant deviations in result: " << errorcounter << endl;
        cout << "Number of warnings about sign (Li-ni): " << signwarnings << endl;
        cout << "L1 = " << L1 << "; L2 = " << L2 << "; L3 = " << L3 << endl;
        cout << "(int)L3/2+1 = " << (int)L3/2+1 << endl;
    }
    else    cout << "No significant deviations. Program OK." << endl;
}

void MonteCarlo::shortsim(double beta)
{
    // For FFTW:
    // Making the plan
    std::vector<double> spins_in_z = std::vector<double>(N);
    vector< complex<double> > qconf(N);  // Output array
    givezplanforFFT(spins_in_z, qconf);

    // Making an array for the correlation function, one for each value of q
    vector<double> correlation_function        = vector<double>(N);
    for(int i=0; i<N; i++)    correlation_function[i] = 0;
    vector<double> correlation_function_av     = vector<double>(N);
    for(int i=0; i<N; i++)    correlation_function_av[i] = 0;

    // Equilibration steps
    for(int i=0; i<eqsteps; i++)        mcstepf_metropolis(beta);

    for(int i=0; i<mcsteps_inbin; i++) // A small number of Monte Carlo steps
    {
        cout << "MCstep " << i << endl;
        // Running the Monte Carlo procedure. One MCstep
        mcstepf_metropolis(beta);

        // Update the spin array
        for(int j=0; j<N; j++)    spins_in_z[j] = mylattice.sites[j].spinz;

        // Then execute the plan
        fftw_execute(pz);
        // For different system types
        if(mylattice.dim==1) // Chain
        {
            int L = mylattice.L;
            cout << "L = " << L << endl;
            for(int n=0; n<N; n++) // Loop over the q-values
            {
                // FFTW procedure
                int index;
                if(n<=(int)N/2)    index = n;
                else               index = N-n;
                correlation_function_av[n] += (qconf[index]*conj(qconf[index])).real()/(N*N);
            }
        }
        if(mylattice.dim==2)
        {
            int L1 = mylattice.L1;
            int L2 = mylattice.L2;
            int elinar = 0; // For retrieving the elements residing in the output array
            for(int n=0; n<N; n++)
            {
                vector<int> cord = mylattice.sitecoordinates[n];
                int n1 = cord[0];
                int n2 = cord[1];
                int index;
                if(n2<=(int)L2/2)
                {   // If our element is stored in the output array, we retrieve it
                    index = elinar;
                    elinar++;        // We move one index forward in the output array
                    //cout << "I have met ms elinar" << endl;
                }
                else
                {   // If our element is not stored in the output array, we retrieve its
                    // complex conjugate. We needn't do anything with it as the result is a
                    // complex number times its complex conjugate
                    index = ((int)L2/2+1)*((L1-n1)%L1)+(L2-n2)%L2;
                    //cout << "Retrieving the index of the cc" << endl;
                } // End if-tests
                correlation_function_av[n] += (qconf[index]*conj(qconf[index])).real()/(N*N);
            }
        }
        if(mylattice.dim==3) // Simple cubic or fcc
        {
            int L1 = mylattice.L1;
            int L2 = mylattice.L2;
            int L3 = mylattice.L3;
            //for(int k=0; k<N; k++)    correlation_function_av_bin[k] = 0;
            int elinar = 0; // For retrieving the elements residing in the output array
            for(int n=0; n<N; n++)
            {   // FFTW procedure
                vector<int> cord = mylattice.sitecoordinates[n];
                int n1 = cord[0];
                int n2 = cord[1];
                int n3 = cord[2];
                int index;
                if(n3<=(int)L3/2)
                {   // If our element is stored in the output array, we retrieve it
                    index = elinar;
                    elinar++;        // We move one index forward in the output array
                }
                else
                {   // If our element is not stored in the output array, we retrieve its
                    // complex conjugate. We needn't do anything with it as the result is a
                    // complex number times its complex conjugate
                    index = L2*((int)L3/2+1)*((L1-n1)%L1)+((int)L3/2+1)*((L2-n2)%L2)+((L3-n3)%L3);
                    //cout << "Retrieving the index of the cc" << endl;
                } // End if-tests
                correlation_function_av[n] += (qconf[index]*conj(qconf[index])).real()/(N*N);
            } // End loop over n
        } // End if-test dimension
    } // End Monte Carlo steps

    for(int n=0; n<N; n++)        correlation_function_av[n] /= mcsteps_inbin;
}

void MonteCarlo::test_couplings_strengths()
{
    // Seeing that our switches work
    if(isotropic)       cout << "bool isotropic true. Heisenberg terms will be considered." << endl;
    if(sianisotropy)    cout << "bool sianisotropy true" << endl;
    if(magfield)        cout << "bool magfield true" << endl;
    if(dm)              cout << "bool dm true" << endl;

    cout << "We are printing the couplings at work here:" << endl;
    int n = 1; // We could just test for all, but that muddles our overview. Do this first, then, maybe, throw the following section into a loop
    int nneighbours;
    if(notperiodic)    nneighbours = mylattice.sites[n].no_of_neighbours_site;
    else               nneighbours = no_of_neighbours;
    if(isotropic)
    {
        cout << "Heisenberg couplings:" << endl;
        for(int j=0; j<nneighbours; j++)
        {
            double J = mylattice.sites[n].bonds[j].J;
            if(!mylattice.extended)            cout << "Neighbour no.: " << j << ", J = " << J << endl;
            else
            {
                string dir = mylattice.sites[n].bonds[j].direction;
                cout << "Neighbour no.: " << j << ", J = " << J << "; Direction: " << dir << endl;
            }
        }
    }
    if(sianisotropy)
    {
        cout << "Single ion anisotropy strengths:" << endl;
        double Dix = mylattice.sites[n].Dix;   // Should these lie in Lattice instead of in Lattice::site?
        double Diy = mylattice.sites[n].Diy;
        double Diz = mylattice.sites[n].Diz;
        cout << "Dix = " << Dix << "; Diy = " << Diy << "; Diz = " << Diz << endl;
    }
    if(magfield)
    {
        cout << "Magnetic field strength:" << endl;
        double hx = mylattice.sites[n].hx;   // Should these lie in Lattice instead of in Lattice::site?
        double hy = mylattice.sites[n].hy;
        double hz = mylattice.sites[n].hz;
        cout << "hx = " << hx << "; hy = " << hy << "; hz = " << hz << endl;
    }
    if(dm)
    {
        cout << "Dzyaloshinskii-Moriya couplings:" << endl;
        double Dx = mylattice.sites[n].bonds[0].Dx; // The same for all bonds, at leat in my case.
        double Dy = mylattice.sites[n].bonds[0].Dy; // Could probably store it in Lattice instead of Bonds
        double Dz = mylattice.sites[n].bonds[0].Dz;
        cout << "Dx = " << Dx << "; Dy = " << Dy << "; Dz = " << Dz << endl;
    }

}

/*
void MonteCarlo::compareFFTW_withmanual_super(double beta)
{
    for(int i=0; i<mcsteps_inbin; i++)
    {
        int L1 = mylattice.L1;
        int L2 = mylattice.L2;
        int L3 = mylattice.L3;
        int elinar = 0; // For retrieving the elements residing in the output array
        for(int n=0; n<N; n++)
        {
            vector<int> cord = mylattice.sitecoordinates[n];
            int n1 = cord[0];
            int n2 = cord[1];
            int n3 = cord[2];
            int index;
            if(n3<=(int)L3/2)
            {   // If our element is stored in the output array, we retrieve it
                index = elinar;
                elinar++;        // We move one index forward in the output array
                //cout << "I have met ms elinar" << endl;
            }
            else
            {   // If our element is not stored in the output array, we retrieve its
                // complex conjugate. We needn't do anything with it as the result is a
                // complex number times its complex conjugate
                index = L2*((int)L3/2+1)*((L1-n1)%L1)+((int)L3/2+1)*((L2-n2)%L2)+((L3-n3)%L3);
                //cout << "Retrieving the index of the cc" << endl;
            } // End if-tests
            if(SCBUG)    cout << "Adding to the spin correlation function" << endl;
            correlation_function_av_bin[n] += (qconf[index]*conj(qconf[index])).real()/(N*N);
    }
}
*/




