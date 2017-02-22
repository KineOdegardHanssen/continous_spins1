#include "montecarlo.h"

MonteCarlo::MonteCarlo()
{
}

MonteCarlo::MonteCarlo(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
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
    mylattice = Lattice(L,isotropic, sianisotropy, magfield, dm);
    mylattice.setstrengths(sitestrengthsin, heisenbergin, dm_in);
    cout << "Instance of class Lattice initialized" << endl;

    // Type of Lattice
    if(periodic)
    {
        if(type_lattice=='F')      mylattice.fcc_helical_initialize();          // F for fcc
        if(type_lattice=='E')      mylattice.fcc_helical_initialize_extended(); // E for extended
        else if(type_lattice=='C') mylattice.cubic_helical_initialize();        // C for cubic
        else if(type_lattice=='Q') mylattice.quadratic_helical_initialize();    // Q for quadratic
        else if(type_lattice=='O') mylattice.chain_periodic_initialize();       // O for one-dimensional
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
        // Setting up file to print to
        char *filename = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename, "%s_spincorrelationfunction.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        spcorFile.open(filename);
        delete filename;

        cout << "File set" << endl;

        if(mylattice.dim>1)  // The chain is trivial
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

            if(mylattice.dim==2)
            {
                cout << "For quadratic, in if-test" << endl;
                for(int i=0; i<N; i++) // Possibly only up to N/2.
                {
                    vector<int> ns = mylattice.sitecoordinates[i];
                    cout << "ns retrieved" << endl;
                    // These must be changed if we change into possibly setting L1, L2, L3 different
                    // Don't really need b1, b2, b3, could just use a1, a2, a3 multiplied by 2*M_PI...
                    // Could be more general
                    qx = ns[0]*mylattice.b1[0]/L1; // Not so sure abour this, just need it to compile
                    qy = ns[1]*mylattice.b2[1]/L2; // Double check
                    //qz = ns[2]*mylattice.b3[2]/L; // No qz
                    cout << "qvec set" << endl;
                    // Print to file. Site number, qx, qy, qz.
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
                    qx = (ns[0]*mylattice.b1[0] + ns[2]*mylattice.b3[0])/L1; //Not so sure abour this, just need it to compile
                    qy = (ns[0]*mylattice.b1[1] + ns[1]*mylattice.b2[1])/L2; // Double check later
                    qz = (ns[1]*mylattice.b2[2] + ns[2]*mylattice.b3[2])/L3;
                    // Print to file. Site number, qx, qy, qz.
                    qFile << i << " " << qx << " " << qy << " " << qz << endl;
                }
                //cout << "Done printing to qFile" << endl;
            }
        }
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
        if(type_lattice=='F')      mylattice.fcc_helical_initialize();          // F for fcc
        if(type_lattice=='E')      mylattice.fcc_helical_initialize_extended(); // E for extended
        else if(type_lattice=='C') mylattice.cubic_helical_initialize();        // C for cubic
        else if(type_lattice=='Q') mylattice.quadratic_helical_initialize();    // Q for quadratic
        else if(type_lattice=='O') mylattice.chain_periodic_initialize();       // O for one-dimensional
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
        // Setting up file to print to
        char *filename = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename, "%s_spincorrelationfunction.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        spcorFile.open(filename);
        delete filename;

        cout << "File set" << endl;

        if(mylattice.dim>1)  // The chain is trivial
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

            if(mylattice.dim==2)
            {
                cout << "For quadratic, in if-test" << endl;
                for(int i=0; i<N; i++) // Possibly only up to N/2.
                {
                    vector<int> ns = mylattice.sitecoordinates[i];
                    cout << "ns retrieved" << endl;
                    // These must be changed if we change into possibly setting L1, L2, L3 different
                    // Don't really need b1, b2, b3, could just use a1, a2, a3 multiplied by 2*M_PI...
                    // Could be more general
                    qx = ns[0]*mylattice.b1[0]/L1; // Not so sure abour this, just need it to compile
                    qy = ns[1]*mylattice.b2[1]/L2; // Double check
                    //qz = ns[2]*mylattice.b3[2]/L; // No qz
                    cout << "qvec set" << endl;
                    // Print to file. Site number, qx, qy, qz.
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
                    qx = (ns[0]*mylattice.b1[0] + ns[2]*mylattice.b3[0])/L1; //Not so sure abour this, just need it to compile
                    qy = (ns[0]*mylattice.b1[1] + ns[1]*mylattice.b2[1])/L2; // Double check later
                    qz = (ns[1]*mylattice.b2[2] + ns[2]*mylattice.b3[2])/L3;
                    // Print to file. Site number, qx, qy, qz.
                    qFile << i << " " << qx << " " << qy << " " << qz << endl;
                }
                //cout << "Done printing to qFile" << endl;
            }
        }
    }
    if(randomtest)
    {
        char *filename = new char[1000];                                     // File name can have max 1000 characters
        sprintf(filename, "%s_testofrandom.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
        randomtestFile.open(filename);
        delete filename;
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

                    energy_contribution_bonds -= Dx*(sy*szk-syk*sz)+Dy*(sz*sxk-szk*sx)+Dz*(sx*syk-sy*sxk);
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

void MonteCarlo::giveplanforFFT(vector<double>& r, vector<complex<double> >& q)  // Return p? Or have p as a class variable?
{
    int rank = mylattice.dim;               // Dimension of lattice
    vector<int> Ls = mylattice.dimlengths;  // List containing [L], [L1,L2], [L1,L2,L3],
                                            // depending on the lattice
    // p declared as a class variable
    p = fftw_plan_dft_r2c(rank,
                          &Ls[0],
                          &r[0],
                          reinterpret_cast<fftw_complex*>(&q[0]),
                          FFTW_ESTIMATE);
}

void MonteCarlo::runmetropolis(double beta)
{
    if(DEBUG)    cout << "In runmetropolis in MonteCarlo" << endl;
    bool HUMBUG  = false;
    bool LADYBUG = false;

    // Header for spcorFile
    if(calculatespincorrelationfunction)    spcorFile << beta << " " << N << endl;

    // Initializing the energy
    initialize_energy();

    if(DEBUG)    cout << "Done with initialize_energy()" << endl << "Starting equilibration steps" << endl;

    // Equilibration steps
    double starttime = clock();
    for(int i=0; i<eqsteps; i++)
    {
        if(HUMBUG)    cout << "In equilibration steps loop, i = " << i << endl;
        mcstepf_metropolis(beta); //, generator_u, generator_v, generator_n, generator_prob, distribution_prob, distribution_u, distribution_v, distribution_n);
        //if(i<11)      cout << "i = " << i << ", energy: " << energy_old << endl;
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
    std::vector<double> spins_in_z = std::vector<double>(N);
    // Array for the results
    // Determining the length of the array
    int dim = mylattice.dim;
    int qlimit = 1; // To be multiplied;
    for(int l=0; l<(dim-1); l++)
    {    // Looping over all dimensions but the last
        qlimit *= mylattice.dimlengths[l];
    }
    qlimit *= mylattice.dimlengths[dim-1]/2+1;
    vector<double> correlation_function_av = vector<double >(qlimit); // Double check
    for(int l= 0; l<qlimit; l++)    correlation_function_av[l] = 0;

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
    {   // For each bin
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
        if(LADYBUG)
        {
            if(i>0)
            {
                cout << "i = " << i << "; Average energy before loop: " << energy_av/(mcsteps_inbin*i) << endl;
                cout << "i = " << i << "; energies[i]: " << energies[i] << endl;
            }
        }

        // Setting vectors

        for(int j=0; j<mcsteps_inbin; j++)    // Loop over mcsteps in bin
        {   // For each mcstep
            mcstepf_metropolis(beta);
            if(randomtest)
            {
                randomtestFile << ran2(&testseed) << endl;
            }

            // acceptancerate
            acceptancerates[i] += acceptancerate;
            ar_av += acceptancerate;
            // energy
            energies[i]    += energy_old;    // Storing to get the standard deviation
            energies_sq[i] += energy_old*energy_old;
            energy_av      += energy_old;
            energy_sq_av   += energy_old*energy_old;
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
            for(int k=0; k<N; k++)
            {
                mx+= mylattice.sites[k].spinx;
                my+= mylattice.sites[k].spiny;
                mz+= mylattice.sites[k].spinz;
                // For the correlation function
                if(calculatespincorrelationfunction)    spins_in_z[k] = mylattice.sites[k].spinz;
            }
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
                vector< complex<double> > qconf(N);  // Output array
                giveplanforFFT(spins_in_z, qconf);
                fftw_execute(p);

                for(int l=0; l<qlimit; l++)
                {   // Accumulating the average
                    // Should consider whether I actually want an output array of half the length
                    // of the input array.
                    correlation_function_av[l] =  correlation_function_av[l] + (qconf[l]*conj(qconf[l])).real()/(N*N);
                    // Multiplying a complex number by its complex conjugate should yield a real number.
                    // but I call .real() to get the right data type.
                }
            }

            //Print to bigFile
            if(printeveryMCstep)
            {
                bigFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << beta << " " << energy_old << " " << energy_old*energy_old << " " << mx << " " << my << " " << mz << endl;
            }

            // Some sort of measurement of the magnetization... How to do this when we have a continuous spin?
        }  // End loop over mcsteps
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

        if(calculatespincorrelationfunction)
        {   // Take the average and print to file
            // Make the averages (is this more efficient than just calculating the average?)
            for(int l = 0; l<qlimit; l++)
            {
                correlation_function_av[l] = correlation_function_av[l]/(mcsteps_inbin);
                spcorFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << correlation_function_av[l] << " ";  // Should I include a beta, just in case?
            }
            spcorFile << endl; // End line to get ready for new result
        }

    }  // End loops over bins
    //----------------------// Acceptance rate //-----------------------//
    ar_av = ar_av/(mcsteps_inbin*no_of_bins);

    // Standard deviation //
    double ar_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    ar_stdv +=(acceptancerates[l]-ar_av)*(acceptancerates[l]-ar_av);
    ar_stdv = sqrt(ar_stdv/(no_of_bins*(no_of_bins-1)));

    //----------------------// Energy //-----------------------//
    energy_av = energy_av/(mcsteps_inbin*no_of_bins);
    energy_sq_av = energy_sq_av/(mcsteps_inbin*no_of_bins);

    // Error in the energy //
    double E_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    E_stdv += (energies[l]-energy_av)*(energies[l]-energy_av);
    E_stdv = sqrt(E_stdv/(no_of_bins*(no_of_bins-1)));

    double Esq_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    Esq_stdv += (energies_sq[l]-energy_sq_av)*(energies_sq[l]-energy_sq_av);
    Esq_stdv = sqrt(Esq_stdv/(no_of_bins*(no_of_bins-1)));

    //---------------------//Heat capacity//---------------------//
    double cv = beta*beta*(energy_sq_av-energy_av*energy_av);

    // Approximate error in the heat capacity:
    cv_average = cv_average/(mcsteps_inbin*no_of_bins);
    //cout << "cv_average: " << cv_average << endl;

    double cv_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    cv_stdv += (cvs[l]-cv_average)*(cvs[l]-cv_average);
    cv_stdv = sqrt(cv_stdv/(no_of_bins*(no_of_bins-1)));

    //-----------------------//Magnetization//----------------------//
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
}


void MonteCarlo::mcstepf_metropolis(double beta) //, std::default_random_engine generator_u, std::default_random_engine generator_v, std::default_random_engine generator_n, std::default_random_engine generator_prob,  std::uniform_real_distribution<double> distribution_prob, std::uniform_real_distribution<double> distribution_u, std::uniform_real_distribution<double> distribution_v, std::uniform_int_distribution<int> distribution_n)
{   // Include a counter that measures how many 'flips' are accepted. But what to do with it? Write to file?

    bool HUMBUG = false;  // The humbug is defeated.
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
        if(sianisotropy)
        {
            if(HUMBUG)    cout << "In sianisotropy in mcstepf" << endl;
            double Dix = mylattice.sites[k].Dix;
            double Diy = mylattice.sites[k].Diy;
            double Diz = mylattice.sites[k].Diz;
            //cout << "Dix : " << Dix << "; Diy : " << Diy << "; Diz : " << Diz << endl;
            energy_diff += (Dix*(sx_t*sx_t -sx*sx) + Diy*(sy_t*sy_t-sy*sy)+ Diz*(sz_t*sz_t -sz*sz));
            //cout << "Contr from sian: " << (Dix*(sx_t*sx_t -sx*sx) + Diy*(sy_t*sy_t-sy*sy)+ Diz*(sz_t*sz_t -sz*sz)) << endl;
        }
        if(magfield)
        {
            if(HUMBUG)    cout << "In magfield in mcstepf" << endl;
            double hx = mylattice.sites[k].hx;
            double hy = mylattice.sites[k].hy;
            double hz = mylattice.sites[k].hz;
            energy_diff += hx*(sx-sx_t) + hy*(sy-sy_t) + hz*(sz-sz_t);
        }
        if(isotropic)
        {
            if(HUMBUG)    cout << "In isotropic in mcstepf" << endl;
            if(HUMBUG)    cout << "no_of_neighbours = " << no_of_neighbours << endl;
            int nneighbours;
            double partnerspinx = 0;
            double partnerspiny = 0;
            double partnerspinz = 0;
            //cout << "Partnerspins before loop: partnerspinx = " << partnerspinx << " " << "; partnerspiny = " << partnerspiny << " " << "; partnerspinz = " << partnerspinz << endl;
            // Determining the number of neighbours
            if(notperiodic)    nneighbours = mylattice.sites[k].no_of_neighbours_site;
            else               nneighbours = no_of_neighbours;
            //cout << nneighbours;
            for(int j=0; j<nneighbours; j++)
            {
                // Picking out the neighbour
                int l = mylattice.sites[k].bonds[j].siteindex2;
                if(HUMBUG)    cout << "Spin no. " << l << " chosen." << endl;

                // Picking out the J each time (may vary depending on bond type)
                double J = mylattice.sites[k].bonds[j].J;
                //cout << "j = " << j << "; J = " << J << endl;
                double sxk = mylattice.sites[l].spinx;
                double syk = mylattice.sites[l].spiny;
                double szk = mylattice.sites[l].spinz;
                partnerspinx += J*sxk;
                partnerspiny += J*syk;
                partnerspinz += J*szk;
                //cout << "j = " << j << "; partnerspinx = " << partnerspinx << " " << "; partnerspiny = " << partnerspiny << " " << "; partnerspinz = " << partnerspinz << endl;
            }
            if(HUMBUG)    cout << "Out of that blasted loop!" << endl;
            energy_diff += partnerspinx*(sx_t-sx) + partnerspiny*(sy_t-sy) + partnerspinz*(sz_t-sz);
            //cout << "Contr from Heisenberg: " <<  partnerspinx*(sx_t-sx) + partnerspiny*(sy_t-sy) + partnerspinz*(sz_t-sz) << endl;
        }
        if(dm)
        {
            if(HUMBUG)    cout << "In dm in mcstepf" << endl;
            // Determining the number of neighbours
            int nneighbours;
            if(notperiodic)    nneighbours = mylattice.sites[k].no_of_neighbours_site;
            else               nneighbours = no_of_neighbours;
            for(int j=0; j<nneighbours; j++)
            {
                int l = mylattice.sites[k].bonds[j].siteindex2;
                if(HUMBUG)    cout << "Spin no. " << l << " chosen." << endl;

                double Dx = mylattice.sites[k].bonds[j].Dx;
                double Dy = mylattice.sites[k].bonds[j].Dy;
                double Dz = mylattice.sites[k].bonds[j].Dz;
                if(HUMBUG)    cout << "Bonds accessed" << endl;

                double sxk = mylattice.sites[l].spinx;
                double syk = mylattice.sites[l].spiny;
                double szk = mylattice.sites[l].spinz;
                if(HUMBUG)    cout << "Components of spin no. " << l << " accessed." << endl;

                energy_diff += Dx*((sy-sy_t)*szk-syk*(sz-sz_t))+Dy*((sz-sz_t)*sxk-szk*(sx-sx_t))+Dz*((sx-sx_t)*syk-(sy-sy_t)*sxk);
            }

            if(HUMBUG)    cout << "Finding the energy difference from dm" << endl;

        }
        if(HUMBUG)    cout << "Done with dm in mcstepf" << endl;

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

            // Updating the energy
            energy_old = energy_new;

            //cout << "New energy less than old energy. MOVE ACCEPTED!" << endl;

            // Updating changes to get the acceptance rate
            changes+=1;

            // Misc
            if(HUMBUG)   cout << "ENERGY DECREASED!" << endl;            
            if(MAJORDEBUG)    debug1d2p();
        }
        else
        {
            double prob = exp(-beta*(energy_new-energy_old));
            double drawn = distribution_prob(generator_prob);
            //double drawn = ran2(&seed2);

            //cout << "MCProb: " << prob << "; Number drawn: " << drawn << endl;
            if(HUMBUG)    cout << "Suggesting energy increase. Probability of success: " << prob << "; number drawn: " << drawn << endl;
            if(drawn<prob)
            {
                //cout << "drawn<prob. MOVE ACCEPTED!" << endl;
                // Updating the spin
                mylattice.sites[k].spinx = sx_t;
                mylattice.sites[k].spiny = sy_t;
                mylattice.sites[k].spinz = sz_t;

                // Misc
                if(HUMBUG)    cout << "ENERGY INCREASED! energy_old = " << energy_old << "; energy_diff = " << energy_diff << "; energy_new = " << energy_new << endl;
                if(HUMBUG)    cout << "Success" << endl;
                if(MAJORDEBUG)    debug1d2p();

                // Updating the energy
                energy_old = energy_new;

                // Updating changes to get the acceptance rate
                changes+=1;
            }
        }
    }
    acceptancerate = changes/N; // Write the percentage of hits to file.
}

void MonteCarlo::endsims()
{
    //print.closeAllFiles();
    if(allFile.is_open())      allFile.close();
    if(printeveryMCstep)       if(bigFile.is_open())      bigFile.close();
}


// Debugging/testing functions
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
        //spinthing += 0.1;
        //spins_in_z[i] = spinthing;
        //spins_in_z[i] = 1.0; // or 1.0
        spins_in_z[i] = pow(-1.0,i); // or 1.0
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
    giveplanforFFT(spins_in_z, qconf);
    fftw_execute(p);

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

    //cout << "Time spent on FFT: " << time_realtocomplex << " s" << endl;
    /*
    cout << "Complex variant" << endl;
    vector< complex<double> > spins_in_z_complex;
    vector< complex<double> > qconfc(N);  // Output array
    vector<double> correlation_function_avc = vector<double >(N);

    for(int i=0;i<N;i++)
    {
        cout << "i = " << i << endl;
        complex<double> spinz(1.0,0.0);
        spins_in_z_complex[i] = spinz;
    }
    cout << "Hey, girl" << endl;
    for(int i=0;i<N;i++)    cout << spins_in_z_complex[i] << " " << endl;


    time_start = clock();
    int rank = mylattice.dim;               // Dimension of lattice
    vector<int> Ls = mylattice.dimlengths;  // List containing [L], [L1,L2], [L1,L2,L3],
                                            // depending on the lattice
    // p declared as a class variable
    p = fftw_plan_dft_r2c(rank,
                          &Ls[0],
                          reinterpret_cast<fftw_complex*>(&spins_in_z_complex[0]),
                          reinterpret_cast<fftw_complex*>(&qconfc[0]),
                          FFTW_ESTIMATE);

    fftw_execute(p);
    time_end = clock();
    double time_complex = (time_end-time_start)/CLOCKS_PER_SEC;

    cout << "Printing the spin correlation function" << endl;
    for(int l=0; l<N; l++)
    {   // Accumulating the average
        // Should consider whether I actually want an output array of half the length
        // of the input array.
        correlation_function_av[l] = (qconfc[l]*conj(qconfc[l])).real()/(N*N); //Should I divide by N? Or sqrt(N)?
        // Multiplying a complex number by its complex conjugate should yield a real number.
        // but I call .real() to get the right data type.
        cout <<  correlation_function_avc[l] << " ";
    }

    cout << "Time, real to complex Fourier transform: " << time_realtocomplex << endl;
    cout << "Time, complex to complex Fourier transform: " << time_complex << endl;
    */

}
