#include "montecarlo.h"

MonteCarlo::MonteCarlo()
{
}

MonteCarlo::MonteCarlo(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, bool isotropic, bool sianisotropy, bool magfield, bool dm, char type_lattice, string filenamePrefix)
{
    // Initializing class for printing to file
    //print(filenamePrefix);
    //print(filenamePrefix);
    //print.givePrefix(filenamePrefix);
    // For now: Just print to all files. Have other options, of course.
    //print.open_all();

    // Handling runningints (wouldn't want to vary this in one class instance, I guess.)
    this->eqsteps = eqsteps;
    this->mcsteps_inbin = mcsteps_inbin;
    this->no_of_bins = no_of_bins;

    // Setting the bools
    this->isotropic = isotropic;
    this->sianisotropy = sianisotropy;
    this->magfield = magfield;
    this->dm = dm;

    // Initializing some other quantities
    acceptancerate = 0;
    DEBUG = false;

    // Setting up files to print to. Might want to allow more flexibility here
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_everybeta.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    allFile.open(filename);
    delete filename;

    char *filenameb = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenameb, "%s_everyMCstep.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    bigFile.open(filenameb);
    delete filenameb;
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


void MonteCarlo::latticetype(int L, char type_lattice)
{
    if(DEBUG)    cout << "In latticetype in class MonteCarlo. L =" << L << endl;
    // Making the Lattice
    double starttime = clock();
    mylattice = Lattice(L, isotropic, sianisotropy, magfield, dm);
    cout << "Instance of class Lattice initialized" << endl;

    // Type of Lattice
    if(type_lattice=='F')      mylattice.fcc_helical_initialize();
    else if(type_lattice=='C') mylattice.cubic_helical_initialize();
    else if(type_lattice=='Q') mylattice.quadratic_helical_initialize();
    double endtime = clock();
    double total_time = (endtime - starttime)/(double) CLOCKS_PER_SEC;
    if(DEBUG)    cout << "Lattice set, time: " << total_time << endl;

    // The lattice  // Could just extract the bools from here...
    this->N = mylattice.N;
    this->no_of_neighbours = mylattice.no_of_neighbours;
    if(DEBUG)    cout << "Number of sites and neighbours retrieved to MonteCarlo." << endl;
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
        // Contribution from sites
        if(sianisotropy)
        {
            if(BEDBUG)    cout << "In sianisotropy" << endl;
            double Dix = mylattice.sites[i].Dix;
            double Diy = mylattice.sites[i].Diy;
            double Diz = mylattice.sites[i].Diz;
            double sx = mylattice.sites[i].spinx;
            double sy = mylattice.sites[i].spiny;
            double sz = mylattice.sites[i].spinz;
            energy_contribution_sites -= Dix*sx*sx + Diy*sy*sy+ Diz*sz*sz;
            //energy_contribution_sites += sianisotropy_energy(i, mylattice);
        }
        if(magfield)
        {
            if(BEDBUG)    cout << "In magfield" << endl;
            double hx = mylattice.sites[i].hx;
            double hy = mylattice.sites[i].hy;
            double hz = mylattice.sites[i].hz;
            double sx = mylattice.sites[i].spinx;
            double sy = mylattice.sites[i].spiny;
            double sz = mylattice.sites[i].spinz;
            energy_contribution_sites -= hx*sx + hy*sy + hz*sz;
        }
        // Contribution from bonds
        if(isotropic)
        {
            if(BEDBUG)   cout << "In isotropic" << endl;
            // Declare no_of_neighbours here in case
            double partnerspinx = 0;
            double partnerspiny = 0;
            double partnerspinz = 0;
            for(int j=0; j<no_of_neighbours; j++)
            {
                if(BEDBUG)    cout << "in loop in isotropic, j = " << j << endl;
                int k = mylattice.sites[i].bonds[j].siteindex2; // Hope I can actually get to this value.
                //                                          // I could, alternatively, just store the index
                //                                          // of bonds. But then I need to organize the bonds
                //                                          // and that is such a hassle.
                if(BEDBUG)    cout << "Have set k in loop, k = " << k << "; N = " << N << endl;
                double J = mylattice.sites[i].bonds[j].J;
                if(BEDBUG)    cout << "Have accessed J in bond between spins" << endl;
                double sx = mylattice.sites[k].spinx;
                double sy = mylattice.sites[k].spiny;
                double sz = mylattice.sites[k].spinz;
                if(BEDBUG)    cout << "Have accessed the components of the spin on the other end" << endl;
                partnerspinx += J*sx;
                partnerspiny += J*sy;
                partnerspinz += J*sz;
                if(BEDBUG)    cout << "Have gathered this contribution into partnerspin" << endl;
            }
            if(BEDBUG)   cout << "Done with the loop in isotropic" << endl;
            double sx = mylattice.sites[i].spinx;
            double sy = mylattice.sites[i].spiny;
            double sz = mylattice.sites[i].spinz;
            energy_contribution_bonds += partnerspinx*sx + partnerspiny*sy + partnerspinz*sz;
            // half this thing. Or find a reasonable way to not double count.
            if(BEDBUG)     cout << "Done with isotropic" << endl;
        }
        if(dm)
        {
            if(BEDBUG)    cout << "In dm" << endl;
            // Double loops and stuff. Could maybe make this more efficient
            double sxi = mylattice.sites[i].spinx;
            double syi = mylattice.sites[i].spiny;
            double szi = mylattice.sites[i].spinz;
            for(int j=0; j<no_of_neighbours; j++)
            {

                int k = mylattice.sites[i].bonds[j].siteindex2; // Hope I can actually get to this value.
                double Dx = mylattice.sites[i].bonds[j].Dx;
                double Dy = mylattice.sites[i].bonds[j].Dy;
                double Dz = mylattice.sites[i].bonds[j].Dz;

                double sxk = mylattice.sites[k].spinx;
                double syk = mylattice.sites[k].spiny;
                double szk = mylattice.sites[k].spinz;

                energy_contribution_bonds -= Dx*(syi*szk-syk*szi)+Dy*(szi*sxk-szk*sxi)+Dz*(sxi*syk-syi*sxk);
            }
        }
        if(BEDBUG) cout << "Done with one, onto the others" << endl;
    }
    energy_old = energy_contribution_sites + energy_contribution_bonds/2.0;
    if(DEBUG)    cout << "Energy initialized, energy = " << energy_old << endl;
}

void MonteCarlo::reset_energy()
{  // In case it gets stuck in some region, I guess...
    double originalvalue = 1/sqrt(3);
    for(int i=0; i<N; i++)
    {   // Reset all spins
        mylattice.sites[i].spinx = originalvalue;
        mylattice.sites[i].spiny = originalvalue;
        mylattice.sites[i].spinz = originalvalue;
    }

    initialize_energy();
}

void MonteCarlo::runmetropolis(double beta)
{
    if(DEBUG)    cout << "In runmetropolis in MonteCarlo" << endl;
    bool HUMBUG  = false;
    bool LADYBUG = false;

    // Random generators
    std::default_random_engine generator_u;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_u(-1,1);

    std::default_random_engine generator_v;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_v(-1,1);

    std::default_random_engine generator_prob;                    // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_prob(0,1);

    // For index. This is given helical boundary conditions, then I only need one index
    std::default_random_engine generator_n;
    std::uniform_int_distribution<int> distribution_n(0,N-1);

    // Initializing the energy
    initialize_energy();

    if(DEBUG)    cout << "Done with initialize_energy()" << endl << "Starting equilibration steps" << endl;

    // Equilibration steps
    double starttime = clock();
    for(int i=0; i<eqsteps; i++)
    {
        if(HUMBUG)    cout << "In equilibration steps loop, i = " << i << endl;
        mcstepf_metropolis(beta, generator_u, generator_v, generator_n, generator_prob, distribution_prob, distribution_u, distribution_v, distribution_n);
        if(i<11)      cout << "i = " << i << ", energy: " << energy_old << endl;
    }
    double endtime = clock();
    double total_time = (endtime - starttime)/(double) CLOCKS_PER_SEC;
    cout << "Time equilibration steps: " << total_time << endl;


    if(DEBUG)    cout << "Done with equilibration steps" << endl << "Starting the Monte Carlo steps and measurements" << endl;

    // Monte Carlo steps and measurements
    starttime = clock();
    std::vector<double> acceptancerates    = std::vector<double>(no_of_bins);
    std::vector<double> energies    = std::vector<double>(no_of_bins);
    std::vector<double> energies_sq = std::vector<double>(no_of_bins);
    std::vector<double> cvs         = std::vector<double>(no_of_bins);
    std::vector<double> mxs = std::vector<double>(no_of_bins);
    std::vector<double> mys = std::vector<double>(no_of_bins);
    std::vector<double> mzs = std::vector<double>(no_of_bins);
    // Resetting quantities
    double ar_av        = 0;
    double energy_av    = 0;
    double energy_sq_av = 0;
    double cv_average   = 0;
    double mx_av = 0;
    double my_av = 0;
    double mz_av = 0;
    for(int i=0; i<no_of_bins; i++)  // Loop over the bins
    {   // For each bin
        energies[i]    = 0;
        energies_sq[i] = 0;
        cvs[i]         = 0;
        mxs[i]         = 0;
        mys[i]         = 0;
        mzs[i]         = 0;
        if(LADYBUG)    cout << "i = " << i << "; energy_av before loop: " << energy_av << endl;
        // Setting vectors

        for(int j=0; j<mcsteps_inbin; j++)    // Loop over mcsteps in bin
        {   // For each mcstep
            mcstepf_metropolis(beta, generator_u, generator_v, generator_n, generator_prob, distribution_prob, distribution_u, distribution_v, distribution_n);
            // acceptancerate
            acceptancerates[i] += acceptancerate;
            ar_av += acceptancerate;
            // energy
            energies[i]    += energy_old;    // Storing to get the standard deviation // Something odd here
            energies_sq[i] += energy_old*energy_old;
            energy_av      +=energy_old;
            energy_sq_av   +=energy_old*energy_old;
            // Magnetization
            double mx = 0;
            double my = 0;
            double mz = 0;
            for(int k=0; k<N; k++)
            {
                mx+= mylattice.sites[k].spinx;
                my+= mylattice.sites[k].spiny;
                mz+= mylattice.sites[k].spinz;
            }
            mx = mx/N;
            my = my/N;
            mz = mz/N;

            mx_av += mx;
            my_av += my;
            mz_av += mz;

            mxs[i] += mx;
            mys[i] += my;
            mzs[i] += mz;

            //Print to bigFile
            bigFile << beta << " " << energy_old << " " << energy_sq_av << " " << mx << " " << my << " " << mz << endl;

            // Some sort of measurement of the magnetization... How to do this when we have a continuous spin?
        }  // End loop over mcsteps
        // For every bin, we find the following quantities:
        mxs[i]   = mxs[i]/mcsteps_inbin;
        mys[i]   = mys[i]/mcsteps_inbin;
        mzs[i]   = mzs[i]/mcsteps_inbin;
        energies[i]         = energies[i]/mcsteps_inbin;
        energies_sq[i]      = energies_sq[i]/mcsteps_inbin;
        cvs[i]   = beta*beta*(energies_sq[i]-energies[i]*energies[i]);

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

    double cv_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    cv_stdv += (cvs[l]-cv_average)*(cvs[l]-cv_average);
    cv_stdv = sqrt(cv_stdv/(no_of_bins*(no_of_bins-1)));

    //-----------------------//Magnetization//----------------------//
    mx_av = mx_av/(mcsteps_inbin*no_of_bins);
    my_av = my_av/(mcsteps_inbin*no_of_bins);
    mz_av = mz_av/(mcsteps_inbin*no_of_bins);

    // Error in the magnetization //
    double mx_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    mx_stdv += (mxs[l]-mx_av)*(mxs[l]-mx_av);
    mx_stdv = sqrt(mx_stdv/(no_of_bins*(no_of_bins-1)));

    double my_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    my_stdv += (mys[l]-my_av)*(mys[l]-my_av);
    my_stdv = sqrt(my_stdv/(no_of_bins*(no_of_bins-1)));

    double mz_stdv = 0;
    for(int l=0; l<no_of_bins; l++)    mz_stdv += (mzs[l]-mz_av)*(mzs[l]-mz_av);
    mz_stdv = sqrt(mz_stdv/(no_of_bins*(no_of_bins-1)));

    // Printing
    allFile << beta << " " << energy_av << " " << E_stdv << " " << energy_sq_av << " " << Esq_stdv << " " << cv << " " << cv_stdv << " " <<  mx_av ;
    allFile << " " << mx_stdv << " " << my_av << " " << my_stdv << " " << mz_av << " " << mz_stdv << " " << ar_av << " " << ar_stdv << endl;
    //print.printing_everybin(beta, energy_av, E_stdv, energy_sq_av, Esq_stdv, cv, cv_stdv, mx_av, mx_stdv, my_av, my_stdv, mz_av, mz_stdv);

    // Guess I should have stuff here instead. Print once for every beta.
    endtime = clock();
    total_time = (endtime - starttime)/(double) CLOCKS_PER_SEC;
    cout << "Time MC steps and measurements: "  << total_time << endl;
    cout << "Done with the Monte Carlo procedure this time around" << endl;
}


void MonteCarlo::mcstepf_metropolis(double beta, std::default_random_engine generator_u, std::default_random_engine generator_v, std::default_random_engine generator_n, std::default_random_engine generator_prob,  std::uniform_real_distribution<double> distribution_prob, std::uniform_real_distribution<double> distribution_u, std::uniform_real_distribution<double> distribution_v, std::uniform_int_distribution<int> distribution_n)
{   // Include a counter that measures how many 'flips' are accepted. But what to do with it? Write to file?

    bool HUMBUG = false;  // The humbug is defeated.
    if(HUMBUG)    cout << "In mcstepf_metropolis" << endl;
    double changes = 0;
    if(HUMBUG)   cout << "In mcstepf. Looping over spins now" << endl;
    for(int n=0; n<N; n++)
    {
        if(HUMBUG)    cout << "Inside loop in mcstepf. n = " << n << endl;
        double energy_diff = 0;

        int k = distribution_n(generator_n);
        if(HUMBUG)    cout << "Random spin k drawn. k = " << k << endl;

        double sx = mylattice.sites[k].spinx;
        double sy = mylattice.sites[k].spiny;
        double sz = mylattice.sites[k].spinz;

        if(HUMBUG)    cout << "Components of spin " << k << " accessed" << endl;

        // Changing the spin (tentatively):
        /*
        double u = distribution_u(generator_u);
        double v = distribution_v(generator_v);
        if(HUMBUG)    cout << "Have drawn random numbers in mcstepf" << endl;


        double theta = acos(1.0-2.0*u);
        double phi = 2.0*M_PI*v;

        double sintheta = sin(theta);
        double sx_t = sintheta*cos(phi);
        double sy_t = sintheta*sin(phi);
        //double sz_t = sqrt(1-sintheta*sintheta);
        double sz_t = cos(theta);
        */

        double zetasq = 2.0;  // zeta squared. Value to initialize loop.
        double zeta1;
        double zeta2;
        while(zetasq>=1.0) // Is this correct? Could equal 1 too, right?
        {
            double r_1 = distribution_u(generator_u);
            double r_2 = distribution_v(generator_v);

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
        // Check normalization
        //if(DEBUG)    cout << "Normalized? S^2 = " << (sx_t*sx_t + sy_t*sy_t + sz_t*sz_t) << endl;

        if(HUMBUG)    cout << "Have made a uniform spherical distribution using them" << endl;
        // Energy contribution after spin change

        if(sianisotropy)
        {
            if(HUMBUG)    cout << "In sianisotropy in mcstepf" << endl;
            double Dix = mylattice.sites[k].Dix;
            double Diy = mylattice.sites[k].Diy;
            double Diz = mylattice.sites[k].Diz;
            energy_diff += Dix*(sx*sx - sx_t*sx_t) + Diy*(sy*sy-sy_t*sy_t)+ Diz*(sz*sz-sz_t*sz_t); // This is - originally
            //energy_diff -= sianisotropy_energy(k, sx, sy, sz, mylattice);
        }
        if(magfield)
        {
            if(HUMBUG)    cout << "In magfield in mcstepf" << endl;
            double hx = mylattice.sites[k].hx;
            double hy = mylattice.sites[k].hy;
            double hz = mylattice.sites[k].hz;
            energy_diff += hx*(sx-sx_t) + hy*(sy-sy_t) + hz*(sz-sz_t);
            //energy_diff -= magfield_energy(k, sx, sy, sz, mylattice);
        }
        if(isotropic)
        {
            if(HUMBUG)    cout << "In isotropic in mcstepf" << endl;
            if(HUMBUG)    cout << "no_of_neighbours = " << no_of_neighbours << endl;
            double partnerspinx = 0;
            double partnerspiny = 0;
            double partnerspinz = 0;
            for(int j=0; j<no_of_neighbours; j++)
            {
                int l = mylattice.sites[k].bonds[j].siteindex2;
                if(HUMBUG)    cout << "Spin no. " << l << " chosen." << endl;

                double J = mylattice.sites[k].bonds[j].J;
                double sxk = mylattice.sites[l].spinx;
                double syk = mylattice.sites[l].spiny;
                double szk = mylattice.sites[l].spinz;
                partnerspinx += J*sxk;
                partnerspiny += J*syk;
                partnerspinz += J*szk;
            }
            if(HUMBUG)    cout << "Out of that blasted loop!" << endl;
            energy_diff += partnerspinx*(sx_t-sx) + partnerspiny*(sy_t-sy) + partnerspinz*(sz_t-sz);
        }
        if(dm)
        {
            if(HUMBUG)    cout << "In dm in mcstepf" << endl;
            for(int j=0; j<no_of_neighbours; j++)
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
            //energy_diff -= dm_energy(k, sx, sy, sz, mylattice);
        }
        if(HUMBUG)    cout << "Done with dm in mcstepf" << endl;
        //cout << "Contribution from energy before: " << energy_diff << endl;

        // This should work, but there is probably some error here...
        double energy_new = energy_old + energy_diff;
        //cout << "energy_diff: " << energy_diff << ";  Difference in spin: [" << (sx-sx_t) << "," << (sy-sy_t) << "," << (sz-sz_t) << "]" << endl;
        // Or should I just test if energy_new < 0? ... May have to find energy_new anyways...

        // Updating the energy and the state according to Metropolis
        if(energy_new <= energy_old)
        {
            mylattice.sites[k].spinx = sx_t;
            mylattice.sites[k].spiny = sy_t;
            mylattice.sites[k].spinz = sz_t;
            //cout << "Energy decreased. deltaS = [" << mylattice.sites[k].spinx-sx << "," << mylattice.sites[k].spiny-sy << "," << mylattice.sites[k].spinz-sz << "]." << " deltaE =  " << energy_diff  << endl;
            changes+=1;
            //cout << "Percentage of hits: " << changes/(n+1)  << endl;
            //cout << "ENERGY DECREASED! energy_old = " << energy_old << "; energy_diff = " << energy_diff << "; energy_new = " << energy_new << endl;
            if(HUMBUG)    cout << "ENERGY DECREASED!" << endl;
            //cout << "Resetting energy. Energy_old = " << energy_old << "; energy_new = " << energy_new << endl;
            energy_old = energy_new; // Update energy
            //cout << "Energy should be reset now. energy_old = " << energy_old << "; energy_new = " << energy_new << endl;
        }
        else
        {
            double prob = exp(-beta*(energy_new-energy_old));
            double drawn = distribution_prob(generator_prob);
            if(HUMBUG)    cout << "Suggesting energy increase. Probability of success: " << prob << "; number drawn: " << drawn << endl;
            if(drawn<prob)
            {
                mylattice.sites[k].spinx = sx_t;
                mylattice.sites[k].spiny = sy_t;
                mylattice.sites[k].spinz = sz_t;
                //cout << "Energy increased. deltaS = [" << mylattice.sites[k].spinx-sx << "," << mylattice.sites[k].spiny-sy << "," << mylattice.sites[k].spinz-sz << "]" << energy_diff << endl;
                //cout << "ENERGY INCREASED! energy_old = " << energy_old << "; energy_diff = " << energy_diff << "; energy_new = " << energy_new << endl;
                changes+=1;
                if(HUMBUG)    cout << "Success" << endl;
                //cout << "Percentage of hits: " << changes/(n+1)  << endl;
                energy_old = energy_new;  // Update energy
            }
        }
    }
    acceptancerate = changes/N; // Write the percentage of hits to file.
}

void MonteCarlo::endsims()
{
    //print.closeAllFiles();
    if(allFile.is_open())      allFile.close();
    if(bigFile.is_open())      bigFile.close();
}
