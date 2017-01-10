#include "montecarlo.h"

MonteCarlo::MonteCarlo(bool isotropic, bool sianisotropy, bool magfield, bool dm, Lattice mylattice)
{
    // Setting the bools
    this->isotropic = isotropic;
    this->sianisotropy = sianisotropy;
    this->magfield = magfield;
    this->dm = dm;

    // The lattice  // Could just extract the bools from here...
    this->mylattice = mylattice;
    this->N = mylattice.N;

    acceptancerate = 0;
}

void MonteCarlo::initialize_energy()
{
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
}

void MonteCarlo::runmetropolis(string filenamePrefix)
{
    bool HUMBUG  = false;
    bool LADYBUG = false;
    // Printing //Send prefix in?
    // Opening file to print to
    ofstream printFile;
    //string filenamePrefix = "test";
    //string filenamePrefix = "fcc10t10t10_sian1_1_1_beta2p5";
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_cspinMC.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    printFile.open(filename);
    delete filename;

    ofstream bigFile;
    char *filenameb = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenameb, "%s_dev_energyav.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    bigFile.open(filenameb);
    delete filenameb;

    // File for storing the acceptance rates for each MC-step
    ofstream arFile;
    char *filenamea = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamea, "%s_acceptancerate.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    arFile.open(filenamea);

    // Random generators
    std::default_random_engine generator_u;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_u(0,1);

    std::default_random_engine generator_v;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_v(0,1);

    std::default_random_engine generator_prob;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution_prob(0,1);

    // For index. This is given helical boundary conditions, then I only need one index
    std::default_random_engine generator_n;
    std::uniform_int_distribution<int> distribution_n(0,N-1);

    initialize_energy();

    // Equilibration steps
    for(int i=0; i<eqsteps; i++)
    {
        if(HUMBUG)    cout << "In equilibration steps loop, i = " << i << endl;
        mcstepf_metropolis();
        if(i<11)      cout << "i = " << i << ", energy: " << energy_old << endl;
    }

    // Monte Carlo steps and measurements
    for(int i=0; i<no_of_bins; i++)
    {   // For each bin
        double energy_av = 0;
        if(LADYBUG)    cout << "i = " << i << "; energy_av before loop: " << energy_av << endl;
        std::vector<double> energies = std::vector<double>(mcsteps_inbin);
        for(int j=0; j<mcsteps_inbin; j++)
        {    // For each mcstep, acceptancerate
            mcstepf_metropolis(generator_u, generator_v, generator_n, generator_prob, distribution_prob, distribution_u, distribution_v, distribution_n);
            // energy
            energies[j] = energy_old;    // Storing to get the standard deviation
            energy_av +=energy_old;
            arFile << acceptancerate << endl;  // Maybe I should change what I send in to this one. Possibly change how I handle acceptancerate as well.
            bigFile << i << " " << energy_av/(j+1) << endl;

            // Some sort of measurement of the magnetization... How to do this when we have a continuous spin?
        }
        // Energy
        energy_av = energy_av/mcsteps_inbin;

        // Error in the energy //
        double E_stdv = 0;
        for(int i=0; i<no_of_bins; i++)    E_stdv += (energies[i]-energy_av)*(energies[i]-energy_av);
        E_stdv = sqrt(E_stdv/(no_of_bins*(no_of_bins-1)));

        printFile << energy_av << " " << E_stdv << endl;
        // Print to file
    }
}


void MonteCarlo::mcstepf_metropolis( std::default_random_engine generator_u, std::default_random_engine generator_v, std::default_random_engine generator_n, std::default_random_engine generator_prob,  std::uniform_real_distribution<double> distribution_prob, std::uniform_real_distribution<double> distribution_u, std::uniform_real_distribution<double> distribution_v, std::uniform_int_distribution<int> distribution_n)
{   // Include a counter that measures how many 'flips' are accepted. But what to do with it? Write to file?
    bool DEBUG = false;
    bool HUMBUG = false;  // The humbug is defeated. I think...
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
        double u = distribution_u(generator_u);
        double v = distribution_v(generator_v);
        if(HUMBUG)    cout << "Have drawn random numbers in mcstepf" << endl;

        double theta = acos(1.0-2.0*u);
        double phi = 2.0*M_PI*v;

        double sx_t = sin(theta)*cos(phi);
        double sy_t = sin(theta)*sin(phi);
        double sz_t = cos(theta);

        if(HUMBUG)    cout << "Have made a uniform spherical distribution using them" << endl;
        // Energy contribution after spin change

        if(sianisotropy)
        {
            if(HUMBUG)    cout << "In sianisotropy in mcstepf" << endl;
            double Dix = mylattice.sites[k].Dix;
            double Diy = mylattice.sites[k].Diy;
            double Diz = mylattice.sites[k].Diz;
            energy_diff += Dix*(sx*sx - sx_t*sx_t) + Diy*(sy*sy-sy_t*sy_t)+ Diz*(sz*sz-sy_t*sy_t); // This is - originally
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
                // I could, alternatively, just store the index
                // of bonds. But then I need to organize the bonds
                // and coordinate them with the sites. List of
                //sites that points to the bonds and vice versa.
                double J = mylattice.sites[k].bonds[j].J;
                double sxk = mylattice.sites[l].spinx;
                double syk = mylattice.sites[l].spiny;
                double szk = mylattice.sites[l].spinz;
                partnerspinx += J*sxk;
                partnerspiny += J*syk;
                partnerspinz += J*szk;
            }
            if(HUMBUG)    cout << "Out of that blasted loop!" << endl;
            energy_diff = partnerspinx*(sx_t-sx) + partnerspiny*(sy_t-sy) + partnerspinz*(sz_t-sz);
        }
        if(dm)
        {
            if(HUMBUG)    cout << "In dm in mcstepf" << endl;
            for(int j=0; j<no_of_neighbours; j++)
            {
                int l = mylattice.sites[k].bonds[j].siteindex2; // Hope I can actually get to this value.
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
            if(DEBUG)    cout << "ENERGY DECREASED!" << endl;
            energy_old = energy_new; // Update energy
        }
        else
        {
            double prob = exp(-beta*(energy_new-energy_old));
            double drawn = distribution_prob(generator_prob);
            if(DEBUG)    cout << "Suggesting energy increase. Probability of success: " << prob << "; number drawn: " << drawn << endl;
            if(drawn<prob)
            {
                mylattice.sites[k].spinx = sx_t;
                mylattice.sites[k].spiny = sy_t;
                mylattice.sites[k].spinz = sz_t;
                //cout << "Energy increased. deltaS = [" << mylattice.sites[k].spinx-sx << "," << mylattice.sites[k].spiny-sy << "," << mylattice.sites[k].spinz-sz << "]" << energy_diff << endl;
                //cout << "ENERGY INCREASED! energy_old = " << energy_old << "; energy_diff = " << energy_diff << "; energy_new = " << energy_new << endl;
                changes+=1;
                if(DEBUG)    cout << "Success" << endl;
                //cout << "Percentage of hits: " << changes/(n+1)  << endl;
                energy_old = energy_new;  // Update energy
            }
        }
    }
    acceptancerate = changes/N; // Write the percentage of hits to file.
    return mcstepf_return;
}
