#ifndef MONTECARLO_H
#define MONTECARLO_H
#include <fstream>
#include <vector>
#include <lattice.h>
#include <site.h>    // Are these really neccessary... ?
#include <bond.h>

using namespace std;

class MonteCarlo
{
public:
    MonteCarlo(int eqsteps, int mcsteps_inbin, int no_of_bins, bool isotropic, bool sianisotropy, bool magfield, bool dm, Lattice mylattice);

    int N, eqsteps, mcsteps_inbin, no_of_bins;
    double energy_old, acceptancerate;
    bool isotropic, sianisotropy, magfield, dm;

    Lattice mylattice;

    void initialize_energy();
    void runmetropolis(string filenamePrefix);
    void mcstepf_metropolis(std::default_random_engine generator_u, std::default_random_engine generator_v, std::default_random_engine generator_n, std::default_random_engine generator_prob,  std::uniform_real_distribution<double> distribution_prob, std::uniform_real_distribution<double> distribution_u, std::uniform_real_distribution<double> distribution_v, std::uniform_int_distribution<int> distribution_n); //

};

#endif // MONTECARLO_H
