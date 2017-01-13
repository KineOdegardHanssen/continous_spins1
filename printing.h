#ifndef PRINTING_H
#define PRINTING_H
#include <fstream>
#include <vector>
#include <string>
using std::ofstream; using std::string;

class Printing
{
public:
    ofstream    allFile;
    ofstream    bigFile;
    ofstream    arFile;
    string      filenamePrefix;

    // Intitialization and destructor
    Printing();
    Printing(string filenamePrefix);
    ~Printing();
    void closeAllFiles();

    void open_allFile();
    void open_bigFile();
    void open_arFile();
    void open_all();

    // Printing a header line on the files. Still prototype/placeholder
    void print_header_allFile(double N);

    // Printing in the Monte Carlo procedure
    void printing_everybin(double energy_av, double E_stdv, double energy_sq_av, double Esq_stdv, double cv, double cv_stdv, double mx_av, double mx_stdv, double my_av, double my_stdv, double mz_av, double mz_stdv);
    void printing_everystep(double beta, double energy_old, double energy_sq_av, double cv, double mx, double my, double mz);
    void printing_acceptancerates(double acceptancerate);

};   // End of Printing class declaration


#endif // PRINTING_H
