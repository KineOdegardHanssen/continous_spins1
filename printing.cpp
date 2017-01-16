#include "printing.h"


Printing::Printing()
{
}

Printing::Printing(string filenamePrefix)
{ // Initialisation of a printer object - it needs a prefix
    this->filenamePrefix = filenamePrefix;
    open_all();
} // End initialization


void Printing::givePrefix(string filenamePrefix)
{ // Initialisation of a printer object - it needs a prefix
    this->filenamePrefix = filenamePrefix;
} // End initialization
/*
Printing::~Printing()
{ // Destructor works by closing all files
    closeAllFiles();
} // End destructor
*/

void Printing::open_allFile()
{
    // Should I have this here? Or open as I please?
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_cspinMC.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    allFile.open(filename);
    delete filename;
}

void Printing::open_bigFile()
{
    char *filenameb = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenameb, "%s_dev_energyav.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    bigFile.open(filenameb);
    delete filenameb;
}

void Printing::open_arFile()
{
    // File for storing the acceptance rates for each MC-step
    char *filenamea = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamea, "%s_acceptancerate.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    arFile.open(filenamea);
    delete filenamea;

    // Should I also have a header line as not to confuse files?
}

void Printing::open_all()
{   // Call functions here as this is only to be done once per class instance
    open_allFile();
    open_bigFile();
    open_arFile();
}


void Printing::print_header_allFile(double N)
{   // Prototype for printing headers
    // Check this out more thouroughly later
    allFile << N << endl;
}


void Printing::printing_everybin(double beta, double energy_av, double E_stdv, double energy_sq_av, double Esq_stdv, double cv, double cv_stdv, double mx_av, double mx_stdv, double my_av, double my_stdv, double mz_av, double mz_stdv)
{
    allFile << energy_av << " " << E_stdv << " " << energy_sq_av << " " << Esq_stdv << " " << cv << " " << cv_stdv << " " <<  mx_av ;
    allFile << " " << mx_stdv << " " << my_av << " " << my_stdv << " " << mz_av << " " << mz_stdv << " " << beta << endl;
} // End printingPosition-function

void Printing::printing_everystep(double beta, double energy_old, double energy_sq_av, double mx, double my, double mz)
{
    bigFile << beta << " " << energy_old << " " << energy_sq_av << " " << mx << " " << my << " " << mz << endl;
}

void Printing::printing_acceptancerates(double beta, double acceptancerate)
{   // Decide how often I shoud print this.
    arFile << beta << " " << acceptancerate << endl;
}

void Printing::closeAllFiles()
{ // Function closing all open files
    if(allFile.is_open())      allFile.close();
    if(bigFile.is_open())      bigFile.close();
    if(arFile.is_open())       arFile.close();
} // End closeAllFiles-function
