#include <iostream>
#include "spinsystem.h"
#include <armadillo>

using namespace std;

int main()
{
    //Needed parameters are;
    //char outfilename, mcs, initial_temp, final_temp, temp_step

    string outfilename = "prettywoman_.txt";
    double n_cycl = 1E5;
    double initial_temp = 2.4;
    double final_temp = 2.4;
    double temp_step = 1;
    int n_spins = 20;

    //Open output file


    ofstream myfile;
    myfile.open(outfilename.c_str(), ofstream::out);
    if (!myfile.is_open()){
        cout << "FILE " << outfilename.c_str() << " not opened, aborting!" << endl;
        terminate();
    }


    spinsystem mysystem(n_spins);
    mysystem.go(outfilename, n_cycl, initial_temp, final_temp, temp_step);

    return 0;
}
