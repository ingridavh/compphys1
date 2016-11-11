#include <iostream>
#include "spinsystem.h"
#include <armadillo>

using namespace std;

int main()
{
    //Needed parameters are;
    //char outfilename, mcs, initial_temp, final_temp, temp_step

    string outfilename = "prettywoman.txt";
    int mcs = 1e6;
    double initial_temp = 1.0;
    double final_temp = 1.0;
    double temp_step = 1;
    int n_spins = 20;

//    spinsystem A(n_spins);
//    A.go(outfilename, mcs, initial_temp,final_temp, temp_step);

    //Function of number of Monte Carlo cycles

    int n_tests = 8;
    arma::vec cycles = arma::zeros<arma::vec>(n_tests);
    for (int i = 0; i < n_tests; i++) cycles(i) = i;

    for (int i = 1; i < 8; i++) cycles[i] = i;
    //Open output file

    string filename = "exp_mcs.txt";

    ofstream myfile;
    myfile.open(filename.c_str(), ofstream::out);
    if (!myfile.is_open()){
        cout << "FILE " << filename.c_str() << " not opened, aborting!" << endl;
        terminate();
    }

    for (int i=0; i<n_tests; i++)
    {
        double n_cycl = pow(10, cycles[i]);
        spinsystem mysystem(n_spins);
        mysystem.go(outfilename, n_cycl, initial_temp, final_temp, temp_step);
        arma::vec ex = mysystem.getExpecs();
        cout << "The number of cycles is " << n_cycl << ", and energy is " << ex(0)  << endl;
        myfile << n_cycl << " , " <<  ex(0) << " , " << ex(4) << endl;

    }

    return 0;
}
