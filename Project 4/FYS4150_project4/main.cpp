#include <iostream>
#include "spinsystem.h"

using namespace std;

//For command line change to (int argc, char* argv[])

int main()
{
    cout << "This is a spin system" << endl;
    //Needed parameters are;
    //char outfilename, mcs, initial_temp, final_temp, temp_step

    string outfilename = "results";

    //WHAT IS IDUM AND MCS?
    long idum;
    int mcs = 10;
    double initial_temp = 1.0;
    double final_temp = 4.0;
    double temp_step = 0.001;

    int n_spins = 2;

    spinsystem A(n_spins);
    A.go(outfilename, mcs, initial_temp,final_temp, temp_step);


//---------------Read input parameters from command line (later)

//    char *outfilename;
//    long idum;
//    int **spin_matrix, n_spins, mcs;


//    double w[17], average[5], initial_temp, final_temp, E, M, temp_step;

//    // Read in output file, abort if there are too few command-line arguments
//    if( argc <= 1)
//    {
//        cout << "Bad usage" << argv[0] <<
//                " read also output file on same line" << endl;
//        exit(1);
//    }
//    else
//    {
//        outfilename=argv[1];
//    }

//----------------End-----------------------------------------------

    cout << "This is the end" << endl;
    return 0;
}

