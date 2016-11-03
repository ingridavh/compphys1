#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
    char *outfilename;
    long idum;
    int **spin_matrix, n_spins, mcs;
    double w[17], average[5], initial_temp, final_temp, E, M, temp_step;

    // Read in output file, abort if there are too few command-line arguments
    if( argc <= 1)
    {
        cout << "Bad usage" << argv[0] <<
                " read also output file on same line" << endl;
        exit(1);
    }
    else
    {
        outfilename=argv[1];
    }








    return 0;
}

