#include <iostream>
#include "spinsystem.h"
#include <armadillo>
#include "mpi.h"
#include <string>
#include <sstream>

using namespace std;

void calc_PE()
    {

    int myrank = 0;
    int nproc = 1;

    int n_spins = 20;
    int mcs = 1E5;
    double itemp1 = 1.0;
    double ftemp1 = 1.0;
    double itemp24 = 2.4;
    double ftemp24 = 2.4;

    double temp_step = 1.0;

    ofstream ofile1;
    ofstream ofile24;
    string filename1= "PE1.txt";
    string filename24 = "PE24.txt";

    spinsystem mysystem(n_spins);

    if(!ofile1.is_open()) {
            ofile1.open(filename24.c_str(), ofstream::out);
            if(!ofile1.good()) {
                cout << "Error opening file " << filename24 << ". Aborting!" << endl;
                terminate();
            }
    }

    mysystem.go(filename24, mcs, itemp24, ftemp24, temp_step, myrank, nproc);
    ofile1.close();
}


int main(int argc, char* argv[])
{
    int nproc;
    int myrank;

    ofstream ofile;
    string filename=argv[1];

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

    if (myrank == 0 && argc <= 5) {
        cout << "Bad Usage: " << argv[0] <<
                " read output file, Number of spins, mc cycles, initial temp, final temp and temp step" << endl;
        exit(1);
    }

    //Tell all nodes about initial variables
    int n_spins = atoi(argv[2]);
    int mcs = atoi(argv[3]);
    double initial_temp = atof(argv[4]);
    double final_temp = atof(argv[5]);
    double temp_step = atof(argv[6]);

    //Testing different method

    double dT = (final_temp - initial_temp) / (double) nproc;

    //Assign temperatures to nodes
    double initial_my_temp = initial_temp + myrank*dT;
    double final_my_temp = initial_my_temp + dT - temp_step;

    cout << "Worker bee number " << myrank << " of " << nproc << " will work from temp " <<
            initial_my_temp << " to " << final_my_temp << endl;

    if (myrank == nproc -1) initial_my_temp += temp_step;

    //Open one file per node
    ostringstream oss;
        oss << filename << n_spins << "_rank" << myrank << ".txt";
    string myfilename = oss.str();

    ofile.open(myfilename.c_str(), ofstream::out);

    if(!ofile.is_open()) {
            ofile.open(myfilename.c_str(), ofstream::out);
            if(!ofile.good()) {
                cout << "Error opening file " << myfilename << ". Aborting!" << endl;
                terminate();
            }
    }

    cout << "Bzz, rank " << myrank << " reporting filename " << myfilename << endl;

    spinsystem mysystem(n_spins);

    mysystem.go(myfilename, mcs, initial_my_temp, final_my_temp, temp_step, myrank, nproc);

    cout << "Bzzz, rank " << myrank << " has finished writing to file " << myfilename << endl;

    ofile.close();

    MPI_Finalize();

//    calc_PE();
    return 0;
}














