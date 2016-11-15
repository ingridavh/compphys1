#include <iostream>
#include "spinsystem.h"
#include <armadillo>
#include "mpi.h"
#include <string>
#include <sstream>

using namespace std;


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

    return 0;

//    //Testing MPI------------------
//    spinsystem mysystem(n_spins);

//    double time_start, time_end, total_time;
//    time_start = MPI_Wtime();

//    //Run the simulation
//    //char outfilename, mcs, initial_temp, final_temp, temp_step
//    mysystem.go(filename, mcs, initial_temp, final_temp, temp_step, myrank, nproc);

//    time_end = MPI_Wtime();
//    total_time = time_end - time_start;

//    if (myrank == 0){
//        cout << "Time = " << total_time << " on number of processors: " <<
//                nproc << endl;
//    }

//    MPI_Finalize();



    return 0;
}



