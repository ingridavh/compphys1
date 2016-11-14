#include <iostream>
#include "spinsystem.h"
#include <armadillo>
#include "mpi.h"
#include <string>

using namespace std;


int main(int argc, char* argv[])
{
    int nproc;
    int myrank;

    ofstream ofile;
    string filename=argv[1];
    ofile.open(filename.c_str(), ofstream::out);

    if(!ofile.is_open()) {
            ofile.open(filename.c_str(), ofstream::out);
            if(!ofile.good()) {
                cout << "Error opening file " << filename << ". Aborting!" << endl;
                terminate();
            }
    }




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

    spinsystem mysystem(n_spins);

    arma::mat info_matrix = mysystem.go(filename, mcs, initial_my_temp, final_my_temp, temp_step, myrank, nproc);


    int size = (final_my_temp - initial_my_temp)/temp_step;
    //Write to file in the right order
    MPI_Barrier (MPI_COMM_WORLD);
    for (int i = 0; i < nproc; i++) {
        if (i == myrank) {
            cout << "Hello world, I have  rank " << myrank <<
                    " out of " << nproc << endl;
            for (int i=0; i < size; i++)
            {
                for (int j = 0; j < 5; j++){
                    ofile << info_matrix(i,j) << " , " ;
                }
                ofile << endl;
            }


        MPI_Finalize();
        }

    }

    ofile.close();


    return 0;

//    //Testing MPI------------------
//    if (myrank==0)
//    {
//        cout << "HEY! We have " << nproc << " available processors" << endl;
//    }
//    cout << "Worker bee no " << myrank << " says the initial temperature is before 'go' " <<
//            initial_temp << endl;

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



