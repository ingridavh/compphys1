#include <iostream>
#include "spinsystem.h"
#include <armadillo>
#include "mpi.h"
#include <string>

using namespace std;


int main(int argc, char* argv[])
{
    int NProcesses;
    int rankProcess;
    string filename;
    int N_spins;
    int mcs;
    double initial_temp;
    double final_temp;
    double temp_step;

    // MPI initializations

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &NProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &rankProcess);
    if (rankProcess == 0 && argc <= 5) {
        cout << "Bad Usage: " << argv[0] <<
                " read output file, NUmber of spins, mc cycles, initial temp, final temp and temp step" << endl;
        exit(1);
    }
    if ((rankProcess==0) && (argc > 1)){
        filename=argv[1];
        N_spins = atoi(argv[2]);
        mcs = atoi(argv[3]);
        initial_temp = atof(argv[4]);
        final_temp = atof(argv[5]);
        temp_step = atof(argv[6]);
    }

    // Broadcast to all nodes common variables since only
    //master node reads from command line

    MPI_Bcast(&mcs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&initial_temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&final_temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temp_step, 1, MPI_INT, 0, MPI_COMM_WORLD);

    spinsystem mysystem(N_spins);

    double time_start, time_end, total_time;
    time_start = MPI_Wtime();


    //Run the simulation
    mysystem.go(filename, mcs, initial_temp, final_temp, temp_step, rankProcess, NProcesses);


    time_end = MPI_Wtime();
    total_time = time_end - time_start;


    if (rankProcess == 0){
        cout << "Time = " << total_time << " on number of processors: " <<
                NProcesses << endl;
    }

    MPI_Finalize();

    //Needed parameters are;
    //char outfilename, mcs, initial_temp, final_temp, temp_step

//    string outfilename = "prettywoman_.txt";
//    double n_cycl = 1E5;
//    double initial_temp = 2.4;
//    double final_temp = 2.4;
//    double temp_step = 1;
//    int n_spins = 20;

//    //Open output file
//    ofstream myfile;
//    myfile.open(outfilename.c_str(), ofstream::out);
//    if (!myfile.is_open()){
//        cout << "FILE " << outfilename.c_str() << " not opened, aborting!" << endl;
//        terminate();
//    }


//    spinsystem mysystem(n_spins);
//    mysystem.go(outfilename, n_cycl, initial_temp, final_temp, temp_step);

    return 0;
}



