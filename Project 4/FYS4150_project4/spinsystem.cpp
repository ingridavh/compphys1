#include "spinsystem.h"
#include "mpi.h"
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <cstdlib>

using namespace std;

spinsystem::spinsystem(int n_spins) :
    m_J(0),
    m_spin_matrix(arma::mat(n_spins, n_spins)),
    m_n_spins(n_spins),
    m_no_accept(0)
{
}

//"Body" of the class which uses the Metropolis algorithm
//to calculate expectation values of the energy and magnetization

void spinsystem::go(string outfilename, int mcs, double initial_temp,
                    double final_temp, double temp_step, int myrank, int nproc){
    m_mcs = mcs;
    m_finaltemp = final_temp;

    if(!m_ofile.is_open()) {
            m_ofile.open(outfilename.c_str(), ofstream::out);
            if(!m_ofile.good()) {
                    cout << "Error opening file " << outfilename << ". Aborting!" << endl;
                terminate();
            }
        }

    m_ofile << "Totally awesome! The program is working! Now numbers: " << endl;
    m_ofile << "Temp \t Energy \t Cv \t abs(Magnetization) \t X" << endl;

    double *w = new double[17];
    m_no_accept = 0;

    for (double temp = initial_temp; temp <= final_temp; temp += temp_step){

        //Reset average values
        for (int i=0; i<5; i++){
            m_average[i] = 0;
        }

        // Initialize energy and magnetization
        double E = 0;
        double M = 0;

        //Set up array for possible energy changes
        for( int de = -8; de <= 8; de++) w[de+8] = 0;
        for( int de = -8; de <= 8; de+=4) w[de+8] = exp(-de/temp);

        //Initialize array for expectation values
        initialize(E, M, temp);


        //Start Monte Carlo computation
        for (int cycles = 1; cycles <= m_mcs; cycles++){
            Metropolis(E, M, w, temp);

            //update expectation values
            m_average[0] += E;
            m_average[1] += E*E;
            m_average[2] += M;
            m_average[3] += M*M;
            m_average[4] += fabs(M);

//            m_ofile << E << endl;

            //For monitoring values as functions of Monte Carlo cycles

//            if (cycles >= 100 && cycles%100 ==0 ){
//                for (int i = 0; i < 5; i++) m_average[i] /= cycles;
//                m_ofile << cycles << " , " << m_average[0] << " , " << m_average[4] << " , " << m_no_accept << endl;
//                for (int i = 0; i < 5; i++) m_average[i] *= cycles;
//                }
//            }
        }

        for (int i =0; i<5; i++) m_average[i] /= m_mcs;
        output(mcs, temp);

    }
}

//set up spin matrix and initial magnetization
void spinsystem::initialize(double& E, double& M, double temp){
    arma::vec spins = {1, -1};
    for (int y= 0; y < m_n_spins; y++){
        for (int x=0; x < m_n_spins; x++){
            //Set the spin orientation for ground state. If the temp is low enough we set it to ground state
            if (temp > 1.5){
                double r = rand() %2;
                m_spin_matrix(y,x) = spins[r];
            }
            else{
                m_spin_matrix(y,x) = 1;
            }

            M += (double) m_spin_matrix(y,x);
        }
    }

    //set up initial energy
    for (int y = 0; y < m_n_spins; y++){
        for (int x = 0; x < m_n_spins; x++){
            E -= (double) m_spin_matrix(y,x)*
                    (m_spin_matrix(periodic(y, m_n_spins, -1),x) +
                    m_spin_matrix(y,periodic(x, m_n_spins, -1)));
        }
    }
    m_average[0] += E;
    m_average[1] += E*E;
    m_average[2] += M;
    m_average[3] += M*M;
    m_average[4] += fabs(M);
}

//Use periodic boundary conditions
double spinsystem::periodic(int i, int limit, int add){
    return (i+limit+add) % (limit);
}

//The Metropolis algorithm for calculating
//expectation values
void spinsystem::Metropolis(double& E, double &M, double *w, double temp){
    //loop over all spins
    for (int y = 0; y < m_n_spins; y++){
        for (int x=0; x < m_n_spins; x++){
            //Find random position
            int ix = (int) (rand() % m_n_spins);
            int iy = (int) (rand() % m_n_spins);
            m_spin_matrix(iy, ix) *= -1;

            int deltaE = - 2*m_spin_matrix(iy,ix)*
                    (m_spin_matrix(iy,periodic(ix, m_n_spins, -1))+
                    m_spin_matrix(periodic(iy, m_n_spins, -1),ix)+
                    m_spin_matrix(iy,periodic(ix, m_n_spins, 1))+
                    m_spin_matrix(periodic(iy, m_n_spins,1),ix));

            //Perform the Metropolis test with a random number between 0 and 1
            double r = rand()/ (double)RAND_MAX;
            if (r <= exp(-deltaE/temp))
            {
                //flip one spin and accept new spin configuration
                //update energy and magnetization
                M += (double) 2*m_spin_matrix(iy,ix);
                E += (double) deltaE;
                m_no_accept += 1;
            } else {
                m_spin_matrix(iy, ix) *= -1;
            }
        }
    }
}

//Function that writes data to output file
void spinsystem::output(int mcs, double temperature){
    //Per spin
    double temp = temperature;

    double Energy = m_average[0];
    double Energy2 = m_average[1];
    double Magnetization = m_average[2];
    double Magnetization2 = m_average[3];
    double Magnetization_abs = m_average[4];
    double E_variance = (Energy2 - Energy*Energy)/(m_n_spins*m_n_spins);
    double M_variance = (Magnetization2 - Magnetization_abs*Magnetization_abs)/(m_n_spins*m_n_spins);
    double Cv = (1./(temp*temp))*E_variance;
    double X = (1./temp)*M_variance;

    //Print to an outfile
    m_ofile << temp << " \t " << Energy/(m_n_spins*m_n_spins) << " \t " << Cv << " \t "
            << Magnetization_abs/(m_n_spins*m_n_spins) << " \t " << X <<endl;
}
