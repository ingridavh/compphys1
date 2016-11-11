#include "spinsystem.h"
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <cstdlib>

using namespace std;

spinsystem::spinsystem(int n_spins) :
    m_J(0),
    m_spin_matrix(arma::mat(n_spins, n_spins)),
    m_n_spins(n_spins)
{

}

//"Body" of the class which uses the Metropolis algorithm
//to calculate expectation values of the energy and magnetization

void spinsystem::go(string outfilename, int mcs, double initial_temp, double final_temp, double temp_step)
{
    //Open the output file
    //m_ofile.open(outfilename.c_str(), ofstream::out);
    m_mcs = mcs;
    m_finaltemp = final_temp;

    if(!m_ofile.is_open()) {
            m_ofile.open(outfilename.c_str(), ofstream::out);
            if(!m_ofile.good()) {
                    cout << "Error opening file " << outfilename << ". Aborting!" << endl;
                terminate();
            }
        }

    double *w = new double[17];

    for (double temp = initial_temp; temp <= final_temp; temp += temp_step){
        // Initialize energy and magnetization
        double E = 0;
        double M = 0;

        //Set up array for possible energy changes
        for( int de = -8; de <= 8; de++) w[de+8] = 0;
        for( int de = -8; de <= 8; de+=4) w[de+8] = exp(-de/temp);

        //Initialize array for expectation values
        for (int i= 0; i<5; i++) m_average[i] = 0;
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
        }
        // print results to file
        output(mcs, temp);
    }
    m_ofile.close();
    for (int i = 0; i < 5; i++) m_average[i] /= mcs;
    double Cv = (1./(final_temp*final_temp))*(m_average[1] - m_average[0]*m_average[0]);
    double X = (1./final_temp)*(m_average[3] - m_average[2]*m_average[2]);

//    cout << "The average energy is " << m_average[0] << endl;
//    cout << "Cv " << Cv << endl;
//    cout << "The average magnetization is " << m_average[2] << endl;
//    cout << "X " << X << endl;
//    cout << "abs(M) " << m_average[4] << endl;

}

arma::vec spinsystem::getExpecs()
{
    arma::vec A;
    A.zeros(5);
    for (int i=0; i<5; i++) A[i] = m_average[i];
    return A;
}


//set up spin matrix and initial magnetization
void spinsystem::initialize(double& E, double& M, double temp){
    for (int y= 0; y < m_n_spins; y++){
        for (int x=0; x < m_n_spins; x++){
            //Set the spin orientation for ground state
            m_spin_matrix(y,x) = 1;
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
}


//Use periodic boundary conditions
double spinsystem::periodic(int i, int limit, int add)
{
    return (i+limit+add) % (limit);
}


//The Metropolis algorithm for calculating
//expectation values
void spinsystem::Metropolis(double& E, double &M, double *w, double temp)
{
    //loop over all spins
    for (int y = 0; y < m_n_spins; y++){
        for (int x=0; x < m_n_spins; x++){
            //Find random position
            int ix = (int) (rand() % m_n_spins);
            int iy = (int) (rand() % m_n_spins);
            m_spin_matrix(iy, ix) *= -1;
//            cout << "The random numbers iy and ix are " << iy << " " << ix << endl;
            int deltaE = - 2*m_spin_matrix(iy,ix)*
                    (m_spin_matrix(iy,periodic(ix, m_n_spins, -1))+
                    m_spin_matrix(periodic(iy, m_n_spins, -1),ix)+
                    m_spin_matrix(iy,periodic(ix, m_n_spins, 1))+
                    m_spin_matrix(periodic(iy, m_n_spins,1),ix));
            //Here we perform the Metropolis test
            //with a random number between 0 and 1
            double r = rand()/ (double)RAND_MAX;

            if (r <= exp(-deltaE/temp))
            {
                //flip one spin and accept new spin configuration
                //m_spin_matrix(iy,ix) *= -1;
                //update energy and magnetization
                M += (double) 2*m_spin_matrix(iy,ix);
                E += (double) deltaE;
            } else {
                m_spin_matrix(iy, ix) *= -1;
            }
        }
    }
}

void spinsystem::read_input (int&, int&, double&, double&, double&){
    //Function to read in data from screen
}


//Calculate the magnetization
double spinsystem::CalculateMagnetization(){
    return arma::accu(m_spin_matrix);
}

//Calculate the energy
double spinsystem::CalculateEnergy()
{
    double Energy = 0;
    for (int i=0; i<m_n_spins; i++){
        for (int j=0; j<m_n_spins; j++){
            //Periodic boundary conditions
            if (i<1) Energy -= m_J*(m_spin_matrix(m_n_spins-1,j));
            else Energy -= m_J*(m_spin_matrix(i-1,j));

            if (j<1) Energy -=m_J*( m_spin_matrix(i,m_n_spins-1));
            else Energy -= m_J*(m_spin_matrix(i,j-1));
        }
    }
    return Energy;

}


//Function that writes data to output file
void spinsystem::output(int mcs, double temperature)
{
    //Per spin
    for (int i = 0; i < 5; i++) m_average[i] /= mcs*m_n_spins*m_n_spins;
    double Energy = m_average[0];
    double Magnetization = m_average[2];
    double Magnetization_abs = m_average[4];
    double Cv = (1./(m_finaltemp*m_finaltemp))*(m_average[1] - m_average[0]*m_average[0]);
    double X = (1./m_finaltemp)*(m_average[3] - m_average[2]*m_average[2]);

    //Print to an outfile
    m_ofile << "Energy, Cv, Magnetization, abs(Magnetization), X" << endl;
    m_ofile << Energy << " , " << Cv/m_n_spins << " , " << Magnetization << " , " << Magnetization_abs << " , " << X <<endl;
}
