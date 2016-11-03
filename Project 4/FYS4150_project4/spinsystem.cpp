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
    if(!m_ofile.good()) {
            m_ofile.open(outfilename.c_str(), ofstream::out);
            if(!m_ofile.good()) {
                cout << "Error opening file " << outfilename << ". Aborting!" << endl;
                terminate();
            }
        }

    //HVA ER DENNE EGENTLIG?
    double *w = new double[17];
    double average[5];

    for (double temp = initial_temp; temp <= final_temp; temp += temp_step){
        // Initialize energy and magnetization
        double E = 0;
        double M = 0;

        //Set up array for possible energy changes
        for( int de = -8; de <= 8; de++) w[de+8] = 0;
        for( int de = -8; de <= 8; de+=4) w[de+8] = exp(-de/temp);

        //Initialize array for expectation values
        for( int i = 0; i < 5; i++) average[i] = 0.;
        initialize(E, M, temp);

        //Random starting point
        long idum = -1;

        //Start Monte Carlo computation
        for (int cycles = 1; cycles <= mcs; cycles++){
            Metropolis(idum, M, E, w);

            //update expectation values
            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);
        }
        // print results to file
        output(mcs, temp, average);
    }
    m_ofile.close();
    cout << "The average energy is " << average[0] << endl;
    cout << "The average magnetizatoin is " << average[2] << endl;

}


//set up spin matrix and initial magnetization
void spinsystem::initialize(double& E, double& M, double temp){
    for (int y= 0; y < m_n_spins; y++){
        for (int x=0; x < m_n_spins; x++){
            //Set the spin orientation for ground state
            if (temp < 1.5) m_spin_matrix(y,x) = 1;
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
void spinsystem::Metropolis(long& idum, double& E, double &M, double *w)
{
    //loop over all spins
    for (int y = 0; y < m_n_spins; y++){
        for (int x=0; x < m_n_spins; x++){
            //Find random position
            int ix = (int) (rand() % m_n_spins);
            int iy = (int) (rand() % m_n_spins);
            int deltaE = 2*m_spin_matrix(iy,ix)*
                    (m_spin_matrix(iy,periodic(ix, m_n_spins, -1))+
                    m_spin_matrix(periodic(iy, m_n_spins, -1),ix)+
                    m_spin_matrix(iy,periodic(ix, m_n_spins, 1))+
                    m_spin_matrix(periodic(iy, m_n_spins,1),ix));
            //Here we perform the Metropolis test
            //with a random number between 0 and 1
            if (rand()/RAND_MAX <= w[deltaE+8])
            {
                //flip one spin and accept new spin configuration
                m_spin_matrix(iy,ix) *= -1;
                //update energy and magnetization
                M += (double) 2*m_spin_matrix(iy,ix);
                E += (double) deltaE;
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
            if (i<1) Energy +=-m_J*(m_spin_matrix(i+1,j) + m_spin_matrix(m_n_spins-1,j));
            else Energy += -m_J*(m_spin_matrix(i-1,j) + m_spin_matrix(i+1,j));

            if (j<1) Energy +=-m_J*(m_spin_matrix(i,j+1) + m_spin_matrix(i,m_n_spins-1));
            else Energy += -m_J*(m_spin_matrix(i,j-1) + m_spin_matrix(i,j+1));
        }
    }
    return Energy;

}


//Function that writes data to output file
void spinsystem::output(int mcs, double temperature, double *average)
{
    //Function that prints results to file
    //divided by total number of cycles
    double norm = 1/((double) (mcs));
    double Eaverage = average[0]*norm;
    double E2average = average[1]*norm;
    double Maverage = average[2]*norm;
    double M2average = average[3]*norm;
    double Mabsaverage = average[4]*norm;

    //All expectation values are per spin, divide by 1/(m_n_spins*m_n_spins)
    double Evariance = (E2average - Eaverage*Eaverage)/(m_n_spins*m_n_spins);
    double Mvariance = (M2average - Maverage*Maverage)/(m_n_spins*m_n_spins);
    double M2variance = (M2average - Mabsaverage*Mabsaverage)/(m_n_spins*m_n_spins);
    //Why is this different?
    double Mvariance_weird = (M2average - Mabsaverage*Mabsaverage)/(m_n_spins*m_n_spins);

    //Print to an outfile
    m_ofile << "Energy, E variance, Magnetization, M variance" << endl;
    m_ofile << Eaverage << " , " << Eaverage << " , " << Maverage << " , " << Mvariance <<endl;
}
