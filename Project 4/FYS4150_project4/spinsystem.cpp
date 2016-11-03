#include "spinsystem.h"
#include <iomanip>
#include <armadillo>
#include <fstream>
//Hva er dette??
#include "lib.h"

using namespace std;

spinsystem::spinsystem(int n_spins) :
    m_J(0),
    m_n_spins(n_spins)
{

}


//"Body" of the class which uses the Metropolis algorithm
//to calculate expectation values of the energy and magnetization

void spinsystem::go(char outfilename)
{
    //Open the output file
    ofile.open(outfilename);
    //Read in initial values: size of lattice, temp and cycles
    read_input(n_spins, mcs, initial_temp, final_temp, temp_step);
    //Initialize the spin matrix
    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    //Random starting point
    idum = -1;

    for (double temp = initial_temp; temp <= final_temp; temp += temp_step)
    {
        // Initialize energy and magnetization
        E = M = 0;
        //Set up array for possible energy changes
        for( int de = -8; de <= 8; de++) w[de+8] = 0;
        for( int de = -8; de <= 8; de+=4) w[de+8] = exp(-de/temp);

        //Initialize array for expectation values
        for( int i = 0; i < 5; i++) average[i] = 0.;
        initialize(E, M, spin_matrix, temp);

        //Start Monte Carlo computation
        for (int cycles = 1; cycles <= mcs; cycles++)
        {
            Metropolis(n_spins, idum, spin_matrix, E, M, w);

            //update expectation values
            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);
        }
        // print results to file
        output(n_spins, mcs, temp, average);
    }






//set up spin matrix and initial magnetization
void spinsystem::initialize(double& E, double& M, int **spin_matrix, temp)
{
    for (int y= 0; y < m_n_spins; y++)
    {
        for (int x=0; x < m_n_spins; x++)
        {
            //Set the spin orientation for ground state
            if (temp < 1.5) m_spin_matrix[y][x] = 1;
            M += (double) m_spin_matrix[y][x];
        }
    }

    //set up initial energy
    for (int y = 0; y < m_n_spins; y++)
    {
        for (int x = 0; x < m_n_spins; x++)
        {
            E -= (double) m_spin_matrix[y][x]*
                    (m_spin_matrix[periodic(y, m_n_spins, -1)][x] +
                    m_spin_matrix[y][periodic(x, m_n_spins, -1)]);
        }
    }

}


//Use periodic boundary conditions
double spinsystem::periodic(int i, int limit, int add)
{
    return (i+limit+add) % (limit);
}


//The Metropolis algorithm for calculating
void spinsystem::Metropolis(long& idum, int **m_spin_matrix,
                double& E, double &M, double *w)
{
    //loop over all spins
    for (int y = 0; y < m_n_spins; y++)
    {
        for (int x=0; x < m_n_spins; x++)
        {
            //Find random position
            int ix = (int) (ran1(&idum)*(double)m_n_spins);
            int iy = (int) (ran1(&idum)*(double)m_n_spins);
            int deltaE = 2*m_spin_matrix[ix][iy]*
                    (m_spin_matrix[iy][periodic(ix, m_n_spins, -1)]+
                    m_spin_matrix[periodic(iy, m_n_spins, -1)][ix]+
                    m_spin_matrix[iy][periodic(ix, m_n_spins, 1)]+
                    m_spin_matrix[periodic(iy, m_n_spins,1)][ix]);
            //Here we perform the Metropolis test
            if (ran1(&idum) <= w[deltaE+8])
            {
                //flip one spin and accept new sin configuration
                m_spin_matrix[iy][ix] *= -1;
                //update energy and magnetization
                M += (double) 2*m_spin_matrix[iy][ix];
                E += (double) deltaE;
            }
        }
    }
}

void spinsystem::read_input (int&, int&, double&, double&, double&)
{
    //Function to read in data from screen
}


//Calculate the magnetization
double spinsystem::CalculateMagnetization()
{
    M = arma::accu(m_spin_matrix);
    return M;
}

//Calculate the energy
double spinsystem::CalculateEnergy()
{
    Energy = 0;
    for (i=0; i<m_n_spins; i++)
    {
        for (j=0; j<m_n_spins; j++)
        {
            //Periodic boundary conditions
            if (i<1) Energy +=-J*(m_spin_matrix[i+1][j] + m_spin_matrix[m_n_spins-1][j]);
            else Energy += -J*(m_spin_matrix[i-1][j] + m_spin_matrix[i+1][j]);

            if (j<1) Energy +=-J*(m_spin_matrix[i][j+1] + m_spin_matrix[i][m_n_spins-1]);
            else Energy += -J*(m_spin_matrix[i][j-1] + m_spin_matrix[i][j+1]);
        }
    }

}




void spinsystem::go(char outfilename)
{
    //Open the output file
    ofile.open(outfilename);
    //Read in initial values: size of lattice, temp and cycles
    read_input(n_spins, mcs, initial_temp, final_temp, temp_step);
    //Initialize the spin matrix
    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    //Random starting point
    idum = -1;

    for (double temp = initial_temp; temp <= final_temp; temp += temp_step)
    {
        // Initialize energy and magnetization
        E = M = 0;
        //Set up array for possible energy changes
        for( int de = -8; de <= 8; de++) w[de+8] = 0;
        for( int de = -8; de <= 8; de+=4) w[de+8] = exp(-de/temp);

        //Initialize array for expectation values
        for( int i = 0; i < 5; i++) average[i] = 0.;
        initialize(E, M, spin_matrix, temp);

        //Start Monte Carlo computation
        for (int cycles = 1; cycles <= mcs; cycles++)
        {
            Metropolis(n_spins, idum, spin_matrix, E, M, w);

            //update expectation values
            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);
        }
        // print results to file
        output(n_spins, mcs, temp, average);
    }

    //WHAT IS THIS?
    free_matrix((void **) spin_matrix); // free memory

    //Close the output file
    ofile.close();
}




//Function that writes data to output file
void spinsystem::output(int m_n_spins, int mcs, double temperature, double *average)
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
    double Mvariance = (M2average - Maverage*Maverage)/(m_n_spins*spins);
    double M2variance = (M2average - Mabsaverage*Mabsaverage)/(m_n_spins*m_n_spins);
    //Why is this different?
    double Mvariance = (M2average - Mabsaverage*Mabsaverage)/(m_n_spins*m_n_spins);

    //Print to an outfile
}













