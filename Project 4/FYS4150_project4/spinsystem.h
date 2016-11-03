#ifndef SPINSYSTEM_H
#define SPINSYSTEM_H
#include <armadillo>

using namespace std;


class spinsystem
{
public:
    spinsystem(int n_spins);
    double CalculateEnergy();
    double CalculateMagnetization();
    void Metropolis(long& idum, double& E, double &M, double *w);
    double periodic(int i, int limit, int add);
    void go(string outfilename, int mcs, double initial_temp, double final_temp, double temp_step);
    void output(int mcs, double temperature, double *average);
private:
    void initialize(double& E, double& M, double temp);
    void read_input (int&, int&, double&, double&, double&);
    double m_J;
    arma::mat m_spin_matrix;
    int m_n_spins;
    ofstream m_ofile;
};

#endif // SPINSYSTEM_H
