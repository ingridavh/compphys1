#ifndef SPINSYSTEM_H
#define SPINSYSTEM_H
#include <armadillo>

using namespace std;


class spinsystem
{
public:
    spinsystem(int n_spins);
    void Metropolis(double& E, double &M, double *w, double temp);
    double periodic(int i, int limit, int add);
    void go(string outfilename, int mcs, double initial_temp, double final_temp, double temp_step, int rankProcess, int NProcesses);
    void output(int mcs, double temperature, arma::vec avg);
    arma::vec getExpecs();
private:
    void initialize(double& E, double& M, double temp);
    void read_input (int&, int&, double&, double&, double&);
    double m_J;
    arma::mat m_spin_matrix;
    int m_n_spins;
    double m_average[5];
    double m_finaltemp;
    double m_no_accept;
    int m_mcs;
    ofstream m_ofile;
};

#endif // SPINSYSTEM_H
