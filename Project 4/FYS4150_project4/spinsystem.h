#ifndef SPINSYSTEM_H
#define SPINSYSTEM_H
#include <armadillo>


class spinsystem
{
public:
    spinsystem();
    double CalculateEnergy();
    double CalculateMagnetization();
    double Metropolis(long& idum, int **m_spin_matrix, double& E, double &M, double *w);
    double periodic(int i, int limit, int add);
private:
    double m_J;
    arma::mat m_spin_matrix;
    int m_n_spins;
};

#endif // SPINSYSTEM_H
