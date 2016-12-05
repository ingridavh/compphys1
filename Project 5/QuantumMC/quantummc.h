#ifndef QUANTUMMC_H
#define QUANTUMMC_H
#include <armadillo>
using namespace arma;

class QuantumMC
{
public:
    QuantumMC(int N);
    void add_perturbation();
    void runMCintegration(int T);
    void no_perturbation();
    double minimize_alpha(double *a_list, int num_steps);
    double minimize_beta(double *b_list, int num_steps);
    void get_low(double *a_list, double *b_list, int num_steps);
    int m_perturbation;
    int m_N;
    double m_alpha;
    int m_mcs;
    int m_blomst;
    double m_beta;
    double m_energy;
    double m_T;
    double m_var;
    double m_omega;

private:
    double wavefunc(const mat &r);
    double wavefunc2(const mat &r);
    double localEnergy(const mat &r);
    int m_nDim;
    double m_steplength;
    double m_h;
    double m_h2;
    double m_r12;
    mat m_rOld;
    mat m_rNew;
};

#endif // QUANTUMMC_H
