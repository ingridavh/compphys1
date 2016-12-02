#ifndef QUANTUMMC_H
#define QUANTUMMC_H
#include <armadillo>
using namespace arma;

class QuantumMC
{
public:
    QuantumMC(int N);
    void add_perturbation();
    void runMCintegration();
    int m_perturbation;
    int m_N;
    double m_alpha;
    int m_mcs;

private:
    double wavefunc(const mat &r);
    double localEnergy(const mat &r);
    int m_nDim;
    int m_omega;
    double m_steplength;
    double m_h;
    double m_h2;
    mat m_rOld;
    mat m_rNew;
};

#endif // QUANTUMMC_H
