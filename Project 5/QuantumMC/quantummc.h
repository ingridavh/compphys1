#ifndef QUANTUMMC_H
#define QUANTUMMC_H


class QuantumMC
{
public:
    QuantumMC(int N);
    void perturbation(int M);
    void runMCintegration();
    int m_perturbation;
    int m_N;
    double m_alpha;

private:
    double wavefunc(const mat &r);
    double localEnergy(const mat &r);
    int m_nDim;
    int m_Q;
    double m_steplength;
    double m_h;
    double m_h2;
    long m_idum;
    int m_mcs;
    mat m_rOld;
    mat m_rNew;
};

#endif // QUANTUMMC_H
