#include "quantummc.h"
#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

QuantumMC::QuantumMC(int N) :
    m_nDim(3),
    m_omega(1),
    m_N(N),
    m_perturbation(0),
    m_steplength(1.0),
    m_h(0.001),
    m_h2(1/double(m_h*m_h)),
    m_alpha(1),
    m_mcs(1000000)
{        
}

void QuantumMC::add_perturbation(){
    m_perturbation = 1;
}

void QuantumMC::runMCintegration(){
    //Initialiser med forskjellige tall hver gang
    srand(time(NULL));
    double faktor = 1.0/RAND_MAX;

    mat rOld = zeros<mat>(m_N, m_nDim);
    mat rNew = zeros<mat>(m_N, m_nDim);

    //Set all functions and energies to zero
    double wavefuncOld = 0;
    double wavefuncNew = 0;
    double ESum = 0;
    double ESum2 = 0;

    //Initial trial positions
    for(int i = 0; i < m_N; i++){
        for(int j = 0; j < m_nDim; j++){
            double blomst = double (rand())*faktor;
            cout << "The random number is " << blomst << endl;
            rOld(i,j) = m_steplength * (blomst - 0.5);
        }
    }

    rNew = rOld;

    //Loop over Monte Carlo cycles
    for(int cycle=0; cycle < m_mcs; cycle++)
    {
        //Store the current value of the wave function
        wavefuncOld = wavefunc(rOld);

        //New position for testing
        for(int i = 0; i < m_N; i++){
            for(int j = 0; j < m_nDim; j++){
                rNew(i,j) = rOld(i,j) + m_steplength*(double(rand())*faktor - 0.5);
            }

            //Recalculate the value of the wave function
            wavefuncNew = wavefunc(rNew);

            //Check for step acceptance (yes=update position, no=discard new position)
            if(double (rand())*faktor <= (wavefuncNew*wavefuncNew)/(wavefuncOld*wavefuncOld)) {
                for (int j = 0; j < m_nDim; j++) {
                    rOld(i,j) = rNew(i,j);
                    wavefuncOld = wavefuncNew;
                }

            } else {
                for (int j= 0; j < m_nDim; j++) {
                    rNew(i,j) = rOld(i,j);
                }
            }

            //Update energy values
            double dE = localEnergy(rNew);
            ESum += dE;
            ESum2 += dE*dE;
        }
    }

    double energy = ESum/double (m_mcs*m_N);
    double energy2 = ESum2/double (m_mcs*m_N);
    cout << "Energy: " << energy << " Energy (squared sum): " << energy2 << endl;

}


double QuantumMC::wavefunc(const mat &r){
    double argument = 0;
    for(int i = 0; i < m_N; i++){
        double r_oneparticle = 0;
        for(int j = 0; j < m_nDim; j++){
            r_oneparticle += r(i,j) * r(i, j);
        }

        argument += sqrt(r_oneparticle);
    }

    return exp(-argument * m_alpha);
}


//Compute the energy for a certain configuration of particle positions
double QuantumMC::localEnergy(const mat &r){
    mat rPlus = zeros<mat>(m_N, m_nDim);
    mat rMinus = zeros<mat>(m_N, m_nDim);
    rPlus = rMinus = r;

    double wavefuncminus = 0;
    double wavefuncplus = 0;
    double wavefunccurrent = wavefunc(r);

    //Calculate kinetic energy using brute force

    double EK = 0;

    for(int i = 0; i < m_N; i++){
        for(int j = 0; j < m_nDim; j++){
            rPlus(i,j) += m_h;
            rMinus(i,j) -= m_h;
            wavefuncminus = wavefunc(rMinus);
            wavefuncplus = wavefunc(rPlus);
            EK -= (wavefuncminus + wavefuncplus - 2 * wavefunccurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }

    //EK = 0.5 * m_h2 * EK / wavefunccurrent;
    EK = 0.5 * m_h2 * EK;



    //Calculate potential energy
    double EP = 0;
    double r_oneparticle = 0;

    for(int i = 0; i < m_N; i++){
        r_oneparticle = 0;
        for(int j = 0; j < m_nDim; j++){
            r_oneparticle += r(i,j) * r(i,j);
        }
        EP -= m_omega / sqrt(r_oneparticle);
    }

    //Contribution from electron-electron potential
    if (m_perturbation = 1){
        double r12 = 0;
        for(int i = 0; i < m_N; i++){
            r_oneparticle = 0;
            for(int j= i+1; j < m_N; j++){
                r12 = 0;
                for(int k = 0; k < m_nDim; k++){
                    r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                EP += 1/sqrt(r12);
            }
        }
    }

    return EK+EP;
}
