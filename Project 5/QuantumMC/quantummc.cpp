#include "quantummc.h"
#include <lib.h>
#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

QuantumMC::QuantumMC(int N) :
    m_nDim(3),
    m_Q(2),
    m_N(N),
    m_perturbation(0),
    m_steplength(1.0),
    m_h(0.001),
    m_h2(1000000),
    m_idum(-1),
    m_alpha(0.5*m_Q),
    m_mcs(1000000)
{        
}

void QuantumMC::perturbation(int M){
    m_perturbation = 1;
}

void QuantumMC::energy(){
//If m_perturbation is 1 add the energy from the perturbation
}

void QuantumMC::runMCintegration(){
    rOld = zeros<mat>(m_N, m_nDim);
    rNew = zeros<mat>(m_N, m_nDim);

    //Set all functions and energies to zero
    double wavefuncOld = 0;
    double wavefuncNew = 0;
    double ESum = 0;
    double ESum2 = 0;
    double dE;

    //Initial trial positions
    for(int i = 0; i < m_N; i++){
        for(int j=0; j<m_nDim; j++){
            rNew(i,j) = rOld(i,j) + m_steplength*(ran2(&idum) -0.5);
        }
    }

    rNew = rOld;

    //Loop over Monte Carlo cycles
    for(int cycle=0; cycle < m_mcs; cycle++)
    {
        //Store the current value of the wave function
        wavefuncOld = wavefunc(rOld);

        //New position for testing
        for(int i= 0; i < m_nDim; i++){
            for(int j=0; j < m_N; i++){
                rNew(i, j) = rOld(i, j) + m_steplength*(ran2(&idum) - 0.5);
            }

            //Recalculate the value of the wave function
            wavefuncNew = wavefunc(rNew);

            //Check for step acceptance (yes= update position, no=discard new position)
            if(ran2(&idum) <= )
        }

    }

}






















