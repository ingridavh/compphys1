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
    m_h(0.00001),
    m_h2(1./double(m_h*m_h)),
    m_alpha(1),
    m_mcs(2500000),
    m_blomst(0),
    m_beta(1.0),
    m_r12(0)
{        
}

void QuantumMC::add_perturbation(){
    m_perturbation = 1;
}

void QuantumMC::no_perturbation(){
    m_perturbation = 0;
}

void QuantumMC::runMCintegration(int T){
    //Takes the number T as input. T=1 means that the first test function should be used
    //T != 1 means that the second test function should be used

    m_T = T;
    m_steplength = 1.0;
    m_energy = 0;
    m_blomst = 0;
    m_r12 = 0;

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
            rOld(i,j) = m_steplength * (blomst - 0.5);
        }
    }

    rNew = rOld;
    double r12_avg=0;

    //Loop over Monte Carlo cycles
    for(int cycle=0; cycle < m_mcs; cycle++)
    {
        if (cycle % 1000 == 0)
        {
            if (m_blomst/double(cycle*m_N) < 0.4 ){
                m_steplength -= 1;
            } else if (m_blomst/double(cycle*m_N) > 0.6 ) {
                m_steplength += 1;
            }
        }

        //Store the current value of the wave function
        wavefuncOld = (m_T == 1) ? wavefunc(rOld) : wavefunc2(rOld);


        //New position for testing
        for(int i = 0; i < m_N; i++){
            for(int j = 0; j < m_nDim; j++){
                rNew(i,j) = rOld(i,j) + m_steplength*(double(rand())*faktor - 0.5);
            }

            //Recalculate the value of the wave function
            wavefuncNew = (m_T == 1) ? wavefunc(rNew) : wavefunc2(rNew);


            //Check for step acceptance (yes=update position, no=discard new position)
            if(double (rand())*faktor <= (wavefuncNew*wavefuncNew)/(wavefuncOld*wavefuncOld)) {
                for (int j = 0; j < m_nDim; j++) {
                    rOld(i,j) = rNew(i,j);
                }
                wavefuncOld = wavefuncNew;
                m_blomst += 1;

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

    double energy = ESum / double (m_mcs*m_N);
    double energy2 = ESum2 / double (m_mcs*m_N);
    m_energy = energy;
    m_var = energy2 - energy*energy;
    m_r12 /= (double (m_mcs*m_N));
//    cout << "Energy: " << energy << " Energy (squared sum): " << energy2 << endl;
//    cout << "Variance: " << energy2 - energy*energy << endl;
//    cout << "Percentage of accepted configurations " << m_blomst/double(m_mcs*m_N)*100 << "%" << endl;
    cout << "For omega: " << m_omega << ", the average distance r12 is " << m_r12 << endl;
}

double QuantumMC::wavefunc(const mat &r){
    double argument = 0;
    for(int i = 0; i < m_N; i++){
        double r_oneparticle = 0;
        for(int j = 0; j < m_nDim; j++){
            r_oneparticle += r(i,j) * r(i, j);
        }

        argument += r_oneparticle;
    }
    return exp(-0.5 * argument * m_alpha * m_omega );
}

double QuantumMC::wavefunc2(const mat &r){
    double argument = 0;
    for(int i=0; i < m_N; i++){
        double r_oneparticle = 0;
        for(int j=0; j < m_nDim; j++){
            r_oneparticle += r(i,j)*r(i,j);
        }
        argument += r_oneparticle;
    }
    double argument12 = 0;
    double r12 = 0;

    for(int i = 0; i< m_N; i++){
        for(int j=1+i; j < m_N; j++){
            r12 = 0;
            for(int k = 0; k < m_nDim; k++){
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            argument12 += r12;
        }
    }
    return exp( - 0.5 * m_alpha * m_omega * argument)*exp(argument12* m_beta/(2*(1+m_beta*argument12)));
}

//Compute the energy for a certain configuration of particle positions
double QuantumMC::localEnergy(const mat &r){
    mat rPlus = zeros<mat>(m_N, m_nDim);
    mat rMinus = zeros<mat>(m_N, m_nDim);
    rPlus = rMinus = r;

    double wavefuncminus = 0;
    double wavefuncplus = 0;

    //Compact if-test
    double wavefunccurrent = (m_T == 1) ? wavefunc(r) : wavefunc2(r);

    //Calculate kinetic energy using brute force
    // -1/2 nabla^2 terms

    double EK = 0;

    for(int i = 0; i < m_N; i++){
        for(int j = 0; j < m_nDim; j++){
            rPlus(i,j) += m_h;
            rMinus(i,j) -= m_h;

            if (m_T == 1){
                wavefuncminus = wavefunc(rMinus);
                wavefuncplus = wavefunc(rPlus);
            } else {
                wavefuncminus = wavefunc2(rMinus);
                wavefuncplus = wavefunc2(rPlus);
            }
            EK -= (wavefuncminus + wavefuncplus - 2 * wavefunccurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }

    //EK = 0.5 * m_h2 * EK / wavefunccurrent;
    EK *= 0.5 * m_h2 / wavefunccurrent ;

    //Calculate potential energy
    // 1/2 omega^2 r^2 terms
    double EP = 0;
    double r_oneparticle = 0;

    for(int i = 0; i < m_N; i++){
        r_oneparticle = 0;
        for(int j = 0; j < m_nDim; j++){
            r_oneparticle += r(i,j) * r(i,j);
        }
        EP += r_oneparticle;
    }
    EP *= 0.5 * m_omega * m_omega;

    //Contribution from electron-electron potential
    // 1/r_12 terms
    if (m_perturbation == 1){
        double r12 = 0;
        for(int i = 0; i < m_N; i++){
            r_oneparticle = 0;
            for(int j= i+1; j < m_N; j++){
                r12 = 0;
                for(int k = 0; k < m_nDim; k++){
                    r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                EP += 1/sqrt(r12);
                m_r12 += sqrt(r12);
            }
        }
    } else {
        EP += 0;
    }
    return EK+EP;
}

double QuantumMC::minimize_alpha(double *a_list, int num_steps)
{
    //Find the alpha value that minimizes the energy
    double E_min = 10;
    double a_min = 0;
    double deltas [num_steps];

    for(int i =0; i < num_steps; i++){
        cout << a_list[i] << " , ";
        m_alpha = a_list[i];
        runMCintegration(m_T);
        deltas[i] = m_steplength;
        if(m_energy < E_min){
            E_min = m_energy;
            a_min = a_list[i];
        }
    }
    cout << endl;
    m_energy = E_min;
    cout << "The alpha value " << a_min << " minimizes the energy. Energy: " << m_energy << endl;
    return a_min;
}

double QuantumMC::minimize_beta(double *b_list, int num_steps)
{
    //Find the beta value that minimizes the energy

    //Give error if we are using the first wavefunction
    if (m_T ==1){
        cout << "Error! The wavefunction you are using has no beta values!" << endl;
        terminate();
    }

    double E_min = 10;
    double b_min = 0;
    for(int i = 0; i < num_steps; i++){
        cout << b_list[i] << " , ";
        m_beta = b_list[i];
        runMCintegration(m_T);
        if(m_energy < E_min){
            E_min = m_energy;
            b_min = b_list[i];
        }
    }

    m_energy = E_min;
    cout << endl;
//    cout << "The beta value " << b_min << " minimizes the energy. Energy: " << E_min << endl;
    return b_min;
}

//Function that find alpha and beta to get the best fit
void QuantumMC::get_low(double *a_list, double *b_list, int num_steps){
    if (m_T == 1){
        cout << "Error! The wavefunction you are using has no beta!" << endl;
        terminate();
    }

    cout << "Initial alpha: " << m_alpha << ", initial beta: " << m_beta << endl;

    int runs = 1;
    double E_old; double E_new; double alpha_old; double beta_old;

    //Find lowest energy for alphas, with constant initial beta
    m_alpha = minimize_alpha(a_list, num_steps);

    E_old = m_energy;
    alpha_old = m_alpha;
    beta_old = m_beta;

    //Find lowest energy for betas, with alpha minimized
    m_beta = minimize_beta(b_list, num_steps);
    E_new = m_energy;

    if (E_new > E_old)
        m_beta = beta_old;

    //If minimizing beta gives lower energy, re-evaluate alpha and beta
    while (E_new < E_old){
        runs++;
        E_old = E_new;

        //Minimize alpha for new beta value, and store temporary energy
        m_alpha = minimize_alpha(a_list, num_steps);
        double E_temp = m_energy;

        //If we found a better alfa, minimize beta again
        if (E_temp < E_old){
            E_old = E_temp;
            m_beta = minimize_beta(b_list, num_steps);
            E_new = m_energy;
            if (E_new > E_old)
                E_new = E_temp;
        } else {
            //Otherwise, old energy was lowest, jump out of loop
            E_new = E_old;
        }
    }

    cout << "For omega " << m_omega << ",the best fit comes from alpha= " << m_alpha << ", and beta= "
         << m_beta << " after " << runs << " runs. The energy is " << E_old << endl;
}
