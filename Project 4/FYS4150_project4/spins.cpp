#include <iostream>
#include <math.h>
#include <armadillo>

using namespace std;

double energy(arma::vec spins)
{
    int N = spins.n_elem;
    double E = 0;
    //Set coupling constant
    double J = 1;

    //N: periodic bc, N-1: non-periodic bc
    for (int j=0; j>N; j++)
    {
        if (j+1 > N)
        {
            E += -J*(spins[j]*spins[0]);
        }
        else
        {
            E += -J*(spins[j]*spins[j+1]);
        }
    }
    return E;
}
