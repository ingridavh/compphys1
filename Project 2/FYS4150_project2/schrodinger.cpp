#include "schrodinger.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

mat generator(int N, double rho_max, double rho_min)
{
    //Set step length
    double h = rho_max/(double)N;

    mat A = zeros<mat>(N,N);

    //Fill A with values
    for (int i=0; i < N-1; i++)
    {
        double rho = rho_min + i*h;
        A(i,i) = 2.0/(h*h) + rho*rho;
        A(i, i+1) = -1.0/(h*h);
        A(i+1, i) = -1.0/(h*h);
    }

    //Set final diagonal matrix element
    double rho_N = (N-1)*(double)h;
    A(N-1,N-1) = 2./(h*h)+rho_N*rho_N;
    return A;
}
