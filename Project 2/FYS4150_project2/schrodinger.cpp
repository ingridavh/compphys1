#include "schrodinger.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

mat generator(int N, double rho_max)
{
    //Set step length
    double h = rho_max/(double)(N+1);
    double nondiag_const = -1.0/(h*h);
    double diag_const = 2.0/(h*h);

    vec v(N);
    vec rho(N);
    rho(0) = 0;
    v(0) = rho(0)*rho(0);

    for (int i=0; i<N; i++)
    {
        rho(i) = (i+1)*h;
        v(i) = rho(i)*rho(i);
    }

    mat A = zeros<mat>(N,N);

    //Fill A with values

    A.diag(1).fill(nondiag_const);
    A.diag(-1).fill(nondiag_const);
    A.diag().fill(diag_const);

    for (int i=0; i < N; i++)
    {
        A(i,i) = A(i,i) + v(i);
    }
        return A;
}
