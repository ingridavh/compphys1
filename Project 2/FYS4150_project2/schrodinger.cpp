#include "schrodinger.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

mat generator(int N, double rho_max, double omega)
{
    //Set step length
    double h = rho_max/(double)N;
    double nondiag_const = -1.0/(h*h);
    double diag_const = 2.0/(h*h);

    vec v(N);
    vec rho(N);
    rho(0) = 0;
    v(0) = rho(0)*rho(0);

    //Calculate the potenial, if the frequency is
    //not specified in the function call, it is set to
    //default for non-interacting

    for (int i=0; i<N; i++)
    {
        rho(i) = (i+1)*h;
        if (omega == 0.0)
        {
            v(i) = rho(i)*rho(i);
        }
        else
        {
            v(i)=omega*omega*rho(i)*rho(i) + 1.0/rho(i);
        }


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
