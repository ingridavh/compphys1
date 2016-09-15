#include <iostream>
#include <math.h>
#include <armadillo>
#include "my_lu.h"
using namespace std;
using namespace arma;

void my_lu()
{
    int N = 10;
    double x0 = 0;
    double h = 1/double(N+1);
    //Create NxN-matrix with all diagonal elements equal to 2, and vector b
    mat A = 2*eye<mat>(N,N);
    vec b = zeros<vec>(N);
    //Define a vector x with the unknowns
    vec x = zeros<vec>(N);

    clock_t start, finish;
    start = clock();

    //Set elements directly under and over diagonal
    for (int i=0; i<N-1; i++)
    {
        A(i,i+1)=-1;
        A(i+1,i)=-1;
    }
    for (int i=0; i<N; i++)
    {
        double x = x0 + (i+1)*h;
        b(i) = 100*exp(-10*x);
    }

    //find LU decomp of A, if needed P is the permutation matrix
    mat L, U;
    lu(L,U,A);

    //Check that A = LU
    (A-L*U).print("Test of LU decomposition");

    //Calculate first step with inverse of L
    vec w = L.i()*b;
    vec u = U.i()*w;
    u.print("The unknown function u=");

    finish = clock();
    ((finish-start)/CLOCKS_PER_SEC);
    cout << "CPU time for LU-decomposition in Armadillo is " << (float(finish-start)/CLOCKS_PER_SEC) << endl;
    //Calculate second step with inverse of U

}
