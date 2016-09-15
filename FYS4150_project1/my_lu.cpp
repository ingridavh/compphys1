#include <iostream>
#include <armadillo>
#include "my_lu.h"
using namespace std;
using namespace arma;

void my_lu()
{
    int N = 10;
    mat A = 2*eye<mat>(N,N);
    vec b = randu<vec>(N);

    A.print("A=");
    b.print("b=");

    //solve Ax=b
    vec x = solve(A,b);

    //print x
    x.print("x=");

    //find LU decomp of A, if needed P is the permutation matrix
    mat L, U;
    lu(L,U,A);

    //print L
    L.print("L=");
    //print U
    U.print("U=");

    //Check that A = LU
    (A-L*U).print("Test of LU decomposition");

//    return 0;
}
