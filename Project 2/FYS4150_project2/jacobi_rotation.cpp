#include "jacobi_rotation.h"
#include "exe.h"
#include <iostream>
#include <math.h>
#include <armadillo>
using namespace std;
using namespace arma;

int eps_test(mat B, double eps) //default eps-value is set to 1e-8 in header file
{
    //If the off-diagonal elements are small enough,
    //eps_test returns 1, otherwise it returns 0int N,

    vec max = find_max(B);
    int no;

    if (abs(max(0)) >= eps)
    {
        no = 0;
    }
    else
    {
        no = 1;
    }

    return no;
}


//Find the max off-diagonal value of A

vec find_max(mat A)
{

    int N = A.n_cols;

    //Declare output vector
    //max = (matrix element, l, k)
    vec max = zeros<vec>(3);

    for(int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if(i != j)
            {
                if (abs(A(i,j)) >= abs(max(0)))
                {
                    max(0) = A(i,j), max(1) = i, max(2) = j;
                }
            }
        }
    }

    return max;
}

//Find the trigonometric parameters for the rotation

vec find_trig(mat A, int k, int l)
{
    //trig = (cosine, sine)
    vec trig = zeros<vec>(2);
    double tau, t;

    //Calculate trigonometric parameters
    if (A(k,l) != 0.0)
    {
        tau = (A(l,l)- A(k,k))/(float)(2*A(k,l));

        if (tau >= 0)
        {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }

        else
        {
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }

        trig(0) = 1/sqrt(1 + t*t);
        trig(1) = t*trig(0);
    }
    else
    {
        trig(0) = 1.0;
        trig(1) = 0.0;
    }

    return trig;
}

//Rotate the matrix A

void rotate(mat &A, int k, int l)
{
    int N = A.n_cols;

    vec trig = find_trig(A, k, l);

    double s = trig(0);
    double c = trig(1);

    //old elements

    double a_kk = A(k,k); double a_ll = A(l,l);

    //Compute new elements

    A(k,k) = a_kk*c*c - 2.0*A(k,l)*c*s + a_ll*s*s;
    A(l,l) = a_kk*s*s + 2.0*A(k,l)*c*s + a_ll*c*c;

    //Hard-code new non-diagonal l and k
    A(k,l) = 0;
    A(l,k) = 0;

    //Compute new elements

    for ( int i = 0; i < N; i++)
    {
        if (i != k && i != l)
        {
            double a_ik = A(i,k);
            double a_il = A(i,l);

            A(i,k) = a_ik*c - a_il*s;
            A(k,i) = A(i,k);
            A(i,l) = a_il*c + a_ik*s;
            A(l,i) = A(i,l);
        }
    }
}


