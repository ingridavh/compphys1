#include "jacobi_rotation.h"
#include <iostream>
#include <math.h>
#include <armadillo>
using namespace std;
using namespace arma;

//mat jacobi_tol (int N, mat B, double eps)
//{
//    int epsi = 0;
//    int counter = 0;
//    mat Bnew = Jacobi(N, B);

//    while(epsi == 0)
//    {
//        counter++;
//        Bnew = Jacobi(N, Bnew);
//        //Test if off-diagonal elements are small enough
//        epsi = eps_test(N, Bnew, eps);
//        cout << "epsi is " << epsi << endl;
//        cout << Bnew << endl;
//    }
//    cout << counter << " S matrices were needed" << endl;
//    return Bnew;
//}

mat Jacobi (int N, mat A) //Default eps = 1e-8 in header file

//The Jacobi function takes an NxN-matrix A and finds its eigenvalues,
//which are returned as the diagonal elements of the matrix B

{
    //Declare variables
    mat B;
    //find maximum off-diagonal element
    vec ut_max = find_max(A, N);
    int k = (int) ut_max(1);
    int l = (int) ut_max(2);
    cout << "k is " << k << " and l is " << l << " and max value is " << ut_max(0) << endl;

    //calculate trigonometric parameters
    vec ut_trig = find_trig(A,k,l);
    double c = ut_trig(0);
    double s = ut_trig(1);
    cout << "cos is " << c << "sin is " << s << endl;

    //calculate new matrix B
    mat B_ = B_new(A, N, k, l, c, s);
    cout << "B = " << B_ << endl;
    return B_;
}


int eps_test(int N, mat B, double eps) //default eps-value is set to 1e-8 in header file
{
    //If the off-diagonal elements are small enough,
    //eps_test returns 1, otherwise it returns 0
    double max_val = 0;
    int no;
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if(j != i){
                if (abs(B(i,j)) >= abs(max_val)){
                    max_val = B(i,j);
                }
            }
        }
    }
    if (abs(max_val) >= eps){
        no = 0;
    }
    else{
        no = 1;
    }
    return no;
}


vec find_max(mat A,int N)
{
    //Find the max off-diagonal value of A
    double max_val = 0;
    int k, l;
    //Declare output vector
    vec out = zeros<vec>(3);



    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if(i != j){
                //cout << "i er " << i << ", j er " << j << " og B(i,j) er " << B(i,j) << endl;
                if (abs(A(i,j)) >= abs(max_val)){
                    max_val = A(i,j);
                    k = i;
                    l = j;
                }
            }
        }
    }
    //out = (matrix element, l, k)
    out(0) = max_val;
    out(1) = k;
    out(2) = l;
    return out;
}

vec find_trig(mat A, int k, int l)
{
    vec out = zeros<vec>(2);
    double tau, t, s, c;
    //Calculate trigonometric parameters
    cout << "A(l,l) is " << A(l,l) << "A(k,k) is " << A(k,k) << endl;
    tau = (A(l,l)- A(k,k))/(float)(2*A(k,l));
    cout << "tau is " << tau << endl;
    if (tau > 0)
    {
        t = 1.0/(tau + sqrt(1.0 + tau*tau));
    }
    else
    {
        t = -1.0/(-tau + sqrt(1.0 + tau*tau));
    }
    c = 1/(1 + t*t);
    s = t*c;
    out(0) = c;
    out(1) = s;
    return out;
}

mat B_new(mat A, int N, int k, int l, double c, double s)
{
    mat B = zeros<mat>(N,N);

    //Compute the elements in B
    B(k,k) = A(k,k)*c*c - 2.0*A(k,l)*c*s + A(l,l)*s*s;
    B(l,l) = A(l,l)*c*c + 2.0*A(k,l)*c*s + A(k,k)*s*s;
    B(k,l) = 0.0; //hardcode the zero-elements
    B(l,k) = 0.0;

    for (int i=0; i < N; i++)
    {
        if (i != k && i != l)
        {
            B(i,i) = A(i,i);
            B(i,k) = A(i,k)*c - A(i,l)*s;
            B(k,i) = B(i,k);
            B(i,l) = A(i,l)*c + A(i,k)*s;
            B(l,i) = B(i,l);
        }
    }

    return B;
}


