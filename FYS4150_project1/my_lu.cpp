#include <iostream>
#include <math.h>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <string>
#include <sstream>
#include "my_lu.h"
using namespace std;
using namespace arma;

double my_lu( int N, bool finderr, bool CPU, bool fileprint, bool CPU_avg)
{
    double x0 = 0;
    double h = 1/double(N+1);

    //Create NxN-matrix with all diagonal elements equal to 2, and vector b
    mat A = 2*eye<mat>(N,N);
    vec b = zeros<vec>(N);

    //Define v_ex with exact values
    vec v_ex = zeros<vec>(N);

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
        b(i) = h*h*100*exp(-10*x);
        v_ex(i) = 1 - (1- exp (-10))*x - exp (-10*x);
    }

    //find LU decomp of A, if needed P is the permutation matrix
    mat L, U;
    lu(L,U,A);

    //Calculate u
    vec w = L.i()*b;
    vec u = U.i()*w;

    //Print CPU time if requested
    finish = clock();
    ((finish-start)/CLOCKS_PER_SEC);
    if (CPU == true)
    {
        cout << "CPU time for LU-decomposition in Armadillo is " << (float(finish-start)/CLOCKS_PER_SEC) << endl;
    }

    //Return CPU time
    if (CPU_avg == true)
    {
        return (float(finish-start)/CLOCKS_PER_SEC);
    }

    //Calculate and print relative error if requested
    if (finderr == true)
    {
        double *eps = new double[N];
        double eps_max = 0;
        for (int i = 0; i < N; i++)
        {
            double diff = (u(i)-v_ex(i))/v_ex(i);
            eps[i] = log10 ( abs (diff));
            if (abs(eps[i]) > abs(eps_max))
            {
                eps_max = eps[i];
            }
        }
        cout << "LU-decomposition: log10(h) is " << log10(h) << " and the error is log10(eps) " << eps_max << endl;

    }
    if (fileprint == true)
    {
        ofstream myfile;
        ostringstream oss;
        oss << "LU_result_N=" << N << ".txt";
        string var = oss.str();
        myfile.open (var, ofstream::out);
        if(!myfile.good()){
            cout << "Dette gikk galt" << endl;
            return 1;
        }

        //For loop that writes results to file
        for (int i = 0; i < (N); i++)
        {
            myfile << u(i) << " " << v_ex(i) << "\n";
        }
        myfile.close();
    }
    return 0;

}
