#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
#include "time.h"
#include <string>
#include <sstream>

using namespace std;

//Function for f-values
double F (double x)
{
    double y = -10*x;
    double F_value = 100* exp (y);
    return F_value;
}

//Function for exact v-values
double v_exact (double x)
{
    double v = 1 - (1 - exp (-10))*x - exp(-10*x);
    return v;
}


//Step 1: decomposition function
double decomp (double a, double b, double c_old, double beta_old)
{
    double beta = b - (a*c_old)/beta_old;
    return beta;
}

//Step 1: Forward substitution function
double forward (double a, double beta_old, double f, double f_tilde_old)
{
    double f_tilde = f - (a*f_tilde_old)/beta_old;
    return f_tilde;
}

//Step 2: backward substitution function
double backward(double f_tilde_old, double c_old, double u_new, double beta_old)
{
    double u = (f_tilde_old - c_old*u_new)/beta_old;
    return u;
}


//MAIN FUNCTION

//Give N as function argument
//Default is error=false, CPU=false, fileprint=false
//If set to true, the function prints error and CPU and writes result to file

double gaussian(int N, bool finderr, bool CPU, bool fileprint, bool CPU_avg)
{
    //Number of step lengths and h
    double h = 1./(N+1);
    double x_0 = 0;

    //Declare matrix A, A_new and arrays v, f and f_tilde
    double *a = new double[N+2];
    double *b = new double[N+2];
    double *c = new double[N+2];
    double *beta = new double[N+2];
    double *v = new double[N+2];
    double *f = new double[N+2];
    double *f_tilde = new double[N+2];
    double *v_ex = new double[N+2];

    clock_t start, finish;
    start = clock();
    //Fill inn values of A
    for(int i = 0; i < (N+2); i++)
    {
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
        double x = x_0 + i*h;
        f[i] = h*h*F(x);
        v_ex[i] = v_exact(x);
    }

    //Set boundary conditions
    beta[1] = b[1];
    f_tilde[1] = f[1];
    v[0] = 0;
    v[N+1] = 0;

    //Step 1
    for (int i = 2; i < (N+1); i++)
    {
        beta[i] = decomp(a[i], b[i], c[i-1], beta[i-1]);
        f_tilde[i] = forward(a[i], beta[i-1],f[i], f_tilde[i-1]);
    }

    v[N] = f_tilde[N]/beta[N];

    //Step 2
    for(int i = N; i > 1; i--)
    {
        v[i-1] = backward(f_tilde[i-1], c[i-1], v[i], beta[i-1]);
    }

    finish = clock();
    ((finish-start)/CLOCKS_PER_SEC);

    //Print CPU time if requested

    if (CPU == true)
    {
        cout << "Gaussian: CPU time is " << (float(finish-start)/CLOCKS_PER_SEC) << endl;
    }

    //Return CPU time
    if (CPU_avg == true)
    {
        return (float(finish-start)/CLOCKS_PER_SEC);
    }


    //Calculate and print relative error if requested

    if (finderr == true)
    {
        double *eps = new double[N+2];
        double eps_max = 0;
        for (int i = 0; i < N+2; i++)
        {
            double diff = (v[i]-v_ex[i])/v_ex[i];
            eps[i] = log10 ( abs (diff));
            if (abs(eps[i]) > abs(eps_max))
            {
                eps_max = eps[i];
            }
        }
        cout << "Gaussian: log10(h) is " << log10(h) << " and the error is log10(eps) " << eps_max << endl;

    }

    //Write results to a file Gauss_result_N (n-value).txt if requested


    if (fileprint == true)
    {
        ofstream myfile;
        ostringstream oss;
        oss << "Gauss_result_N=" << N << ".txt";
        string var = oss.str();
        myfile.open (var, ofstream::out);
        if(!myfile.good()){
            cout << "Dette gikk galt" << endl;
            return 1;
        }

        //For loop that writes results to file
        for (int i = 0; i < (N+2); i++)
        {
            myfile << v[i] << " " << v_ex[i] << "\n";
        }
        myfile.close();
    }

    return 0;
}
