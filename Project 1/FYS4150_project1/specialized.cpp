#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "gaussian.h"
#include <cmath>

using namespace std;

double specialized (int N, bool finderr, bool CPU, bool fileprint, bool CPU_avg)
{
    //Number of step lengths and h
    double h = 1./(N+1);
    double x_0 = 0;

    //Declare arrays beta, v, f and f_tilde
    double *beta = new double[N+2];
    double *v = new double[N+2];
    double *f = new double[N+2];
    double *f_tilde = new double[N+2];
    double *v_ex = new double[N+2];

    clock_t start, finish;
    start = clock();

    //Fill inn values of beta and calculate x
    for(int i = 1; i < (N+2); i++)
    {
        double x = x_0 + i*h;
        f[i] = h*h*F(x);
        v_ex[i] = v_exact(x);
    }

    //Set boundary conditions
    f_tilde[1] = f[1];
    beta[1] = 2;

    //Step 1
    for(int i = 2; i < (N+2); i++)
    {
        beta[i] = double(i+1)/i;
        f_tilde[i] = f[i] + ((i-1)*f_tilde[i-1])/i;
    }

    //New boundary condition
    v[N] = f_tilde[N]/beta[N];

    //Step 2
    for(int i = N; i > 1; i--)
    {
        v[i-1] = double(i-1)/(i)*(f_tilde[i-1]+v[i]);
    }

    finish = clock();
    ((finish-start)/CLOCKS_PER_SEC);
    if (CPU == true)
    {
        cout << "Specialized: CPU time is " << (float(finish-start)/CLOCKS_PER_SEC) << endl;
    }
    //returns CPU time value
    if (CPU_avg == true)
    {
        return (float(finish-start)/CLOCKS_PER_SEC);
    }




    //Write results to a file spec_result_N= (n-value).txt
    if (fileprint == true)
    {
        ofstream myfile;
        ostringstream oss;
        oss << "spec_result_N=" << N << ".txt";
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

    //Compute and print relative error if requested
    if (finderr == true)
    {
        double *eps = new double[N+2];
        double eps_max = 0;
        for (int i=1; i < (N+1); i++)
        {
            double diff = (v[i] - v_ex[i])/v_ex[i];
            eps[i] = log10 (abs (diff));
            //Find maximum error
            if (abs(eps[i]) > abs(eps_max))
            {
                eps_max = eps[i];
            }
        }
        cout << "Specialized: log10(h) is " << log10(h) << " and the error is log10(eps) " << eps_max << endl;
    }


    return 0;
}
