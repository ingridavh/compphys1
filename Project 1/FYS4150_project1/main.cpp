#include <iostream>
#include <math.h>
#include <fstream>
#include "gaussian.h"
#include "specialized.h"
#include "time.h"
#include "my_lu.h"

using namespace std;

//All funtions have default values (int N=10, finderr = false, CPU = false, fileprint = false, CPU_avg = false)
//finderr=true: calculates and prints the relative error
//fileprint=true: writes calculated and exact values to file named (Name of method)_result_N=(n-value).txt
//CPU = true: prints CPU time for operation
//CPU_avg = true: returns CPU time for operation


int main()
{
    int N = 100;
    //(N, finderr, CPU, fileprint), e.g. (N, true, true)
    gaussian (N, true);
    specialized(N, true);
    my_lu(N, true);

    //Compute average CPU times
//    int n = 100;
//    float t1 = 0;
//    float t2 = 0;
//    float t3 = 0;
//    for (int i = 0; i <= n; i++)
//    {
//        t1 = t1 + specialized(N, false, false, false, true);
//        t2 = t2 + gaussian(N, false, false, false, true);
//        t3 = t3 + my_lu(N, false, false, false, true);
//    }
//    cout << "Average CPU time for Specialized is " << t1/n << endl;
//    cout << "Average CPU time for Gaussian is " << t2/n << endl;
//    cout << "Average CPU time for LU-decomposition is " << t3/n << endl;
//    return 0;
}

