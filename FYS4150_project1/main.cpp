#include <iostream>
#include <math.h>

using namespace std;

//Create n and arrays with n slots

int *n = new int;

double *a = new double[*n];
double *b = new double[*n];
double *c = new double[*n];
double *bt = new double[*n];

int main()
{
    int N = 10;
    *n = N;
    a[0] = 0;
    b[0] = 2;
    c[0] = -1;
    for(int i = 1; i <= *n; i++)
    {
        *a[i] = -1;
        *b[i] = 2;
        *c[i] = -1;
    }

    return 0;
}

