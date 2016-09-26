#include <iostream>
#include "jacobi_rotation.h"
#include <armadillo>
#include <math.h>
using namespace std;
using namespace arma;

int N = 2;
mat A = {
    {2,1},
    {1,2}};

int main()
{
    int counter = 0;
    cout << "A = " << A << endl;
    mat B;
    while (eps_test(N, A) == 0)
    {
        counter++;
        B = Jacobi(N,A);
        cout << "B = " << B << endl;
        A = B;
    }
    cout << "after " << counter << " S transformations, B is " << B << endl;
    return 0;
}

