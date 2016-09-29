#include <iostream>
#include "jacobi_rotation.h"
#include "schrodinger.h"
#include <armadillo>
#include <math.h>
using namespace std;
using namespace arma;


int main()
{
    int N = 100;
    double rho_max = 10;

    mat R = eye<mat>(N, N);
    mat A = generator(N, rho_max);

    int counter = 0;

    while (eps_test(A) == 0)
    {
        counter++;
        vec max = find_max(A);
        rotate(A, R, max(1), max(2));
    }

    cout << counter << " rotations were required" << endl;
    //cout << A << endl;

    vec energies = A.diag();
    energies = sort(energies);

    cout << "The three lowest energy values are " <<
            energies(0) << " , " << energies(1) <<
            " and " << energies(2) << endl;
    return 0;
}

