#include <iostream>
#include <armadillo>

#include "jacobi_rotation.h"
#include "exe.h"

using namespace std;
using namespace arma;

void exe(mat A)
{
    int counter = 0;

    vec max;
    while (eps_test(A) == 0)
    {
        counter++;
        max = find_max(A);
        rotate(A, max(1), max(2));
    }

    cout << counter << " rotations were required" << endl;
    vec energies = A.diag();
    energies = sort(energies);

    cout << "The three lowest energy values are " <<
            energies(0) << " , " << energies(1) <<
            " and " << energies(3) << endl;
}
