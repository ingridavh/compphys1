#include <iostream>
#include <armadillo>

#include "jacobi_rotation.h"
#include "exe.h"
#include "schrodinger.h"

using namespace std;
using namespace arma;

void exe()
{
    int N = 4;
    double rho_max = 10.0;

    mat S = generator(N, rho_max);
    cout << S << endl;

    int counter = 0;

    while (eps_test(S) == 0)
    {
        counter++;
        vec max = find_max(S);
        rotate(S, max(1), max(2));
    }

    cout << counter << " rotations were required" << endl;
    vec energies = S.diag();
    energies = sort(energies);

    cout << "The three lowest energy values are " <<
            energies(0) << " , " << energies(1) <<
            " and " << energies(3) << endl;
}
