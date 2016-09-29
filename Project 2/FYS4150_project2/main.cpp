#include <iostream>
#include "jacobi_rotation.h"
#include "schrodinger.h"
#include "exe.h"
#include <armadillo>
#include <math.h>
using namespace std;
using namespace arma;

//int N = 30;

int main()
{
    //mat S = generator(N, 10);

    exe();

//    int counter = 0;
//    //cout << S << endl;

//    mat B;
//    while (eps_test(N, S) == 0)
//    {
//        counter++;
//        B = Jacobi(S);
//        S = B;
//    }
//    cout << "After " << counter
//         << " rotations "
//         << "the eigenvalues are " << endl;
//    for (int i = 0; i<N; i++)
//    {
//         cout << B(i,i) << endl;
//    }

    return 0;
}

