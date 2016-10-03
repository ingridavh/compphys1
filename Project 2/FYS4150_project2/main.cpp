#include <iostream>
#include "jacobi_rotation.h"
#include "schrodinger.h"
#include <armadillo>
#include <math.h>

#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using namespace arma;


int main()
{
    int N = 100;
    double rho_max = 4;
    vec omegas = {0.01 , 0.5 , 1 , 5};



    //Calculate energies for several omegas

    for (int i = 0; i < 4;i++)
    {
        mat R = eye<mat>(N, N);
        double omega = omegas(i);
        int counter = 0;

        mat A = generator(N, rho_max, omega);

        while (eps_test(A) == 0)
        {
            counter++;
            vec max = find_max(A);
            rotate(A, R, max(1), max(2));
        }

        vec energies = A.diag();

        //Find indices of the 3 lowest elements
        uvec ind = sort_index(energies);

        cout << "For omega = " << omega << " the ground state energy is "
             << energies(ind(0)) << endl;


        //Write results to file
        ofstream myfile;
        ostringstream oss;
        oss << "eigenvec_omega=" << omega << ".txt";
        string var = oss.str();
        myfile.open (var, ofstream::out);
        if(!myfile.good()){
            cout << "Dette gikk galt" << endl;
            return 1;
        }

        //For loop that writes results to file

        for (int i = 0; i < N; i++)
        {
            double r0 = R(i, ind(0))*R(i, ind(0));
            double r1 = R(i, ind(1))*R(i, ind(1));
            double r2 = R(i, ind(2))*R(i, ind(2));
            myfile << r0 << " " << r1 << " " << r2 << endl;

        }
        myfile.close();

    }




    return 0;
}

