#include <iostream>
#include "quantummc.h"
#include <armadillo>
#include <string>

using namespace std;

int main()
{
    //Open the output file
    string filename = "tiffanys_breakfast.txt";
    ofstream ofile;
    if(!ofile.is_open()){
        ofile.open(filename.c_str(), ofstream::out);
        if(!ofile.good()){
            cout << "Error opening file " << filename << " Aborting!" << endl;
            terminate();

        }
    }

    outfile << "Alpha " << ", omega " << ", Ek " << ", Ep" << endl;

    //Set alpha values
    int T1 = 1; int T2 = 2;
    double alpha_min = 0.5;
    double alpha_max = 1.5;
    double alpha_step = 0.05;
    double beta_step = 0.05;
    int num_steps =  (alpha_max-alpha_min)/alpha_step;

    //Initialize instance
    QuantumMC dots(2);

    //Add energy perturbation
    dots.add_perturbation();

    //Generate test values for alpha and beta
    double a_list [num_steps]; double b_list [2*num_steps];

    for(int i = 0; i < num_steps; i++){
        a_list[i] = alpha_min + i*alpha_step;
        b_list[i] = alpha_min + i*beta_step;
    }

    //Run test to find optimal values for alpha and beta
    //dots.get_low(a_list, b_list, num_steps);

    double omegas [3] = {0.01, 0.5, 1.0};

//    for (int i = 0; i < 3; i++){
//        dots.m_omega = omegas[i];
//        cout << "omega " << omegas[i] << endl;
//        dots.get_low(a_list, b_list, num_steps);
//    }

    //Omegas to find kinetic energy as a function of omega
    double omegas_long [10] = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.55, 0.65, 0.8, 1};
    for(int i = 0; i < 10; i++){
        dots.m_T = T1;
        dots.m_omega = omegas_long[i];
        double a_min = dots.minimize_alpha(a_list, num_steps);
        dots.m_alpha = a_min;
        dots.runMCintegration(T1);
        outfile << a_min << " , " <<
    }




//    ofile << "alpha " << " energy " << " variance" <<endl;

    for(double a=alpha_min; a < alpha_max; a += alpha_step){
        dots.m_alpha = a;
        dots.runMCintegration(T1);
        double energy = dots.m_energy;
        double variance = dots.m_var;
        ofile << a << " "<<energy << " " << variance << endl;
    }


    ofile.close();

    return 0;
}
