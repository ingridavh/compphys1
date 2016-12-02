#include <iostream>
#include "quantummc.h"

using namespace std;

int main()
{
    cout << "Hello World!" << endl;

    QuantumMC dots(2);
    dots.add_perturbation();
    dots.runMCintegration();

    return 0;
}

