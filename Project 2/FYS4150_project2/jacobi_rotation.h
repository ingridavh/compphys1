#ifndef JACOBI_ROTATION
#define JACOBI_ROTATION
#include <armadillo>

int eps_test(arma::mat A, double eps=1e-8);
arma::vec find_max(arma::mat A);
arma::vec find_trig(arma::mat A, int k, int l);
void rotate(arma::mat &A, arma::mat &R, int k, int l);

#endif // JACOBI_ROTATION

