#ifndef JACOBI_ROTATION
#define JACOBI_ROTATION
#include <armadillo>

arma::mat Jacobi (int N, arma::mat A);
int eps_test(int N, arma::mat B, double eps=1e-8);
arma::vec find_max(arma::mat A,int N);
arma::vec find_trig(arma::mat A, int k, int l);
arma::mat B_new(arma::mat A, int N, int k, int l, double c, double s);
//arma::mat jacobi_tol (int N, arma::mat B, double eps=1e-8);

#endif // JACOBI_ROTATION

