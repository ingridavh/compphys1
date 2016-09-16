#ifndef GAUSSIAN_H
#define GAUSSIAN_H
double F (double x);
double v_exact (double x);
double decomp (double a, double b, double c_old, double beta_old);
double forward (double a, double beta_old, double f, double f_tilde_old);
double backward(double f_tilde_old, double c_old, double u_new, double beta_old);
double gaussian(int N=10, bool finderr = false, bool CPU = false, bool fileprint = false, bool CPU_avg = false);
#endif // GAUSSIAN_H

