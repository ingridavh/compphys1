
//Function for f-values
double F (double x)
{
    double y = -10*x;
    double F_value = 100* exp (y);
    return F_value;
}

//Function for exact v-values
double v_exact (double x)
{
    double v = 1 - (1 - exp (-10))*x - exp(-10*x);
    return v;
}


//Step 1: decomposition function
double decomp (double a, double b, double c_old, double beta_old)
{
    double beta = b - (a*c_old)/beta_old;
    return beta;
}

//Step 1: Forward substitution function
double forward (double a, double beta_old, double f, double f_tilde_old)
{
    double f_tilde = f - (a*f_tilde_old)/beta_old;
    return f_tilde;
}

//Step 2: backward substitution function
double backward(double f_tilde_old, double c_old, double u_new, double beta_old)
{
    double u = (f_tilde_old - c_old*u_new)/beta_old;
    return u;
}



void gaussian()
{
    //Number of step lengths and h
    int N = 1000;
    double h = 1./(N+1);
    double x_0 = 0;

    //Declare matrix A, A_new and arrays v, f and f_tilde
    double *a = new double[N+2];
    double *b = new double[N+2];
    double *c = new double[N+2];
    double *beta = new double[N+2];
    double *v = new double[N+2];
    double *f = new double[N+2];
    double *f_tilde = new double[N+2];
    double *v_ex = new double[N+2];

    //Fill inn values of A
    for(int i = 0; i < (N+2); i++)
    {
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
        double x = x_0 + i*h;
        f[i] = h*h*F(x);
        v_ex[i] = v_exact(x);
    }

    //Set boundary conditions
    beta[1] = b[1];
    f_tilde[1] = f[1];
    v[0] = 0;
    v[N+1] = 0;

    //Step 1
    for (int i = 2; i < (N+1); i++)
    {
        beta[i] = decomp(a[i], b[i], c[i-1], beta[i-1]);
        f_tilde[i] = forward(a[i], beta[i-1],f[i], f_tilde[i-1]);
    }

    v[N] = f_tilde[N]/beta[N];

    //Step 2
    for(int i = N; i > 1; i--)
    {
        v[i-1] = backward(f_tilde[i-1], c[i-1], v[i], beta[i-1]);
    }


    //Write results to a file p1_result_N (n-value).txt
    ofstream myfile;
    myfile.open ("p1_result_N1000.txt", ofstream::out);
    if(!myfile.good()){
        cout << "Dette gikk galt" << endl;
        return 1;
    }

    //For loop that writes results to file
    for (int i = 0; i < (N+2); i++)
    {
        myfile << v[i] << " " << v_ex[i] << "\n";
    }
    myfile.close();
    return 0;
}


