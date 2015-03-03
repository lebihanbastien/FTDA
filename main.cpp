#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <complex.h>

#include "ofs.h"
#include "otsh.h"
#include "ots.h"
#include "oftsh.h"
#include "ofts.h"
#include "ftda.h"

#include "qbtbp.h"
#include "qbtbp_ofs.h"
#include "tests.h"

//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl_interp.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

//Custom
#include "custom_ode.h"


using namespace std;

typedef Ofs< complex > Ofsc;
typedef Ofs<double>  Ofsd;

void time_mul(int deg_max, int nvar);
void time_mul_FT(int deg_max, int nvar);

#ifdef __APPLE__
time_t TIC_TIME;
#else
struct timespec TIC_TIME;
#endif

void tic(void)
{
#ifdef __APPLE__
    TIC_TIME = clock();
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &TIC_TIME);
#endif
}

double toc(void)
{
#ifdef __APPLE__
    time_t toc_time;
    toc_time = clock();
    return (double)(toc_time - TIC_TIME) / CLOCKS_PER_SEC;
#else
    struct timespec toc_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &toc_time);
    return (double)(toc_time.tv_sec - TIC_TIME.tv_sec) +
           (double)(toc_time.tv_nsec - TIC_TIME.tv_nsec) / 1e9;
#endif
}

#define MP_DOUBLE_TYPE 1
#define MP_COMPLEX_TYPE 2

/*
#define MP_BASE_TYPE MP_COMPLEX_TYPE
typedef complex double varType;
//*/

///*
#define MP_BASE_TYPE MP_DOUBLE_TYPE
typedef double varType;
//*/

typedef Ots<varType> oftsVarType;

#if MP_BASE_TYPE == MP_DOUBLE_TYPE
#define MY_CONST (1.0)
#endif

#if MP_BASE_TYPE == MP_COMPLEX_TYPE
#define MY_CONST (1.0 + _Complex_I*1.0)
#endif





void qbtbp (Ofts<double complex> &ofts_z, Ofts<double complex> &ofts_Z, double n, double ms, double as);


int main()
{

    int nr = 30;
    int nv = 4;

    FTDA::init(nv, nr);

    //----------------------------------------------------------------------------------------------------------
    //CPU Time computations
    //----------------------------------------------------------------------------------------------------------
    /*
    int nvar, deg_max;

    for (nvar = 2; nvar <= 4; nvar++)
    {
        for (deg_max = 2; deg_max <= 16; deg_max++)
        {
            time_mul(deg_max, nvar);
        }
        printf("\n");
    }
    //*/


    //----------------------------------------------------------------------------------------------------------
    //CPU Time computations for Fourier-Taylor series
    //----------------------------------------------------------------------------------------------------------
    //*
    int nvar, deg_max;
    for (nvar = 2; nvar <= 4; nvar++)
    {
        for (deg_max = 2; deg_max <= 6; deg_max++)
        {
            time_mul_FT(deg_max, nvar);
        }
        printf("\n");
    }
    //*/



    //----------------------------------------------------------------------------------------------------------
    // Test of the basic routines in Ofts
    //----------------------------------------------------------------------------------------------------------
    //oneVariableTest(4, 4);
    //oneVariableTest_Ofs(4, 4);


    //----------------------------------------------------------------------------------------------------------
    // Quasi Bicircular TBP resolution (bj's and cj's) - Taylor Implementation
    //----------------------------------------------------------------------------------------------------------
    /*
        Ofts< Ots<varType> > ofts_z(1,30,2,40);
        Ofts< Ots<varType> > ofts_Z(1,30,2,40);

        double as = 388.81114;
        double omegas = 0.925195985520347;
        double ns =1-omegas;
        double ms = ns*ns*pow(as,3.0)-1;
        double n = 1-ns;

        tic();
        qbtbp(ofts_z, ofts_Z, n, ms, as, ns);
        double diff = toc();

        cout << "CPU time: " << diff << endl;
    //*/

    //----------------------------------------------------------------------------------------------------------
    // Quasi Bicircular TBP resolution (bj's and cj's) - Fourier Implementation
    //----------------------------------------------------------------------------------------------------------
    /*
    int nf = 20;
    Ofts< Ofs<varType> > ofts_z(1,30,2,nf);
    Ofts< Ofs<varType> > ofts_Z(1,30,2,nf);

    Ofs<double complex> bjc(nf);
    Ofs<double complex> cjc(nf);

    double as = 388.81114;
    double omegas = 0.925195985520347;
    double ns =1-omegas;
    double ms = ns*ns*pow(as,3.0)-1;
    double n = 1-ns;
    double mu = 0.0121505816;
    tic();
    qbtbp_ofs(bjc, cjc, ofts_z, ofts_Z, n, ms, as, ns);
    double diff = toc();

    cout << "CPU time: " << diff << endl;
    //*/

    //----------------------------------------------------------------------------------------------------------
    // Reading into the students results file
    //----------------------------------------------------------------------------------------------------------
    /*
    double as = 388.81114;
    double omegas = 0.925195985520347;
    double ns =1-omegas;
    double ms = ns*ns*pow(as,3.0)-1;
    double n = 1-ns;
    double mu = 0.0121505816;


    ifstream readStream;
    ofstream writeStream;
    int nf = 30;
    Ofs<double complex> bj_student(nf);
    Ofs<double complex> cj_student(nf);

    //Reading the bjs
    readStream.open("data/bj_student.txt");
    double currentNumber;
    //Order 0
    readStream >> currentNumber; //currentNumber = 0;
    readStream >> currentNumber;
    bj_student.setCoef(currentNumber, 0);
    for(int i=1; i<= nf; i++)
    {
        readStream >> currentNumber; //currentNumber = i;
        readStream >> currentNumber;
        bj_student.setCoef(currentNumber, i);
        readStream >> currentNumber;
        bj_student.setCoef(currentNumber, -i);
    }
    readStream.close();

    //Writing the bjs
    writeStream.open("data/bj_student_copy.txt");
    writeStream << bj_student << endl;
    writeStream.close();


    //Reading the cjs
    readStream.open("data/cj_student.txt");
    //Order 0
    readStream >> currentNumber; //currentNumber = 0;
    readStream >> currentNumber;
    cj_student.setCoef(currentNumber, 0);
    for(int i=1; i<= nf; i++)
    {
        readStream >> currentNumber; //currentNumber = i;
        readStream >> currentNumber;
        cj_student.setCoef(currentNumber, i);
        readStream >> currentNumber;
        cj_student.setCoef(currentNumber, -i);
    }
    readStream.close();

    //Writing the cjs
    writeStream.open("data/cj_student_copy.txt");
    writeStream << cj_student << endl;
    writeStream.close();

    */

    //----------------------------------------------------------------------------------------------------------
    // Integration & Analytical Vs Numerical results
    //----------------------------------------------------------------------------------------------------------
    /*
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    custom_ode_structure ode_s;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent; //Brent-Dekker root finding method
    //Ode solver parameters
    double param[2];
    param[0] = mu;
    param[1] = ms;
    //General structures
    init_ode_structure(&ode_s, T, T_root, 1e-14, 1e-14, 1e-12, 1e-12, 8, 1e-6,  qbtbp_derivatives, NULL, param);
    //Analytical Vs Numerical results
    analyticVsnumeric(2*M_PI, bjc, cjc, ode_s, as, n, ns);
    //analyticVsnumeric(2*M_PI, bj_student, cj_student, ode_s, as, n, ns);  //use it to test the results from Josep's student implementation
    //*/


    return 0;
}


void time_mul(int deg_max, int nvar)
{
    int i, j, N;
    double diff;

    Ots<varType> p (nvar, deg_max);
    Ots<varType> a (nvar, deg_max);
    Ots<varType> b (nvar, deg_max);

    p.setCoefs(MY_CONST);
    a.setCoefs(MY_CONST);
    b.setCoefs(MY_CONST);


    // Call once to get rough estimate of time.
    tic();
    for (i = 0; i <= deg_max; i++)
    {
        p.sprod(a,b,i);
    }
    diff = toc();

    // Select number of times to call.
    if (diff > 10.0)
    {
        N = 2;
    }
    else if (diff > 1.0)
    {
        N = 10;
    }
    else if (diff > 0.1)
    {
        N = 100;
    }
    else
    {
        N = 1000;
    }

    // Perform timing test.
    tic();
    for (j = 0; j < N; j++)
    {
        for (i = 0; i <= deg_max; i++)
        {
            p.sprod(a,b,i);
        }
    }
    diff = toc();

    diff = diff / N;
    printf("nvar = %2d, deg = %2d (1000 calls): %0.6f seconds\n",nvar, deg_max, 1000*diff);
}

void time_mul_FT(int deg_max, int nvar)
{
    int i, j, N;
    double diff;

    Ofts< oftsVarType > p(nvar, deg_max, 2, deg_max);
    Ofts< oftsVarType > a(nvar, deg_max, 2, deg_max);
    Ofts< oftsVarType > b(nvar, deg_max, 2, deg_max);

    p.setCoefs(MY_CONST);
    a.setCoefs(MY_CONST);
    b.setCoefs(MY_CONST);

    // Call once to get rough estimate of time.
    tic();
    for (i = 0; i <= deg_max; i++)
    {
        p.sprod(a,b,i);
    }
    diff = toc();

    // Select number of times to call.
    if (diff > 10.0)
    {
        N = 2;
    }
    else if (diff > 1.0)
    {
        N = 10;
    }
    else if (diff > 0.1)
    {
        N = 100;
    }
    else
    {
        N = 1000;
    }

    // Perform timing test.
    tic();
    for (j = 0; j < N; j++)
    {
        for (i = 0; i <= deg_max; i++)
        {
            p.sprod(a,b,i);
        }
    }
    diff = toc();

    diff = diff / N;
    printf("nvar = %2d, deg = %2d (1000 calls): %0.6f seconds\n",nvar, deg_max, 1000*diff);
}
