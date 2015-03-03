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
void tic(void);
double toc(void);

struct timespec TIC_TIME;

#define MP_DOUBLE_TYPE 1
#define MP_COMPLEX_TYPE 2

// Choice of BASE type (double or complex)
//--------------------------------------------
/*
#define MP_BASE_TYPE MP_COMPLEX_TYPE
typedef complex double varType;
//*/

///*
#define MP_BASE_TYPE MP_DOUBLE_TYPE
typedef double varType;
//*/


// Choice of coefficient type for Fourier-Taylor series
//--------------------------------------------
typedef Ofs<varType> oftsVarType;


int main()
{

    int nr = 30;
    int nv = 4;

    FTDA::init(nv, nr);

    //----------------------------------------------------------------------------------------------------------
    // Test of the basic routines in Ofts
    //----------------------------------------------------------------------------------------------------------
    oneVariableTest(4, 4);
    oneVariableTest_Ofs(4, 4);

    return 0;
}


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
