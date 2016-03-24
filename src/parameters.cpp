#include "parameters.h"
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>

using namespace std;

//------------------------------------------------------------------------------------
//   Global constants
//------------------------------------------------------------------------------------
int OFTS_ORDER;
int OFS_ORDER;
int OTS_ORDER;
int MODEL_TYPE;
int REDUCED_NV;

//------------------------------------------------------------------------------------
//Misc
//------------------------------------------------------------------------------------
const std::string SSEPR = "---------------------------------";
const std::string SEPR = "---------------------------------------------------";

//------------------------------------------------------------------------------------
//   Precisions
//------------------------------------------------------------------------------------
/**
 *  \brief Sets an average precision in cout.
 **/
void coutmp()
{
    cout <<  setw(5) << setprecision(5) << std::showpos  <<  setiosflags(ios::scientific);
}

/**
 *  \brief Sets a big precision in cout.
 **/
void coutlp()
{
    cout <<  setw(5) << setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);
}

/**
 *  \brief Sets a small precision in cout.
 **/
void coutsp()
{
    cout <<  setw(5) << setprecision(3) << resetiosflags(ios::scientific);
}

//------------------------------------------------------------------------------------
//   Print
//------------------------------------------------------------------------------------
/**
 *  \brief Prints an array of double using cout.
 **/
void vector_printf(double *y, int n)
{
    coutmp();
    for(int i = 0; i < n; i ++) cout << i << "  " << y[i] << endl;
}

/**
 *  \brief Prints an array of complex double using cout.
 **/
void vector_complex_printf(cdouble *y, int n)
{
    coutmp();
    for(int i = 0; i < n; i ++) cout << i << "  " << creal(y[i]) << "  " << cimag(y[i]) << endl;
}

//------------------------------------------------------------------------------------
//   Norm
//------------------------------------------------------------------------------------
/**
 *  Euclidian norm computed on the first k components of a complex vector: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(cdouble z0[], int k)
{
    double res = cabs(z0[0])*cabs(z0[0]);
    for(int p = 1; p < k; p++) res += cabs(z0[p])*cabs(z0[p]);
    return(sqrt(res));
}

/**
 *  Euclidian norm computed on the first k components of a double vector: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(double z0[], int k)
{
    double res = cabs(z0[0])*cabs(z0[0]);
    for(int p = 1; p < k; p++) res += cabs(z0[p])*cabs(z0[p]);
    return(sqrt(res));
}

