#ifndef GSLC_H_INCLUDED
#define GSLC_H_INCLUDED

/**
 * \file gslc.h
 * \brief Additional operations on GSL objects.
 * \author BLB
 * \date May 2015
 * \version 1.0
 */

//std
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <complex.h>
#include <sstream>
#include <math.h>

//GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_eigen.h>

//Custom
#include "Constants.h"
#include "Config.h"


//----------------------------------------------------------------------------------------
// Matrix and vectors
//----------------------------------------------------------------------------------------
/**
 * \brief Transform a vector into a matrix with a given shift in the initial vector
 **/
void gslc_vectorToMatrix(gsl_matrix *m, const double y[], int rows, int columns, int shift);

/**
 * \brief Transform a vector into a complex matrix with a given shift in the initial vector
 **/
void gslc_vectorToComplexMatrix(gsl_matrix_complex *m, const double y[], int rows, int columns, int shift);

/**
 * \brief Transform a matrix into a vector with a given shift in the final vector
 **/
void gslc_matrixToVector(double y[], const gsl_matrix *m, int rows, int columns, int shift);

/**
 * \brief Transform a complex matrix into a vector with a given shift in the final vector
 **/
void gslc_complexMatrixToVector(double y[], const gsl_matrix_complex  *m, int rows, int columns, int shift);



//----------------------------------------------------------------------------------------
// Complex numbers
//----------------------------------------------------------------------------------------
/**
 * \brief Create a gsl_complex from 2 doubles
 **/
gsl_complex gslc_complex(double real, double imag);

/**
 * \brief Create a gsl_complex from a complex number
 **/
gsl_complex gslc_complex(cdouble x);

/**
 * \brief Create a complex number from a gsl_complex
 **/
cdouble gslc_complex(gsl_complex c);


//----------------------------------------------------------------------------------------
// Printing a real matrix
//----------------------------------------------------------------------------------------
/**
 * \brief Print a real matrix
 **/
void gslc_matrix_printf(const gsl_matrix *M);

/**
 * \brief Print an approximation of a real matrix: the coefficients under 1e-5 are displayed as null.
 **/
void gslc_matrix_printf_approx(const gsl_matrix *M);

/**
 * \brief Print a real matrix into a txt file.
 **/
void gslc_matrix_fprintf(const gsl_matrix *M, char* fileName);

//----------------------------------------------------------------------------------------
// Printing a complex matrix
//----------------------------------------------------------------------------------------
/**
 * \brief Print a complex matrix.
 **/
void gslc_matrix_complex_printf(const gsl_matrix_complex *M);

/**
 * \brief Print the real part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_printf_real(const gsl_matrix_complex *M);

/**
 * \brief Print the imag part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_printf_imag(const gsl_matrix_complex *M);

/**
 * \brief Print an approximation of a complex matrix: the coefficients under 1e-5 are displayed as null.
 **/
void gslc_matrix_complex_printf_approx(const gsl_matrix_complex *M);

/**
 * \brief Print a complex matrix into a txt file.
 **/
void gslc_matrix_complex_fprintf(const gsl_matrix_complex *M, char* fileName);

/**
 * \brief Print a complex matrix into a txt file, separating the real and imaginary part, in order to be easier to read.
 **/
void gslc_matrix_complex_fprintf_pretty(const gsl_matrix_complex *M, char* fileName);

/**
 * \brief Print the real part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_fprintf_real(const gsl_matrix_complex *M, char* fileName);

//----------------------------------------------------------------------------------------
// Printing a complex vector
//----------------------------------------------------------------------------------------
/**
 * \brief Print a complex vector into a txt file.
 **/
void gslc_vector_complex_fprintf(const gsl_vector_complex *V, char* fileName);

/**
 * \brief Print a complex vector.
 **/
void gslc_vector_complex_printf(const gsl_vector_complex *V);

//----------------------------------------------------------------------------------------
// Printing a vector
//----------------------------------------------------------------------------------------
/**
 * \brief Print a real vector.
 **/
void gslc_vector_printf(const gsl_vector *V);

//----------------------------------------------------------------------------------------
// Printing an eigensystem
//----------------------------------------------------------------------------------------
/**
 * \brief Print an eigensystem with the eigenvalues given as a complex vector.
 **/
void gslc_eigensystem_printf(gsl_vector_complex *eval, gsl_matrix_complex *evec, int M);

/**
 * \brief Print an eigensystem with the eigenvalues given as a complex matrix (on the diagonal).
 **/
void gslc_eigensystem_printf(gsl_matrix_complex *eval, gsl_matrix_complex *evec, int M);

/**
 * \brief Print an eigensystem in a txt file with the eigenvalues given as a complex vector.
 **/
void gslc_eigensystem_fprintf(gsl_vector_complex *eval, gsl_matrix_complex *evec, int M, char* fileName);

/**
 * \brief Print an eigensystem in a txt file with the eigenvalues given as a complex matrix (on the diagonal).
 **/
void gslc_eigensystem_fprintf(gsl_matrix_complex *eval, gsl_matrix_complex *evec, int M, char* fileName);

//----------------------------------------------------------------------------------------
// Reading GSL objects
//----------------------------------------------------------------------------------------
/**
 * \brief Read a complex vector from a txt file, obtained with the routine gslc_matrix_complex_fprintf(const gsl_matrix_complex *M, char* fileName).
 **/
void glsc_matrix_complex_read(gsl_matrix_complex *M, std::string filename);

//----------------------------------------------------------------------------------------
// Views of a matrix
//----------------------------------------------------------------------------------------
/**
 * \brief Set the kth column of the complex matrix M in the complex vector V using the GSL "views" structures.
 **/
void gslc_matrix_complex_column(gsl_vector_complex *V, gsl_matrix_complex *M, int k);

/**
 * \brief Set the kth row of the complex matrix M in the complex vector V using the GSL "views" structures.
 **/
void gslc_matrix_complex_row(gsl_vector_complex *V, gsl_matrix_complex *M, int k);

/**
 * \brief Set the complex vector V in the kth column of the complex matrix M using the GSL "views" structures.
 **/
void gslc_matrix_complex_column_V(gsl_matrix_complex *M, gsl_vector_complex *V, int k);


//----------------------------------------------------------------------------------------
// Misc manipulations
//----------------------------------------------------------------------------------------
/**
 * \brief Gives the infinity norm of a complex matrix M
 **/
double gslc_matrix_complex_infinity_norm(gsl_matrix_complex *M);

/**
 *  \brief Normalization: xm = xm/norm(xm)
 **/
void gslc_vector_complex_normalize(gsl_vector_complex *xm);

/**
 * \brief Isolate the real part of a complex matrix M: Mr = real(M)
 **/
void gslc_matrix_complex_real(gsl_matrix *Mr, gsl_matrix_complex *M);

/**
 * \brief Copy the conjugate of a complex vector xm into xc.
 **/
void gslc_vector_complex_conjugate_memcpy(gsl_vector_complex *xc, gsl_vector_complex *xm);

//----------------------------------------------------------------------------------------
// Specific routines for Monodromy and STM matrices manipulations
//----------------------------------------------------------------------------------------
/**
 * \brief Delete one raw and one column of a given square GSL matrix
 **/
gsl_matrix_complex * gslc_matrix_complex_deleteRC(gsl_matrix_complex *M, int k);

/**
 * \brief Inverse transformation of a vector during Wielandt deflation: w = 1/vw*(w + vx/(vw-vx)*(z.w)*x)
 **/
void gslc_wielandt_inv_trans(gsl_vector_complex *w, gsl_complex vw, gsl_vector_complex const *x, gsl_complex vx, gsl_vector_complex *z);

/**
 *  \brief Inverse a symplectic complex matrix S0 into Sinv:
 *
 *  S0  =  S11 S12
 *         S21 S22
 *
 *  S-1 =  S22 -S12T
 *        -S21T S11T
 **/
void gslc_matrix_complex_symplectic_inverse(const gsl_matrix_complex *S0, gsl_matrix_complex *Sinv);

/**
 *  \brief Inverse a symplectic real matrix S0 into Sinv:
 *
 *  S0  =  S11 S12
 *         S21 S22
 *
 *  S-1 =  S22 -S12T
 *        -S21T S11T
 **/
void gslc_matrix_symplectic_inverse(const gsl_matrix *S0, gsl_matrix *Sinv);

/**
 *  \brief Transpose + conjugate the complex matrix S intro SH: SH = S^H != S^T
 **/
void gslc_matrix_complex_H_memcpy(gsl_matrix_complex *SH, const gsl_matrix_complex *S);

/**
 *  \brief Setting the 2n*2n complex matrix J equal to the fundamental symplectic matrix:
 *
 *  real(J) = | 0  In |   and imag(J) = 0
 *            |-In 0  |
 **/
void glsc_matrix_complex_set_J(gsl_matrix_complex *J);

/**
 * \brief Setting the 2n*2n complex matrix J equal to the fundamental symplectic matrix:
 *
 *  J = | 0  In |
 *      |-In 0  |
 **/
void glsc_matrix_set_J(gsl_matrix *J);

/**
 * \brief Use the symmetry of the QBCP to compute the stable (resp. unstable) from the unstable (resp. stabl) eigenvector
 **/
void gslc_vector_complex_symcpy(gsl_vector_complex *VEP2, gsl_vector_complex *VEP1);

/**
 * \brief Matrix-vector product when the matrix is given as a product of matrices: ym = DAT[1]...DAT[M] * xm, with
 *        ym  a 6x1 vector
 *        xm  a 6x1 vector
 *        DAT a 6x6xM matrix
 *        Careful: the product of matrices DAT is indexed from 1 to M, and not 0 to M-1.
 **/
void gslc_matrix_vector_product(gsl_matrix_complex **DAT, const gsl_vector_complex *xm, gsl_vector_complex *ym, int M);

/**
 * \brief Matrix-matrix product when the matrix is given as a product of matrices: ym = DAT[1]...DAT[M] * xm, with
 *        ym  a 6x6 matrix
 *        xm  a 6x6 matrix
 *        DAT a 6x6xM matrix
 *        Careful: the product of matrices DAT is indexed from 1 to M, and not 0 to M-1.
 **/
void gslc_matrix_matrix_product(gsl_matrix_complex **DAT, const gsl_matrix_complex *xm, gsl_matrix_complex *ym, int M);

/**
 * \brief Matrix inverse-vector product when the matrix is given as a product of matrices: ym = DAT[M]^(-1)...DAT[1]^(-1) * xm, with
 *        ym  a 6x1 vector
 *        xm  a 6x1 vector
 *        DAT a 6x6xM matrix
 *        Careful: the product of matrices DAT is indexed from 1 to M, and not 0 to M-1.
 **/
void gslc_matrix_vector_invproduct(gsl_matrix_complex **DAT,  const gsl_vector_complex *xm, gsl_vector_complex *ym, int M);

/**
 * \brief Initialize and return a matrix as a product of complex matrices (DAT)
 *        CAREFUL: the array of matrices DAT is shifted of one & so that the storage is easier in other routines (e.g. vepro)
 *        ==> DAT has to be used from DAT[1] to DAT[M] and DAT[0] is USELESS.
 **/
gsl_matrix_complex ** gslc_matrix_complex_product_alloc(int size1, int size2, int M);

/**
 * \brief Free a matrix as a product of complex matrices (DAT)
 *        CAREFUL: the array of matrices DAT is shifted of one & so that the storage is easier in other routines (e.g. vepro)
 *        ==> DAT has to be used from DAT[1] to DAT[M] and DAT[0] is USELESS.
 **/
void gslc_matrix_complex_product_free(gsl_matrix_complex ** P, int M);

/**
 * \brief Initialize and return a matrix as a product of real (DAT)
 *        NO shift in this routine, contrary to gslc_matrix_complex_product_alloc
 **/
gsl_matrix ** gslc_matrix_array_alloc(int size1, int size2, int M);

/**
 * \brief Initialize and return an array of GSL matrices
 **/
gsl_matrix** gslc_matrix_array_calloc(int size1, int size2, int M);

/**
 * \brief Free a matrix as a product of real (DAT)
 *        NO shift in this routine, contrary to gslc_matrix_complex_product_alloc
 **/
void gslc_matrix_array_free(gsl_matrix ** P, int M);


//----------------------------------------------------------------------------------------
// Utilitary routines for C vector are also grouped here
//----------------------------------------------------------------------------------------
/**
 *  \brief Prints an array of double using cout.
 **/
void vector_printf(double *y, int n);

/**
 *  \brief Prints an array of complex double using cout.
 **/
void vector_complex_printf(cdouble *y, int n);

/**
 *  Euclidian norm computed on the first k components of a complex vector: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(cdouble z0[], int k);

/**
 *  Euclidian norm computed on the first k components of a double vector: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(double z0[], int k);




#endif // GSLC_H_INCLUDED
