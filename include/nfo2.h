#ifndef NFO2_H_INCLUDED
#define NFO2_H_INCLUDED

/**
 * \file nfo2.h
 * \brief Computes the complete change of coordinates for the Normalized-Centered Hamiltonian of the QBCP.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */

//GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>

//Custom
#include "qbcp.h"
#include "define_env.h"
#include "parameters.h"

#include <stdio.h>
#include <iostream>
using namespace std;

#define INVERSE_SYMP 1
#define INVERSE_GSL 2


//-------------------------------------------------------------------------------------------------------
// Main routine
//-------------------------------------------------------------------------------------------------------
/**
 *  \brief Main routine. Compute the complete change of coordinates.
 */
void nfo2(QBCP_L &qbcp_l, int isStored);

/**
 *  \brief Main routine. Compute the complete change of coordinates.
 */
void nfo2_QBP(QBCP_L &qbcp_l, int isStored);

/**
 *  \brief Continuation routine
 */
void continuation(QBCP_L &qbcp_l, int isStored);

//-------------------------------------------------------------------------------------------------------
//Compute the STM/Monodromy Matrix
//-------------------------------------------------------------------------------------------------------
/**
 *  \brief Get the STM matrix at time \c t1, using GSL tools.
 */
void stmMatrix(const double y[], gsl_odeiv2_driver *d, double t1, gsl_matrix* M);
/**
 *  \brief Get the STM matrix at time \c t1 in complex format, using GSL tools.
 */
void stmComplexMatrix(const double y[], gsl_odeiv2_driver *d, double t1, gsl_matrix_complex* M, double ys[]);
/**
 *  \brief Get the STM matrix on a \c M steps grid in complex format, using GSL tools.
 */
void stmSteppedComplexMatrix(const double y[], gsl_odeiv2_driver *d, double t1, int M, gsl_matrix_complex **DAT);
/**
 *  \brief Get the STM matrix on a \c M steps grid in complex format, using GSL tools. Alternative version that uses the symmetries of the native orbit.
 */
void stmSteppedComplexMatrix_alt(const double y[], gsl_odeiv2_driver *d, double t1, int M, gsl_matrix_complex **DAT);
/**
 *  \brief Computing the Monodromy matrix from the stepped STM stored in DAT[1, ..., M].
 *  \param MMc: the matrix to update.
 *  \param DAT: a double pointer to the gsl_matrix_complex* array to update
 *  \param M: the number of steps
 **/
void monoProd(gsl_matrix_complex* MMc, gsl_matrix_complex** DAT, int M);
/**
 *  \brief Direct diagonalization of the monodromy matrix
 *  \param MM    the monodromy matrix in real format
 *  \param evecd output: the eigenvectors in matrix format
 *  \param evald output: the eigenvalues  in vector format
 *  \param evmd  output: the eigenvalues  in matrix format (along the diagonal)
 **/
void monoDiag(gsl_matrix* MM, gsl_matrix_complex *evecd, gsl_vector_complex *evald, gsl_matrix_complex *evmd);

//-------------------------------------------------------------------------------------------------------
// Matrix decomposition
//-------------------------------------------------------------------------------------------------------
/**
 *  \brief Obtain the monodromy matrix decomposition from a monodromy matrix given as a product of matrices.
 */
void monoDecomp(gsl_odeiv2_driver *d,           //driver for odeRK78
                const double y0[],              //initial conditions
                QBCP_L *qbcp_l,                //Four-Body problem
                int M,                          //Number of steps for Monodromy matrix splitting in a product of matrices
                int STABLE_DIR_TYPE,            //Type of computation of the stable direction
                gsl_matrix_complex* MMc,        //Final monodromy matrix in complex form                                    |
                gsl_matrix* MM,                 //Final monodromy matrix in real form                                       |
                gsl_matrix_complex *Dm,         //The eigenvalues  are stored in Dm(i,i), i = 0,...,5                       |   Outputs
                gsl_matrix_complex *DB,         //DB = 1/T*log(Dm)                                                          |
                gsl_matrix_complex *B,          //B = S*DB*Sinv                                                             |
                gsl_matrix_complex *S,          //The eigenvectors are stored on  S(i,*), i = 0,...,5                       |
                gsl_matrix *JB,                 //Real Jordan form of B                                                     |
                gsl_matrix *Br,                 //Real form of B computed from the Jordan Form                              |
                gsl_matrix *R,                  //B = R*JB*Rinv                                                             |
                gsl_matrix_complex** DAT,       //Contains the STM at each of the M steps                                   |
                int isStored);

/**
 *  \brief Perform a permutation within S and Dm.
 */
void permutationS(gsl_matrix_complex *S,
                 gsl_matrix_complex *Dm,
                 gsl_matrix_complex *evecr,
                 gsl_vector_complex *evalr,
                 gsl_vector_complex *eigenVu,
                 gsl_vector_complex *eigenVs,
                 gsl_complex eigenLu,
                 gsl_complex eigenLs,
                 int *keymap);


/**
 *  \brief Obtain the real Jordan form of the matrix B.
 */
void monoDecompLog(gsl_matrix_complex *S,  gsl_matrix_complex *Dm, double T, double n, gsl_matrix_complex *DB, gsl_matrix_complex *B, gsl_matrix *Br, gsl_matrix *JB, gsl_matrix *R);

/**
 *  \brief Diagonalization of the central part of VECS(*,*,0), wich at this point, has no unstable component under M = DAT[M]*DAT[M-1]*...*DAT[1].
 *
 *    This central part is denoted (v1 v2 v3 v4), with (v1 v2) associated to the xy center, (v3 v4) to the z center.
 *    Providing that (vi, vj) actually spans the 2th dimension subspace (either xy or z center), we have that:
 *
 *        vi = k1 u1 + k2 u1b
 *        vj = k3 u1 + k4 u1b
 *
 *        where (u1, u1b) is the couple of eigenvectors associated with either the xy or z center, with eigenvalues (l1, l1b).
 *
 *        Thus, the routines inverses following code inverses the following system:
 *
 *
 *             |  vi  |   |   I    I     0    0   |  |k1 u1  |       |k1 u1  |
 *             |  vj  |   |   0    0     I    I   |  |k2 u1b |       |k2 u1b |
 *             | M*vi | = | l1*I l1b*I   0    0   |  |k3 u1  | = M64 |k3 u1  |
 *             | M*vj |   |   0    0   l1*I l1b*I |  |k4 u1b |`      |k4 u1b |
 *
 *        with M the monodromy matrix, in order to get u1 and u1b.
 *
 **/
void diagcp(gsl_matrix_complex **DAT, gsl_matrix_complex *VECS_0, gsl_vector_complex *DECSl, gsl_matrix_complex *DECSv, int M);

//-------------------------------------------------------------------------------------------------------
//Resolution of eigensystem: power methods on single matrices
//-------------------------------------------------------------------------------------------------------
/**
 *  \brief Power Method on a single matrix.
 */
void powerMethod(gsl_matrix_complex *Minit, double prec, int sizeM, int version, gsl_vector_complex *eigenV, gsl_complex *eigenL);

/**
 *  \brief Inverse Power Method on a single matrix.
 */
void inversePowerMethod(gsl_matrix_complex *Minit, double prec, int sizeM, int version, gsl_vector_complex *eigenV, gsl_complex *eigenL);

/**
 *  \brief Shifted Inverse Power Method on a single matrix.
 */
void shiftedPowerMethod(gsl_matrix_complex *Minit, double prec, int sizeM, int version, gsl_complex shift, gsl_vector_complex *eigenV, gsl_complex *eigenL);

/**
 *  \brief Shifted Inverse Power Method on a single matrix with alternative update of the eigenvalue
 */
void shiftedPowerMethodEigenvalueUpdate(gsl_matrix_complex *Minit, double prec, int sizeM, gsl_complex shift, gsl_vector_complex *eigenV, gsl_complex *eigenL);

/**
 *  \brief Power Method on a product of matrices. note that VAPL = log10(VAP)
 */
void dipowers(gsl_matrix_complex **DAT, int M, int N, gsl_vector_complex *VEP, gsl_complex *VAP, gsl_complex *VAPL,  int IS);

/**
 *  \brief Decomposition of a product of matrices. Bad precision
 *
 *  Inspired by Josep's Fortran code. WARNING: does NOT work for the central part, monoDecomp is used instead.
 *  Recall that VECS(*,j,k)= DAT(*,*,k)*VECS(*,j,k-1) except for a multiplicative constant for each direction.
 */
void vepro(gsl_matrix_complex **DAT, int M, gsl_matrix_complex **VECS, gsl_matrix_complex *DECS, gsl_matrix_complex *DECSv, gsl_vector_complex *DECSl);

/**
 *  \brief Decomposition of the Monodromy matrix through the use of GSL routines from a monodromy matrix given as a single matrix.
 */
void gslMonoDecomp(gsl_odeiv2_driver *d, const double y0[], double tend, gsl_matrix_complex* gsl_MMc, gsl_matrix* gsl_MM, gsl_vector_complex *eval, gsl_matrix_complex *evec);

/**
 * \brief From J. Masdemont Fortran code. Change the intput VECS in such a way that
 *        the new vectors VECS(*,j,m) j=2,5 span the same 4th dimensional space that
 *        the vectors VECS(*,j,0) j=2,5.
 *
 *  1. The routine updates the matrix DECS so that VECS(*,j,0)*DECS = VECS(*,j,M), for j = 2,5.
 *
 *  2. If isSameBase == true: Given m+1 sets of basis of 6 vectors stored by columns in VECS(6,6,0:m)
 *  assuming that VECS(*,j,k+1)=a_k*VECS(*,j,k) where a_k is a certain
 *  6*6 matrix, and VECS(*,j,m)=di*VECS(*,j,0) for j=0,1, where
 *       - d1=delta(0,0)*delta(0,1)*..*delta(0,m-1) and
 *       - d2=1.d0/(delta(1,0)*delta(1,1)*...*delta(1,m-1))
 *  i.e. VECS(*,0,0) and VECS(*,1,0) are eigenvectors of the product
 *  of the a_k matrices, and delta must be given at the input,
 *  this routine modifies the vectors VECS(*,j,k) j=2,5 k=0,m, keeping
 *  always VECS(*,j,k+1)=a_k*VECS(*,j,k), in such a way that the new
 *  vectors VECS(*,j,m) j=2,5 span the same 4th dimensional space that
 *  the vectors VECS(*,j,0) j=2,5.
 *
 **/
void diahip(gsl_matrix_complex **VECS, double **DELTA, gsl_matrix_complex *DECS, int M, int isSameBase);

/**
 *  \brief Auxiliary routine. It orthonomalizes the bases VECS(*,*,K), for all K=0,M.
 **/
void ortho(gsl_matrix_complex **VECS, int IK, int M);
//-------------------------------------------------------------------------------------------------------
// Matrix integration
//-------------------------------------------------------------------------------------------------------
/**
 *  \brief Obtain the FFT decomposition of various matrices P, Q = inv(P) FT11, FT12, FT21 and FT22.
 */
void nfo2coc(gsl_odeiv2_driver *d, const double y0[], QBCP_L *qbcp_l, gsl_matrix* Br, gsl_matrix* R, gsl_matrix* JB, int N, int isStored);

/**
 *  \brief Integrate the various matrices referenced on a N point grid to seed FFT process.
 */
void nfo2Int(gsl_odeiv2_driver *d,
             const double y0[],
             QBCP_L *qbcp_l,
             gsl_matrix* Br,
             gsl_matrix* R,
             gsl_matrix* JB,
             int N,
             gsl_matrix** P,                //6*6   |
             gsl_matrix** Pfb,              //6*6   |
             gsl_matrix** Q,                //6*6   |
             gsl_matrix** Qfb,              //6*6   |
             gsl_matrix** FT11,             //3*3   |
             gsl_matrix** FT12,             //3*3   | Outputs
             gsl_matrix** FT21,             //3*3   |
             gsl_matrix** FT22,             //3*3   |
             gsl_matrix** G1  ,             //2*2   |
             gsl_matrix** Xe  ,             //2*1   |
             gsl_matrix** Xm  ,             //2*1   |
             gsl_matrix** Xs  );            //2*1   |

/**
 *  \brief FFT of the coefficients of a given matrix P obtained on a N points grid
 */
void nfo2FFT(Ofsc& xFFT, QBCP_L *qbcp_l, gsl_matrix** P, gsl_matrix** Pt, /*gsl_matrix* R,*/ int N, int fftN, int fftPlot, string filename, int flag, int isStored);

/**
 *  \brief Test of the validity of the FFT on a fftPlot points grid.
 */
void nfo2Test(Ofsc& xFFT, QBCP_L *qbcp_l, gsl_matrix** P, /*gsl_matrix* R,*/ int fftPlot, int i0, int j0);


//-------------------------------------------------------------------------------------------------------
// Testing
//-------------------------------------------------------------------------------------------------------
/**
 *  \brief Periodicity test of P.
 */
void nfo2PerTest(gsl_odeiv2_driver *d, const double y0[], QBCP_L *qbcp_l, gsl_matrix* R);

/**
 *  \brief Symmetry test on P.
 */
void nfo2SymTest(gsl_odeiv2_driver *d, const double y0[], QBCP_L *qbcp_l, gsl_matrix* R, int p, int N);

/**
 *  \brief Test of the symplectic character of a given complex matrix.
 */
void symplecticMatrixTest(const gsl_matrix_complex *M, int INVERSE_TYPE);

/**
 *  \brief Test of the symplectic character of a given real matrix.
 */
void realSymplecticMatrixTest(const gsl_matrix *M, int INVERSE_TYPE);

/**
 *  \brief Test of the Monodromy eigensystem.
 */
void eigenSystemTest(gsl_matrix_complex *Dm, gsl_matrix_complex *S, gsl_matrix_complex **DAT, int M);

/**
 *  \brief Symplectic test of P.
 */
void nfo2SympTest(gsl_odeiv2_driver *d, const double y0[], double t1, QBCP_L *qbcp_l, gsl_matrix* R);

//-------------------------------------------------------------------------------------------------------
// Change of base:
//-------------------------------------------------------------------------------------------------------
/**
 *  \brief Change of base: Dm = Sinv*M*S with M given as a product of matrices
 */
void changeOfBase(gsl_matrix_complex *MMc_estimate, gsl_matrix_complex *Dm, gsl_matrix_complex **DAT, const gsl_matrix_complex *S, int INVERSE_TYPE, int M);

/**
 *  \brief Change of base: Dm = Sinv*M*S with M given as a single matrix
 */
void changeOfBaseSingleMatrix(gsl_matrix_complex *Dm, gsl_matrix_complex *M, gsl_matrix_complex *S, int INVERSE_TYPE);

//-------------------------------------------------------------------------------------------------------
// Wielandt's deflation
//-------------------------------------------------------------------------------------------------------
/**
 *  \brief Application of the Wielandt deflation to the monodromy matrix. Algorithm taken from "Introduction to Numerical Analysis with C programs", A. Mate, 2004.
 *  \param MMc      the monodromy matrix to diagonalize.
 *  \param eigenVu  the unstable eigenvector, obtain via power method applied on MMc = DAT[M]*DAT[M-1]*...*DAT[1].
 *  \param eigenLu  the unstable eigenvalue,  obtain via power method.
 *  \param evecr    matrix of eigenvectors (in columns) of MMc (output).
 *  \param evalr    vectors of eigenvalues of MMc (output).
 *
 *  The algorithm is based on the following results:
 *  Given a eigenvector/value couple \f$ (\lambda,\mathbf{x}) \f$ of \f$ \mathbf{M} \f$, one can build the matrix
 *  \f$ \mathbf{B} = \mathbf{M} - \lambda \mathbf{x} \mathbf{z}^T \f$ where
 *  \f$ \mathbf{z} = \frac{1}{\lambda x_r} \mathbf{M}(r,*)^T \f$ where
 *  \f$ x_r \f$ is the component of \f$ \mathbf{x} \f$ of maximum absolute value, which guarantees \f$ x_r \neq 0 \f$.
 *
 *  Then, the if \f$ \rho \neq 0 \f$ is an eigenvalue of \f$ \mathbf{M} \f$ different from \f$ \lambda \f$, then \f$ \rho \f$ is also an eigenvalue of \f$ \mathbf{B} \f$.
 *  Moreover, if \f$ \mathbf{w} \f$ is the eigenvector of \f$ \mathbf{B} \f$ associated to \f$ \rho \f$ , then the eigenvector \f$ \mathbf{y} \f$ of \f$ \mathbf{M} \f$
 *  associated to \f$ \rho \f$ is given by:
 *
 *  \f$
 *  $$
 *  \mathbf{y} = \frac{1}{\rho} \left(\mathbf{w} + \frac{\lambda}{\rho - \lambda} (\mathbf{z}^T \mathbf{w}) \mathbf{x} \right).
 *  $$
 *  \f$
 **/
void wielandt_deflation(gsl_matrix_complex *MMc,
                        gsl_vector_complex *eigenVu,
                        gsl_complex eigenLu,
                        gsl_matrix_complex *evecr,
                        gsl_vector_complex *evalr);

//-------------------------------------------------------------------------------------------------------
// Multimin (DEPRECATED)
//-------------------------------------------------------------------------------------------------------
//Evaluate ||y - My/l||
double eigenfunk(double p[], double b[], double c, void * params);
//deigenfunk/dxi
void deigenfunk(double p[], double df[], double b[], double c, void *params );
#endif // NFO2_H_INCLUDED