#ifndef QBTBP_OFS_H_INCLUDED
#define QBTBP_OFS_H_INCLUDED

/**
 * \file qbtbp.h
 * \brief Routines of the Quasi-Bicircular Three-Body Problem (QBTBP).
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <complex.h>
//#include <fstream>

//GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>

#include <gsl/gsl_sf_bessel.h>

//Custom
#include "ofts.h"
#include "ode.h"
#include "gslc.h"
#include "define_env.h"
#include "config.h"
#include "eminsem.h"

extern "C" {
#include "nrutil.h"
}


using namespace std;
//-----------------------------------------------------------------------------
// Main routine
//-----------------------------------------------------------------------------
/**
 *  \brief Main routine to compute the QBTBP in Ofs format. The iteratives equation of the QBTBP are solved using qbtbp_ofs.
 *   If necessary, the result is tested using qbtbp_test.
 */
void qbtbp(int li_EM, int li_SEM, int isTestOn, int fwrk);

/**
 *  \brief Main routine to compute the Bicircular Three-Body Problem in Ofs format.
 */
void bcp(int li_EM, int li_SEM, int fwrk);

//-----------------------------------------------------------------------------
// Computing the QBTBP
//-----------------------------------------------------------------------------

/**
 *  \brief Computes the QBTBP in Ofs format.
 *  \param zr_ofts: a reference to reduced internal motion \f$ z \f$ so that: \f$ z_r = \sum \limits_{j} b_j e^{ijnt} \f$ and \f$ z = e^{it} z_r \f$ in Earth-Moon units.
 *  \param Zr_ofts: a reference to reduced external  motion \f$ Z \f$ so that: \f$ Z_r = \sum \limits_{j} c_j e^{ijnt} \f$ and \f$ Z = a_s e^{in_st} Z_r \f$ in Earth-Moon units.
 *  \param qbcp_l: a reference to the current qbcp, focused on a given point \f$ L_{1,2} \f$.
 */
void qbtbp_ofs (Ofts< Ofsd > &zr_ofts, Ofts< Ofsd > &Zr_ofts, QBCP_L& qbcp_l, int fwrk);


/**
 *  \brief One step of the recurrence scheme in qbtbp_ofs. See comments in src/h for the signification of the parameters.
 */
void qbtbp_ofs_recurrence(/**reduced internal motion */
                          Ofts< Ofsd > &zt,
                          /**reduced external motion */
                          Ofts< Ofsd > &Zt,
                          /**\f$ \bar{z}_t \f$ */
                          Ofts< Ofsd > &z1,
                          /**\f$ \bar{z}_t^{-3/2} \f$*/
                          Ofts< Ofsd > &z2,
                          /**\f$ z_t^{-1/2} \f$ */
                          Ofts< Ofsd > &z3,
                          /**\f$ z_2 z_3 = z_t^{-1/2}*\bar{z}_t^{-3/2} \f$ */
                          Ofts< Ofsd > &z4,
                          /**\f$ \bar{z}_t^{-1/2} \f$ */
                          Ofts< Ofsd > &z5,
                          /**\f$ z_3 z_5 = (z_t \bar{z}_t)^{-1/2} \f$ */
                          Ofts< Ofsd > &z6,
                          /**\f$ \mu e^{i \theta} z_t \f$ */
                          Ofts< Ofsd > &a1,
                          /** \f$ \epsilon a_1 \f$ */
                          Ofts< Ofsd > &a2,
                          /** \f$ Z_t- \mu e^{i \theta} \epsilon z_t \f$ */
                          Ofts< Ofsd > &a3,
                          /** \f$ a_3^{-1/2}. \f$ */
                          Ofts< Ofsd > &a4,
                          /** \f$ \bar{a}_3 \f$ */
                          Ofts< Ofsd > &a5,
                          /** \f$ a_5^{-3/2} \f$ */
                          Ofts< Ofsd > &a6,
                          /** \f$ a_4 a_6 \f$ */
                          Ofts< Ofsd > &a7,
                          /**\f$  (1-\mu)*e^{i \theta}*z_t \f$ */
                          Ofts< Ofsd > &b1,
                          /** \f$ \epsilon b_1 \f$  */
                          Ofts< Ofsd > &b2,
                          /** \f$ Z_t- (1-\mu) e^{i \theta} \epsilon z_t \f$ */
                          Ofts< Ofsd > &b3,
                          /** \f$ b_3^{-1/2}. \f$ */
                          Ofts< Ofsd > &b4,
                          /**\f$ \bar{b}_3 \f$ */
                          Ofts< Ofsd > &b5,
                          /**\f$ b_5^{-3/2} \f$ */
                          Ofts< Ofsd > &b6,
                          /**\f$ b_4 b_6 \f$ */
                          Ofts< Ofsd > &b7,
                          /**\f$ -\frac{m_s}{a_s^2} e^{-i \theta} b_7 + \frac{m_s}{a_s^2} e^{-i \theta} a_7 - z_4 \f$ */
                          Ofts< Ofsd > &Pm,
                          /**\f$ -n_s^2 \mu b_7 - n_s^2 (1-\mu) a_7 \f$ */
                          Ofts< Ofsd > &Qm,
                          /** \f$ \epsilon \f$ = Ofts with just 1.0 at order 1*/
                          Ofts< Ofsd > &epsilon,
                          ///Pfm = Pm(j) at order n
                          Ofsd  &Pfm,
                          ///Qfm = Qm(j) at order n
                          Ofsd  &Qfm,
                          ///ufm = zt(j) at order n
                          Ofsd  &ufm,
                          ///vfm = Zt(j) at order n
                          Ofsd  &vfm,
                          /** \f$ \sigma_1 = e^{-i\theta} \f$ */
                          Ofsd &sigma1,
                          /** \f$ \sigma_2 = e^{+i\theta} \f$ */
                          Ofsd &sigma2,
                          ///order
                          int m,
                          ///order of the Fourier expansions
                          int nf,
                          ///current QBCP_L
                          QBCP_L &qbcp_l);


/**
 *  \brief Manipulation and Storage in txt files of OFS objects in qbtbp_ofs. See comments in src/h for the signification of the parameters.
 */
void qbtbp_ofs_storage(/**reduced internal motion */
                          Ofts< Ofsd > &zt,
                          /**reduced external motion */
                          Ofts< Ofsd > &Zt,
                          /**\f$ \bar{z}_t \f$ */
                          Ofts< Ofsd > &z1,
                          /**\f$ \bar{z}_t^{-3/2} \f$*/
                          Ofts< Ofsd > &z2,
                          /**\f$ z_t^{-1/2} \f$ */
                          Ofts< Ofsd > &z3,
                          /**\f$ z_2 z_3 = z_t^{-1/2}*\bar{z}_t^{-3/2} \f$ */
                          Ofts< Ofsd > &z4,
                          /**\f$ \bar{z}_t^{-1/2} \f$ */
                          Ofts< Ofsd > &z5,
                          /**\f$ z_3 z_5 = (z_t \bar{z}_t)^{-1/2} \f$ */
                          Ofts< Ofsd > &z6,
                          /**\f$ \mu e^{i \theta} z_t \f$ */
                          Ofts< Ofsd > &a1,
                          /** \f$ \epsilon a_1 \f$ */
                          Ofts< Ofsd > &a2,
                          /** \f$ Z_t- \mu e^{i \theta} \epsilon z_t \f$ */
                          Ofts< Ofsd > &a3,
                          /** \f$ a_3^{-1/2}. \f$ */
                          Ofts< Ofsd > &a4,
                          /** \f$ \bar{a}_3 \f$ */
                          Ofts< Ofsd > &a5,
                          /** \f$ a_5^{-3/2} \f$ */
                          Ofts< Ofsd > &a6,
                          /** \f$ a_4 a_6 \f$ */
                          Ofts< Ofsd > &a7,
                          /**\f$  (1-\mu)*e^{i \theta}*z_t \f$ */
                          Ofts< Ofsd > &b1,
                          /** \f$ \epsilon b_1 \f$  */
                          Ofts< Ofsd > &b2,
                          /** \f$ Z_t- (1-\mu) e^{i \theta} \epsilon z_t \f$ */
                          Ofts< Ofsd > &b3,
                          /** \f$ b_3^{-1/2}. \f$ */
                          Ofts< Ofsd > &b4,
                          /**\f$ \bar{b}_3 \f$ */
                          Ofts< Ofsd > &b5,
                          /**\f$ b_5^{-3/2} \f$ */
                          Ofts< Ofsd > &b6,
                          /**\f$ b_4 b_6 \f$ */
                          Ofts< Ofsd > &b7,
                          /**\f$ -\frac{m_s}{a_s^2} e^{-i \theta} b_7 + \frac{m_s}{a_s^2} e^{-i \theta} a_7 - z_4 \f$ */
                          Ofts< Ofsd > &Pm,
                          /**\f$ -n_s^2 \mu b_7 - n_s^2 (1-\mu) a_7 \f$ */
                          Ofts< Ofsd > &Qm,
                          ///Pfm = Pm(j) at order n
                          Ofsd  &Pfm,
                          ///Qfm = Qm(j) at order n
                          Ofsd  &Qfm,
                          ///ufm = zt(j) at order n
                          Ofsd  &ufm,
                          ///vfm = Zt(j) at order n
                          Ofsd  &vfm,
                          /** \f$ \sigma_1 = e^{-i\theta} \f$ */
                          Ofsd &sigma1,
                          /** \f$ \sigma_2 = e^{+i\theta} \f$ */
                          Ofsd &sigma2,
                          ///order of the Fourier expansions
                          int nf,
                          ///current QBCP_L
                          QBCP_L &qbcp_l);




/**
 *  \brief FFT and Storage in txt files of OFS objects alpha_i in qbtbp_ofs. See comments in src/h for the signification of the parameters.
 */
void qbtbp_ofs_fft_alpha( Ofts< Ofsd > &zt,     //zt = normalized Earth-Moon motion
                          Ofts< Ofsd > &Zt,     //Zt = normalized Sun-(Earth+Moon) motion
                          int nf,                      //order of the Fourier expansions
                          QBCP_L& qbcp_l);            //QBCP
/**
 *  \brief FFT and Storage in txt files of OFS objects delta_i in qbtbp_ofs.
 */
void qbtbp_ofs_fft_delta(Ofts< Ofsd > &zt,    //zt = normalized Earth-Moon motion
                        Ofts< Ofsd > &Zt,     //Zt = normalized Sun-(Earth+Moon) motion
                        int nf,                      //order of the Fourier expansions
                        QBCP_L &qbcp_l);            //QBCP


//-----------------------------------------------------------------------------
// Integrating the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Derivatives of the QBTBP. To plug into GSL integrator.
 */
int qbtbp_derivatives(double t, const double y[], double f[], void *params);

//-----------------------------------------------------------------------------
// Testing the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Test function to compare the analytical solution of the QBTBP to the numerical integration of the equations of motion.
 */
void qbtbp_test(double t1, Ofsc &bjc, Ofsc &cjc, OdeStruct ode_s, QBCP_L &qbcp_l);

/**
 *  \brief Test function to compare the analytical solutions found in IN, EM, SE, and SEM coordinates
 */
void qbtbp_test_IN_EM_SEM(double t1, Ofsc &bjc, Ofsc &cjc);

/**
 *  \brief Comparison of the QBTBP computed via FFT or via OFS computations
 */
void qbtbp_test_FFT_vs_OFS(Ofsc &bjc,    //zt(t)
                      Ofsc &cjc,    //Zt(t)
                      int nf,                      //order of the Fourier expansions
                      int N,                       //Number of points
                      int type,                    //Type of reference
                      OdeStruct ode_s,             //ode structure
                      QBCP_L& qbcp_l);           //QBCP
//-----------------------------------------------------------------------------
// Evaluate the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluate z(t), with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evz(Ofsc& zr, double t, double n, double ni, double ai);

/**
 *  \brief Evaluate dz(t)/dt, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evzdot(Ofsc& zr, Ofsc& ztdot, double t, double n, double ni, double ai);

/**
 *  \brief Evaluate d2z(t)/dt2, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evzddot(Ofsc& zr, Ofsc& ztdot, Ofsc& ztddot, double t, double n, double ni, double ai);


//-----------------------------------------------------------------------------
//Evaluation of parts or all the alpha/beta routines at a given time
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluation of Fourier series given as an array of coefficients, at a given time t.
 */
void evaluateCoef(double *alpha, double t, double omega, int order, double *params, int number);

/**
 *  \brief Evaluation of the time derivatives Fourier series given as an array of coefficients, at a given time t.
 */
void evaluateCoefDerivatives(double *alpha, double t, double omega, int order, double *params, int number);

//-----------------------------------------------------------------------------
//Evalution of individual alpha/beta
//or derivative of alpha/beta of a given type (Odd or Even)
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) cos(k \omega t)  \f$.
 */
double evaluateEven(double t, double omega, int order, double *coef, double *cR);

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N - k \omega coef(k) sin(k \omega t)  \f$.
 */
double evaluateEvenDerivative(double t, double omega,  int order, double *coef, double *sR);

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) sin(k \omega t)  \f$.
 */
double evaluateOdd(double t, double omega,  int order, double *coef, double *sR);

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N  k \omega coef(k) cos(k \omega t)  \f$.
 */
double evaluateOddDerivative(double t, double omega,  int order, double *coef, double *cR);

//-----------------------------------------------------------------------------
// Backup void qbtbp_ofs and qbtbp_ots
//-----------------------------------------------------------------------------
/**
 *  \brief Solution of the qbtbp in Ots format (obsolete)
 */
void qbtbp_ots(Ofts< Ots<double> > &ofts_z, Ofts< Ots<double> > &ofts_Z, double n, double ms, double as, double ns);

/**
 * \brief Calculates the Fourier transform of a set of n real-valued data points. Note currently used.
 *
 * Replaces this data (which
 * is stored in array data[0..n-1] ) by the positive frequency half of its complex Fourier transform.
 * The real-valued first and last components of the complex transform are returned as elements
 * data[1] and data[2] , respectively. n must be a power of 2. This routine also calculates the
 * inverse transform of a complex data array if it is the transform of real data. (Result in this case
 * must be multiplied by 2/n .)
 **/
void realft(double data[], unsigned long n, int isign);

#endif // QBTBP_OFS_H_INCLUDED
