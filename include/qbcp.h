#ifndef QBCP_H_INCLUDED
#define QBCP_H_INCLUDED

/**
 * \file qbcp.h
 * \brief Implementation of various vector fields in the Sun-Earth-Moon QBCP,
 *        as well as additional routines for inner changes of coordinates (Earth-Moon <-> Inertial <-> Sun-(Earth+Moon).
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */


#include <iostream>
#include <time.h>
#include <stdio.h>
#include <complex.h>
#include <gsl_complex_math.h>
#include <fstream>
#include <sstream>

//Custom
#include "ode.h"
#include "qbtbp.h"
#include "diffcorr.h"
#include "define_env.h"
#include "odezero.h"
#include "init.h"
#include "eminsem.h"

extern "C" {
#include "gnuplot_i.h"
}


//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Derivatives for ODE integration
//
//--------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------
// Vector fields
//-------------------------------------------------------
/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates (no variationnal equations).
 **/
int qbfbp_vfn_novar(double t, const double y[], double f[], void *params_void);

/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates.
 **/
int qbfbp_vfn_varnonlin(double t, const double y[], double f[], void *params_void);

/**
 *  \brief Vector field of the QBCP with EM units and EM reference frame.
 **/
int qbfbp_vf(double t, const double y[], double f[], void *params_void);

/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates.
 **/
int qbfbp_vfn_varlin_trans(double t, const double y[], double f[], void *params_void);

/**
 *  \brief Variationnal matrix of the QBCP with EM units and Normalized-Centered coordinates.
 **/
int qbfbp_Dfn_varnonlin(double t, const double y[], double **Df, void *params_void);

//-------------------------------------------------------
// Inner routines for the computation of the vector fields
//-------------------------------------------------------
/**
 *  \brief Build the Variational Equation Matrix Q from the arrays b and alpha
 **/
void set_vareq_matrix(gsl_matrix *Q, double b[], double alpha[]);

/**
 *  \brief Update the vector field of the state in system coordinates (non-normalized)
 **/
int vf_state( const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm );

/**
 *  \brief Update the vector field of the state in Normalized-Centered (NC) coordinates
 **/
int vfn_state(const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm,
              double gamma);

/**
 *  \brief Update the Normalized-Centered variational equation matrix
 **/
int vfn_stm(const double y[], gsl_matrix *Q, double alpha[],
            double ps[], double pe[], double pm[],
            double qps2, double qpe2, double qpm2,
            double ms, double me, double mm,
            double gamma);

/**
 *  \brief Update the Normalized-Centered variational equation matrix - case of a double **Df, instead of a gsl_matrix.
 **/
int vfn_stm(const double y[], double **Df, double alpha[],
            double ps[], double pe[], double pm[],
            double qps2, double qpe2, double qpm2,
            double ms, double me, double mm,
            double gamma);


/**
 *  \brief Update the variational equation matrix in system coordinates (non-normalized)
 **/
int vf_stm(const double y[], gsl_matrix *Q, double alpha[],
           double ps[], double pe[], double pm[],
           double qps2, double qpe2, double qpm2,
           double ms, double me, double mm);

/**
 *  \brief Update the Normalized-Centered linearized variational equation matrix.
 **/
int vfn_stm_lin_trans(const double y[], gsl_matrix *Q, double alpha[],
                     double ps[], double pe[], double pm[],
                     double ms, double me, double mm,
                     double gamma);

/**
 *  \brief Computes the second-order derivatives of the linearized potential of one given primary
 **/
double Uijlin(double pc[], double qpc2, double mc, double factor, int i, int j);

/**
 *  \brief Computes the second-order derivatives of the potential of one given primary
 **/
double Uij(const double y[], double pc[], double qpc2, double mc, double factor, int i, int j);

//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Derivatives for ODE integration in continuation procedures
//
//--------------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Vector field for continuation procedure. The result is f = (1.0 - epsilon)*f1 + espilon*f2. Normalized-Centered coordinates. Variational equations included.
 **/
int qbfbp_vfn_cont(double t, const double y[], double f[], void *params_void);

/**
 *  \brief Vector field for continuation procedure. The result is f = (1.0 - epsilon)*f1 + espilon*f2. System coordinates (non-normalized).
 **/
int qbfbp_vf_cont(double t, const double y[], double f[], void *params_void);

//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Hamiltonians
//
//--------------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Hamiltonian of the QBCP with EM units and EM coordinates
 **/
double qbfbp_H(double t, const double y[], void *params_void);

/**
 *  \brief Hamiltonian of the QBCP with EM units and Normalized-Centered coordinates
 **/
double qbfbp_Hn(double t, const double y[], void *params_void);

//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Floquet analysis
//
//--------------------------------------------------------------------------------------------------------------------------------------------
//
/**
 *  \brief Same as the routine qbfbp_vfn + Integration of the matrix Pbar that appears in the Floquet analysis of the hamiltonian of order 2 in the neighborhood of L1,2.
 *         (see nfo2.h)
 **/
int qbfbp_vfn_varlin_trans_P(double t, const double y[], double f[], void *params_void);

/**
 *  \brief QBCP vector field in the form of the matrix Q1(t), Q2(t) and Q3(t) defined by:
 *              x' =  Q2^T*x + 2Q3*Y
 *              y' = -2Q1 *x -  Q2*y
 *          Used in the Floquet analysis (see nfo2.h)
 **/
int qbfbp_Q(double t, const double y[], gsl_matrix *Q1, gsl_matrix *Q2, gsl_matrix *Q3, void *params_void);

//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Change of coordinates: SEM <-> IN <-> EM
//
//--------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// COC: NC <--> EM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from Normalized-Centered coordinates to Earth-Moon coordinates
 **/
void NCtoEM(double t, const double yNC[], double yEM[], FBPL *qbp);

/**
 *  \brief COC: from Earth-Moon coordinates to Normalized-Centered coordinates
 **/
void EMtoNC(double t, const double yEM[], double yNC[], FBPL *qbp);

//-----------------------------------------------------------------------------
// COC: NC <--> SYS
//-----------------------------------------------------------------------------
/**
 *  \brief COC: Normalized-Centered coordinates to system coordinates. Use in priority instead of NCtoEM or NCtoSEM.
 **/
void NCtoSYS(double t, const double yNC[], double yEM[], FBPL *qbp);

/**
 *  \brief COC: from system coordinates to Normalized-Centered coordinates
 **/
void SYStoNC(double t, const double yEM[], double yNC[], FBPL *qbp);

//-----------------------------------------------------------------------------
// COC: NC <--> SEM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from Sun-Earth-Moon coordinates to Normalized-Centered coordinates
 **/
void SEMtoNC(double t, const double ySEM[], double yNC[], FBPL *qbp);

/**
 *  \brief COC: from Normalized-Centered coordinates to Sun-Earth-Moon coordinates
 **/
void NCtoSEM(double t, const double yNC[], double ySEM[], FBPL *qbp);

//-----------------------------------------------------------------------------
// COC: tests & plots
//-----------------------------------------------------------------------------
/**
 *  \brief Plot the Primaries' Three-Body motion of the QBCP.
 **/
void QBTBP_IN();

/**
 *  \brief Computes a solution in NC EM coordinates and plot it in NC SEM coordinates in a temporary gnuplot window.
 **/
int odePlot_EMtoSEM(const double y[],
                      double t1,
                      FBPL &fbpl,
                      gsl_odeiv2_driver *d,
                      gnuplot_ctrl  *h1,
                      int Npoints,
                      char const *title,
                      char const *ls,
                      char const *lt,
                      char const *lw,
                      int lc,
                      string folder);

/**
 *  \brief Computes a solution in NC SEM coordinates and plot it in NC EM coordinates in a temporary gnuplot window.
 **/
int odePlot_SEMtoEM(const double y[],
                      double t1,
                      gsl_odeiv2_driver *d,
                      gnuplot_ctrl  *h1,
                      int Npoints,
                      char const *title,
                      char const *ls,
                      char const *lt,
                      char const *lw,
                      int lc,
                      string folder);

/**
 *  \brief Test of the dynamics: SEM to EM system.
 **/
void dynTest_SEMtoEM();

/**
 *  \brief Test of the dynamics: EM to SEM system.
 **/
void dynTest_EMtoSEM();




#endif // QBCP_H_INCLUDED
