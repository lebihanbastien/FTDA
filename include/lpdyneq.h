#ifndef LPDYNEQ_H_INCLUDED
#define LPDYNEQ_H_INCLUDED

#include "init.h"
#include "eminsem.h"
extern "C" {
#include "gnuplot_i.h"
}

//==============================================================================
//Dynamical equivalents to the Libration points
//==============================================================================
/**
 *  \brief Computation of the dynamical equivalents of the libration points. Test function.
 *
 *         - This function makes use of the subroutine lpdyneq that is the true heart of the computation. In particular, it contains the initialization of the first guess (the geometrical positions of the CRTBP libration points) + the differential corrector.
 *         - After lpdyneq, the resulting initial conditions are integrated and plotted on a full orbit.
 *         - Finally, a test of the periodicity of the orbit is performed, via the computation of the error |y(0) - y(T)|.
 *         - If isStored is true, the results (x(t), y(t)) in synodical coordinates are stored in txt files of the form: "./plot/QBCP/DYNEQ/DYNEQ_QBCP_EM_L1.txt".
 **/
void compute_dyn_eq_lib_point(FBPL &fbpl, int isStored);

/**
 *  \brief Main routine for the computation of the dynamical equivalent to the Libration points.
 *         The results are given in the form of the initial conditions (42 dimensions) updated in the array y0[].
 *
 *          €€TODO: BCP case is a work in progress. works only for EML1, in NC coordinates.
 **/
void lpdyneq(gsl_odeiv2_driver *d, double y0[]);

//==============================================================================
// Dynamical equivalents via continuation procedures
//==============================================================================
/**
 *  \brief Computation of lpdyneq with continuation between two models. WORK IN PROGRESS.
 **/
void lpdyneq_cont(gsl_odeiv2_driver *d, gsl_odeiv2_control * loose_control, gsl_odeiv2_control * hard_control, double y0[], gnuplot_ctrl *h1);

/**
 *  \brief Computation of lpdyneq with continuation between two models. WORK IN PROGRESS.
 **/
void lpdyneq_cont_2(gsl_odeiv2_driver *d, gsl_odeiv2_control * loose_control, gsl_odeiv2_control * hard_control, double y0[]);

/**
 *  \brief Continuation routine for the dynamical equivalents of the libration points.
 */
void continuation_dyn_eq_lib_point(FBPL &fbpl, int from_model, int to_model);

/**
 *  \brief Continuation routine for other orbits than the libration points.
 */
void continuation_res_orbit(FBPL &fbpl, int from_model, int to_model);

//------------------------------------------------------------------------------------------------------------
// Periodicity Condition
//------------------------------------------------------------------------------------------------------------
/**
 *  \brief Test the periodicity condition y(0) = y(t1), with initial condition y, for the vector field contained in the driver d.
 *
 *  Note that the integer N is the number of variables associated to the driver d, and should be also the number of variables in y. However, the periodicity condition is tested only on
 *  NvarTest variables (a usual example is Nvar = 42 but NvarTest = 6).
 **/
int periodicity_condition(const double y[], int Nvar, int NvarTest, double t1, gsl_odeiv2_driver *d, int isNorm);


//==============================================================================
// Subroutines
//==============================================================================
/**
 *   \brief Computes M+1 patch points, indexed between 0 and M,
 *          along the T-periodic orbit starting at y0.
 *          Each STM is initialized with the identity matrix.
 **/
void lpdyneq_patch_points(const double *y0, gsl_odeiv2_driver *d, double T, double **ym, double *tm, int M);

/**
 *  \brief Inverse a 2x2 matrix. Used in lpdyneq_cont.
 **/
void invMat22(double A[2][2], double invA[2][2]);

/**
 *  \brief Arclength vector field. For all i in [[1,8]], fi = Ai/sum(Aj^2).
 **/
void arclengthvf(double *f, double *Aj);

#endif // LPDYNEQ_H_INCLUDED
